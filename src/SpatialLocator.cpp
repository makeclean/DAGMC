#include "moab/SpatialLocator.hpp"
#include "moab/Interface.hpp"
#include "moab/ElemEvaluator.hpp"
#include "moab/AdaptiveKDTree.hpp"

// include ScdInterface for box partitioning
#include "moab/ScdInterface.hpp"

#ifdef USE_MPI
#include "moab/ParallelComm.hpp"
#endif

bool debug = true;

namespace moab 
{

    SpatialLocator::SpatialLocator(Interface *impl, Range &elems, Tree *tree, ElemEvaluator *eval) 
            : mbImpl(impl), myElems(elems), myDim(-1), myTree(tree), elemEval(eval), iCreatedTree(false)
    {
      if (!myTree) {
        myTree = new AdaptiveKDTree(impl);
        iCreatedTree = true;
      }
      if (!elems.empty()) {
        myDim = mbImpl->dimension_from_handle(*elems.rbegin());
        ErrorCode rval = myTree->build_tree(myElems);
        if (MB_SUCCESS != rval) throw rval;
      }
    }

    ErrorCode SpatialLocator::add_elems(Range &elems) 
    {
      if (elems.empty() ||
          mbImpl->dimension_from_handle(*elems.begin()) != mbImpl->dimension_from_handle(*elems.rbegin()))
        return MB_FAILURE;
  
      myDim = mbImpl->dimension_from_handle(*elems.begin());
      myElems = elems;
      return MB_SUCCESS;
    }
    
#ifdef USE_MPI
    ErrorCode SpatialLocator::initialize_global_box(ParallelComm *pc) 
    {
      if (!pc) return MB_FAILURE;
      
      BoundBox box, gbox;
      ErrorCode rval = myTree->get_bounding_box(box);
      
        //step 2
        // get the global bounding box
      double sendbuffer[6];
      double rcvbuffer[6];

      box.get(sendbuffer); //fill sendbuffer with local box, max values in [0:2] min values in [3:5]
      sendbuffer[3] *= -1;
      sendbuffer[4] *= -1; //negate Xmin,Ymin,Zmin to get their minimum using MPI_MAX
      sendbuffer[5] *= -1; //to avoid calling MPI_Allreduce again with MPI_MIN

      int mpi_err = MPI_Allreduce(sendbuffer, rcvbuffer, 6, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      if (MPI_SUCCESS != mpi_err)	return MB_FAILURE;

      rcvbuffer[3] *= -1;
      rcvbuffer[4] *= -1;  //negate Xmin,Ymin,Zmin again to get original values
      rcvbuffer[5] *= -1;

      globalBox.update_max(&rcvbuffer[0]); //saving values in globalBox
      globalBox.update_min(&rcvbuffer[3]);

        // compute the alternate decomposition; use ScdInterface::compute_partition_sqijk for this
      ScdParData spd;
      spd.partMethod = ScdParData::SQIJK;
      spd.gPeriodic[0] = spd.gPeriodic[1] = spd.gPeriodic[2] = 0;
      double lg = log10((box.bMax - box.bMin).length());
      double mfactor = pow(10.0, 6 - lg);
      int ldims[3], lper[3];
      double dgijk[6];
      box.get(dgijk);
      for (int i = 0; i < 6; i++) spd.gDims[i] = dgijk[i] * mfactor;
      rval = ScdInterface::compute_partition(pc->size(), pc->rank(), spd,
                                             ldims, lper, regNums);
      if (MB_SUCCESS != rval) return rval;
        // all we're really interested in is regNums[i], #procs in each direction
      
      for (int i = 0; i < 3; i++)
        regDeltaXYZ[i] = (globalBox.bMax[i] - globalBox.bMin[i])/double(regNums[i]); //size of each region

      return MB_SUCCESS;
    }

//this function sets up the TupleList TLreg_o containing the registration messages
// and sends it
    ErrorCode SpatialLocator::register_with_forwarders(ParallelComm *pc, TupleList &TLreg_o)
    {

        //step 4
        //set up TLreg_o
      TLreg_o.initialize(1,0,0,6,0);
        // TLreg_o (int destProc, real Xmin, Ymin, Zmin, Xmax, Ymax, Zmax)

      int dest;
      double boxtosend[6];

      BoundBox box;
      ErrorCode result = myTree->get_bounding_box(box);
      if (result != MB_SUCCESS)
        return result;

      box.get(boxtosend);

        //iterate over all regions overlapping with my bounding box using the computerd corner IDs
      for (int k = cornerIJK[2]; k <= cornerIJK[5]; k++) {
        for (int j = cornerIJK[1]; j <= cornerIJK[4]; j++) {
          for (int i = cornerIJK[0]; i <= cornerIJK[3]; i++) {
            dest = k * regNums[0]*regNums[1] + j * regNums[0] + i;
            TLreg_o.push_back(&dest, NULL, NULL, boxtosend);
          }
        }
      }
	
        //step 5
        //send TLreg_o, receive TLrequests_i
      if (pc) pc->proc_config().crystal_router()->gs_transfer(1, TLreg_o, 0);

        //step 6
        //Read registration requests from TLreg_o and add to list of procs to forward to
        //get number of tuples sent to me

        //read tuples and fill processor list;
      int NN = TLreg_o.get_n();
      for (int i=0; i < NN; i++)
          //TLreg_o is now TLrequests_i
        srcProcBoxes[TLreg_o.vi_rd[i]] = BoundBox(TLreg_o.vr_rd+6*i);

      return MB_SUCCESS;
    }

    ErrorCode SpatialLocator::par_locate_points(ParallelComm */*pc*/,
                                                Range &/*vertices*/,
                                                const double /*rel_iter_tol*/, const double /*abs_iter_tol*/,
                                                const double /*inside_tol*/)
    {
      return MB_UNSUPPORTED_OPERATION;
    }

    bool is_initialized(int i) {return (i == -1);}
      
    ErrorCode SpatialLocator::par_locate_points(ParallelComm *pc,
                                                const double *pos, int num_points,
                                                const double rel_iter_tol, const double abs_iter_tol,
                                                const double inside_tol)
    {
      ErrorCode rval;
        //TUpleList used for communication 
      TupleList TLreg_o;   //TLregister_outbound
      TupleList TLquery_o; //TLquery_outbound
      TupleList TLforward_o; //TLforward_outbound
      TupleList TLsearch_results_o; //TLsearch_results_outbound

      
        // steps 1-2 - initialize the alternative decomposition box from global box
      rval = initialize_global_box(pc);
      
        //step 3 - compute the IDs of the regions which contain each corner of local box
      rval = compute_corner_ijks();
      if (rval != MB_SUCCESS) return rval;

        //steps 4-6 - set up TLreg_o, gs_transfer, gather registrations
      rval = register_with_forwarders(pc, TLreg_o);
      if (rval != MB_SUCCESS) return rval;

        // actual parallel point location using intermediate partition

        // target_pts: TL(to_proc, tgt_index, x, y, z): tuples sent to source mesh procs representing pts to be located
        // source_pts: TL(from_proc, tgt_index, src_index): results of source mesh proc point location, ready to send
        //             back to tgt procs; src_index of -1 indicates point not located (arguably not useful...)

      unsigned int my_rank = (pc? pc->proc_config().proc_rank() : 0);

        //TLquery_o: Tuples sent to forwarder proc 
        //TL (toProc, OriginalSourceProc, targetIndex, X,Y,Z)

        //TLforw_req_i: Tuples to forward to corresponding procs (forwarding requests)
        //TL (sourceProc, OriginalSourceProc, targetIndex, X,Y,Z)

      TLquery_o.initialize(3,0,0,3,0);

      int iargs[3];

      for (int pnt=0; pnt < 3*num_points; pnt+=3)
      {
        int forw_id = proc_from_point(pos+pnt); //get ID of proc resonsible of the region the proc is in

        iargs[0] = forw_id; 	//toProc
        iargs[1] = my_rank; 	//originalSourceProc
        iargs[2] = pnt/3;    	//targetIndex 	

        TLquery_o.push_back(iargs, NULL, NULL, const_cast<double*>(pos+pnt));
      }

        //send point search queries to forwarders
      if (pc)
        pc->proc_config().crystal_router()->gs_transfer(1, TLquery_o, 0);

        //now read forwarding requests and forward to corresponding procs
        //TLquery_o is now TLforw_req_i

        //TLforward_o: query messages forwarded to corresponding procs
        //TL (toProc, OriginalSourceProc, targetIndex, X,Y,Z)

      TLforward_o.initialize(3,0,0,3,0);

      int NN = TLquery_o.get_n();

      for (int i=0; i < NN; i++) {
        iargs[1] = TLquery_o.vi_rd[3*i+1];	//get OriginalSourceProc
        iargs[2] = TLquery_o.vi_rd[3*i+2];	//targetIndex
        CartVect tmp_pnt(TLquery_o.vr_rd+3*i);

          //compare coordinates to list of bounding boxes
        for (std::map<int, BoundBox>::iterator mit = srcProcBoxes.begin(); mit != srcProcBoxes.end(); mit++) {
          if ((*mit).second.contains_point(tmp_pnt.array(), abs_iter_tol)) {
            iargs[0] = (*mit).first;
            TLforward_o.push_back(iargs, NULL, NULL, tmp_pnt.array());
          }
        }

      }

      if (pc)
        pc->proc_config().crystal_router()->gs_transfer(1, TLforward_o, 0);

        //step 12
        //now read Point Search requests
        //TLforward_o is now TLsearch_req_i
        //TLsearch_req_i: (sourceProc, OriginalSourceProc, targetIndex, X,Y,Z)
							  
      NN = TLforward_o.get_n();

        //TLsearch_results_o
        //TL: (OriginalSourceProc, targetIndex, sourceIndex, U,V,W);
      TLsearch_results_o.initialize(3,0,0,3,0);

        //step 13 is done in test_local_box

      std::vector<double> params(3*NN);
      std::vector<int> is_inside(NN, 0);
      std::vector<EntityHandle> ents(NN, 0);
      
      rval = locate_points(TLforward_o.vr_rd, TLforward_o.get_n(), 
                           &ents[0], &params[0], &is_inside[0], 
                           rel_iter_tol, abs_iter_tol, inside_tol);
      if (MB_SUCCESS != rval)
        return rval;
      
      for (int i = 0; i < NN; i++) {
        if (is_inside[i]) {
          iargs[0] = TLforward_o.vi_rd[3*i+1];
          iargs[1] = TLforward_o.vi_rd[3*i+2];
          iargs[2] = locTable.get_n();
          TLsearch_results_o.push_back(iargs, NULL, NULL, NULL);
          locTable.push_back(const_cast<int*>(&TLforward_o.vi_rd[3*i+1]), NULL, &ents[i], &params[3*i]);
        }
      }

        //step 14: send TLsearch_results_o and receive TLloc_i
      if (pc)
        pc->proc_config().crystal_router()->gs_transfer(1, TLsearch_results_o, 0);


        // store proc/index tuples in parLocTable
      parLocTable.initialize(2, 0, 0, 0, num_points);
      parLocTable.enableWriteAccess();
      std::fill(parLocTable.vi_wr, parLocTable.vi_wr + 2*num_points, -1);
      
      for (unsigned int i = 0; i < TLsearch_results_o.get_n(); i++) {
        int idx = TLsearch_results_o.vi_rd[3*i+1];
        parLocTable.vi_wr[2*idx] = TLsearch_results_o.vi_rd[3*i];
        parLocTable.vi_wr[2*idx+1] = TLsearch_results_o.vi_rd[3*i+2];
      }

      if (debug) {
        int num_found = num_points - 0.5 * 
            std::count_if(parLocTable.vi_wr, parLocTable.vi_wr + 2*num_points, is_initialized);
        std::cout << "Points found = " << num_found << " out of " << num_points << "." << std::endl;
      }
      
      return MB_SUCCESS;
    }

#endif

    ErrorCode SpatialLocator::locate_points(Range &verts,
                                            const double rel_iter_tol, const double abs_iter_tol, 
                                            const double inside_tol) 
    {
      assert(!verts.empty() && mbImpl->type_from_handle(*verts.rbegin()) == MBVERTEX);
      std::vector<double> pos(3*verts.size());
      ErrorCode rval = mbImpl->get_coords(verts, &pos[0]);
      if (MB_SUCCESS != rval) return rval;
      rval = locate_points(&pos[0], verts.size(), rel_iter_tol, abs_iter_tol, inside_tol);
      if (MB_SUCCESS != rval) return rval;
      
      return MB_SUCCESS;
    }
    
    ErrorCode SpatialLocator::locate_points(const double *pos, int num_points,
                                            const double rel_iter_tol, const double abs_iter_tol, 
                                            const double inside_tol) 
    {
        // initialize to tuple structure (p_ui, hs_ul, r[3]_d) (see header comments for locTable)
      locTable.initialize(1, 0, 1, 3, num_points);
      locTable.enableWriteAccess();

        // pass storage directly into locate_points, since we know those arrays are contiguous
      ErrorCode rval = locate_points(pos, num_points, locTable.vul_wr, locTable.vr_wr, NULL, rel_iter_tol, abs_iter_tol,
                                     inside_tol);
      std::fill(locTable.vi_wr, locTable.vi_wr+num_points, 0);
      locTable.set_n(num_points);
      if (MB_SUCCESS != rval) return rval;
      
      return MB_SUCCESS;
    }
      
    ErrorCode SpatialLocator::locate_points(Range &verts,
                                            EntityHandle *ents, double *params, int *is_inside,
                                            const double rel_iter_tol, const double abs_iter_tol, 
                                            const double inside_tol)
    {
      assert(!verts.empty() && mbImpl->type_from_handle(*verts.rbegin()) == MBVERTEX);
      std::vector<double> pos(3*verts.size());
      ErrorCode rval = mbImpl->get_coords(verts, &pos[0]);
      if (MB_SUCCESS != rval) return rval;
      return locate_points(&pos[0], verts.size(), ents, params, is_inside, rel_iter_tol, abs_iter_tol, inside_tol);
    }

    ErrorCode SpatialLocator::locate_points(const double *pos, int num_points,
                                            EntityHandle *ents, double *params, int *is_inside,
                                            const double rel_iter_tol, const double abs_iter_tol, 
                                            const double inside_tol)
    {
      double tmp_abs_iter_tol = abs_iter_tol;
      if (rel_iter_tol && !tmp_abs_iter_tol) {
          // relative epsilon given, translate to absolute epsilon using box dimensions
        BoundBox box;
        myTree->get_bounding_box(box);
        tmp_abs_iter_tol = rel_iter_tol * box.diagonal_length();
      }
  
      EntityHandle closest_leaf;
      std::vector<double> dists;
      std::vector<EntityHandle> leaves;
      ErrorCode rval = MB_SUCCESS;

      for (int i = 0; i < num_points; i++) {
        int i3 = 3*i;
        ents[i] = 0;
        if (tmp_abs_iter_tol) {
          rval = myTree->distance_search(pos+i3, tmp_abs_iter_tol, leaves, tmp_abs_iter_tol, inside_tol, &dists);
          if (MB_SUCCESS != rval) return rval;
          if (!leaves.empty()) {
              // get closest leaf
            double min_dist = *dists.begin();
            closest_leaf = *leaves.begin();
            std::vector<EntityHandle>::iterator vit = leaves.begin()+1;
            std::vector<double>::iterator dit = dists.begin()+1;
            for (; vit != leaves.end() && min_dist; vit++, dit++) {
              if (*dit < min_dist) {
                min_dist = *dit;
                closest_leaf = *vit;
              }
            }
            dists.clear();
            leaves.clear();
          }
        }
        else {
          rval = myTree->point_search(pos+i3, closest_leaf);
          if (MB_ENTITY_NOT_FOUND == rval) closest_leaf = 0;
          else if (MB_SUCCESS != rval) return rval;
        }

          // if no ElemEvaluator, just return the box
        if (!elemEval) {
          ents[i] = closest_leaf;
          params[i3] = params[i3+1] = params[i3+2] = -2;
          if (is_inside && closest_leaf) is_inside[i] = true;
          continue;
        }
    
          // find natural coordinates of point in element(s) in that leaf
        CartVect tmp_nat_coords; 
        Range range_leaf;
        rval = mbImpl->get_entities_by_dimension(closest_leaf, myDim, range_leaf, false);
        if(rval != MB_SUCCESS) return rval;

          // loop over the range_leaf
        int tmp_inside;
        int *is_ptr = (is_inside ? is_inside+i : &tmp_inside);      
        *is_ptr = false;
        EntityHandle ent = 0;
        for(Range::iterator rit = range_leaf.begin(); rit != range_leaf.end(); rit++)
        {
          rval = elemEval->set_ent_handle(*rit); 
          if (MB_SUCCESS != rval) return rval;
          rval = elemEval->reverse_eval(pos+i3, tmp_abs_iter_tol, inside_tol, params+i3, is_ptr);
          if (MB_SUCCESS != rval) return rval;
          if (*is_ptr) {
            ent = *rit;
            break;
          }
        }
        if (debug && !ent) {
          std::cout << "Point " << i << " not found; point: (" 
                    << pos[i3] << "," << pos[i3+1] << "," << pos[i3+2] << ")" << std::endl;
          std::cout << "Source element candidates: " << std::endl;
          range_leaf.print("   ");
          for(Range::iterator rit = range_leaf.begin(); rit != range_leaf.end(); rit++)
          {
            std::cout << "Candidate " << CN::EntityTypeName(mbImpl->type_from_handle(*rit)) << " " << mbImpl->id_from_handle(*rit) << ": ";
            rval = elemEval->set_ent_handle(*rit); 
            if (MB_SUCCESS != rval) return rval;
            rval = elemEval->reverse_eval(pos+i3, tmp_abs_iter_tol, inside_tol, params+i3, is_ptr);
            if (MB_SUCCESS != rval) return rval;
            std::cout << "Parameters: (" << params[i3] << "," << params[i3+1] << "," << params[i3+2] << ")" 
                      << " inside = " << *is_ptr << std::endl;
          }
        }
        ents[i] = ent;
      }

      return MB_SUCCESS;
    }
    
        /* Count the number of located points in locTable
         * Return the number of entries in locTable that have non-zero entity handles, which
         * represents the number of points in targetEnts that were inside one element in sourceEnts
         *
         */
    int SpatialLocator::local_num_located() 
    {
      int num_located = locTable.get_n() - std::count(locTable.vul_rd, locTable.vul_rd+locTable.get_n(), 0);
      if (num_located != (int)locTable.get_n()) {
        unsigned long *nl = std::find(locTable.vul_rd, locTable.vul_rd+locTable.get_n(), 0);
        if (nl) {
          int idx = nl - locTable.vul_rd;
          if (idx) {}
        }
      }
      return num_located;
    }

        /* Count the number of located points in parLocTable
         * Return the number of entries in parLocTable that have a non-negative index in on a remote
         * proc in parLocTable, which gives the number of points located in at least one element in a
         * remote proc's sourceEnts.
         */
    int SpatialLocator::remote_num_located()
    {
      int located = 0;
      for (unsigned int i = 0; i < parLocTable.get_n(); i++)
        if (parLocTable.vi_rd[2*i] != -1) located++;
      return located;
    }
} // namespace moab

