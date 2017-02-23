typedef std::vector< std::pair< EntityHandle, double > > Ray_History;
typedef std::vector< std::pair< EntityHandele, double> >::iterator Ray_History_iter;

class Ray_Urchin {

private:
  DagMC* dagmc;

public:
  Ray_Urchin(const std::string cad_file,
             const std::string ray_file,
             const CartVect start_pt_in,
             const int start_vol_id);
  
};


