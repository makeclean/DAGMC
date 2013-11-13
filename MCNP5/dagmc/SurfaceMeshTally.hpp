// MCNP5/dagmc/TrackLengthMeshTally.hpp

#include <string>
#include <cassert>
#include <set>

#include "moab/Interface.hpp"
#include "moab/CartVect.hpp"
#include "moab/Range.hpp"

#include "MeshTally.hpp"
#include "TallyEvent.hpp"

namespace moab{

/* Forward Declarations */
class AdaptiveKDTree;
class OrientedBoxTreeTool;

//===========================================================================//
/**
 * \class TrackLengthMeshTally
 * \brief Represents an unstructured mesh tally based on particle tracks
 *
 * TrackLengthMeshTally is a concrete class derived from MeshTally that
 * implements the Tally interface to tally particle tracks on an unstructured
 * mesh as part of a Monte Carlo particle transport simulation.  If a
 * TrackLengthMeshTally object receives a TallyEvent type that is not
 * TallyEvent::TRACK, then no scores are computed.
 *
 * ==========
 * TallyInput
 * ==========
 *
 * The TallyInput struct needed to construct a TrackLengthMeshTally object is
 * defined in Tally.hpp and is set through the TallyManager when a Tally is
 * created.  Options that are currently available for TrackLengthMeshTally
 * objects include
 *
 * 1) "inp"="input_filename", "out"="output_filename"
 * --------------------------------------------------
 * These two options are processed through the MeshTally constructor.  The
 * "inp" key is REQUIRED for TrackLengthMeshTally objects, whereas the "out"
 * key is optional.  See MeshTally.hpp for more information.
 *
 * 2) "convex"="t/f", "convex"="true/fale"
 * ----------------------------------------------
 * If the convex option is set to true, the user is asserting that the input
 * mesh has convex geometry.  Single particle tracks that leave a convex mesh
 * are no longer tracked as they cannot re-enter the mesh.  The default value
 * for this option is false.
 *
 * 3) "conformal"="values"
 * -----------------------
 * If the input mesh is exactly conformal to the geometry, then conformal
 * logic should be used for scoring particle tracks across mesh cells or
 * errors may occur. This option allows the user to identify what geometric
 * cells are conformal to the input mesh.  These values can be defined as a
 * list of cells separated by commas (i.e. 1, 2, 5), and/or a range of cells
 * (i.e. 2-5).  Like the convex case, particle tracks that leave a conformal
 * mesh are no longer tracked.
 *
 * 4) "conf_surf_src"="t/f", "conf_surf_src"="true/false"
 * ------------------------------------------------------
 * If the conformal surface source option is set to true, then the user is
 * asserting that a surface source has been defined that is conformal to a
 * geometric surface.  This means that new particles that are born from
 * this source will use conformal logic to determine the first mesh cell
 * that is crossed.  The default value for this option is false.
 *
 * 5) "tag"="name", "tagval"="value"
 * ---------------------------------
 * The "tag" and "tagval" options can be used to define a subset of the mesh
 * elements in the input file that will be tallied.  These tags are defined
 * on the input mesh itself using the MOAB tagging feature.  Note that "tag"
 * name can only be set once, whereas multiple "tagval" values can be added.
 * This option is only used during setup to define the set of tally points.
 */
//===========================================================================//
class SurfaceMeshTally : public MeshTally
{ 
  public:
    /**
     * \brief Constructor
     * \param input user-defined input parameters for this SurfaceMeshTally
     */
    explicit SurfaceMeshTally(const TallyInput& input);

    /**
     * \brief Virtual destructor
     */
    virtual ~SurfaceMeshTally();

    // >>> DERIVED PUBLIC INTERFACE from Tally.hpp

    /**
     * \brief Computes scores for this SurfaceMesh
     * \param event the parameters needed to compute the scores
     */
    virtual void compute_score(const TallyEvent& event);

    /**
     * \brief Updates TrackLengthMeshTally when a particle history ends
     *
     * Calls MeshTally::end_history() and if input mesh is conformal sets the
     * last_cell to -1 to indicate that a new particle will be born.
     */
    virtual void end_history();

    /**
     * \brief Write results to the output file for this SurfaceMeshTally
     * \param num_histories the number of particle histories tracked
     *
     * The write_data() method writes the current tally and relative standard
     * error results for all of the mesh cells defined as tally points to the
     * output_filename set for this SurfaceMeshTally.  These values are
     * normalized by both the number of particle histories that were tracked
     * and the volume of the mesh cell for which the results were computed.
     */
    virtual void write_data(double num_histories);

  protected:
    /// Copy constructor and operator= methods are not implemented
    SurfaceMeshTally(const SurfaceMeshTally& obj);
    SurfaceMeshTally& operator=(const SurfaceMeshTally& obj);

  protected:
    // MOAB instance that stores all of the mesh data
    moab::Interface* mb;

    // KD-Tree used with unstructured meshes to get starting tetrahedron
    moab::AdaptiveKDTree* kdtree;  
    moab::EntityHandle kdtree_root;

    // Stores normal data for each facet
    std::vector<CartVect> surface_normals;
    // Stores the surface area for each facet
    std::vector<double> surface_areas;

    // >>> PROTECTED METHODS

    /**
     * \brief Parse the TallyInput options for this TrackLengthMeshTally
     */
  //    void parse_tally_options();

    /**
     * \brief Loads and sets the mesh data used by this TrackLengthMeshTally
     *
     * Mesh data is loaded from the input_filename provided in the TallyInput
     * options.  If no tag name and values are found, the set of tally points
     * will be all of the mesh cells defined in input_filename.  Otherwise,
     * the set of tally points stored in MeshTally::tally_mesh_set will only
     * contain the mesh cells that have the given tag name and tag values.
     */
    void set_tally_meshset();

    /**
     * \brief Computes the normal vector for every facet in the mesh
     * \param all_facets, all facets in the problem
     * \return the MOAB ErrorCode value
     */
    ErrorCode compute_surface_normals(const Range &all_facets);

    /**
     * \brief Computes the surface areas for each facet in the problem
     * \param all_facets, all facets in the problem
     * \return the MOAB ErrorCode value
     */
    ErrorCode compute_surface_areas(const Range &all_facets);

    /**
     * \brief Constructs the KDtrees from the mesh data
     * \param all_facets the set of facets extracted from the input mesh
     */
    void build_trees (Range& all_tets);

    /**
     * \brief gets the hits on the kdtree 
     * \param direction, the direction of the ray
     * \param position, the location of the origin of the ray
     * \param track_length, the length of the ray
     */
    std::vector<EntityHandle> get_intersections(const double direction[3],
						const double position[3],
						const double track_length);

    /**
     * \brief computes the surface normal given an entity handle
     * \param triangle, a given MBEntityHandle 
     */
    CartVect surface_normal(const EntityHandle triangle); 

    /**
     * \brief computes the surface area of the facet
     * \param triangle, a given MBEntityHandle 
     */
    double surface_area(const EntityHandle triangle); 

    /**
     * \brief computes the angle between two vectors 
     * \param normal, normal vector of a surface
     * \param direction, direction vector or a ray
     */ 
     double get_angle(const CartVect normal, const CartVect direction);

};

} // end namespace moab

// end of MCNP5/dagmc/TrackLengthMeshTally.hpp
