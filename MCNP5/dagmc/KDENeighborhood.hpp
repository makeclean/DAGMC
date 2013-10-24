// MCNP5/dagmc/KDENeighborhood.hpp

#ifndef DAGMC_KDE_NEIGHBORHOOD_HPP
#define DAGMC_KDE_NEIGHBORHOOD_HPP

#include <set>

#include "moab/Interface.hpp"

#include "TallyEvent.hpp"

// forward declarations
namespace moab {
  class AdaptiveKDTree;
  class CartVect;
}

//===========================================================================//
/**
 * \class KDENeighborhood
 * \brief Defines a neighborhood region for a KDE mesh tally event
 *
 * KDENeighborhood is a class that defines the neighborhood region for a
 * KDEMeshTally as part of a Monte Carlo particle transport simulation.  The
 * neighborhood region is formally defined as the region in space for which
 * the kernel function produces a non-trivial result for any mesh node that
 * exists inside that region.  This set of mesh nodes is also known as the
 * set of calculation points for the KDEMeshTally.
 *
 * In general, it is not always easy to define the exact neighborhood region.
 * Therefore, the default behavior of KDENeighborhood is to use a kd-tree
 * search method to locate all possible calculation points for each TallyEvent.
 * This kd-tree approach produces an exact neighborhood region for collision
 * events, but only an approximation for track-based events.
 *
 * =============================
 * KDENeighborhood Functionality
 * =============================
 *
 * Once a KDENeighborhood has been created, there are two steps that are
 * needed to get its calculation points.  Since the dimensions of the exact
 * neighborhood region usually changes with each TallyEvent, it is first
 * necessary to call update_neighborhood().  Once the neighborhood has been
 * updated, then the set of calculation points associated with that event can
 * be obtained by get_points().
 */
//===========================================================================//
class KDENeighborhood
{
  public:
    /**
     * \brief Constructor
     * \param[in] mbi pointer to a pre-loaded MOAB instance
     * \param[in] mesh_set the MOAB mesh set containing all of the mesh data
     * \param[in] build_kd_tree if true, creates kd-tree based on mesh nodes
     *
     * Note that setting build_kd_tree to false forces the KDENeighborhood to
     * always use all calculation points with every TallyEvent that occurs.
     */
    KDENeighborhood(moab::Interface* mbi,
                    const moab::EntityHandle& mesh_set,
                    bool build_kd_tree = true);

    /**
     * \brief Destructor
     */
    ~KDENeighborhood();

    // >>> PUBLIC INTERFACE

    /**
     * \brief Gets the calculation points for this neighborhood region
     * \return set of calculation points currently in the neighborhood region
     */
    std::set<moab::EntityHandle> get_points() const;

    /**
     * \brief Updates the neighborhood region based on the given tally event
     * \param[in] event the tally event for which the neighborhood is desired
     * \param[in] bandwidth the bandwidth vector (hx, hy, hz)
     */
    void update_neighborhood(const TallyEvent& event,
                             const moab::CartVect& bandwidth);

    /**
     * \brief checks if a point exists within this neighborhood region
     * \param[in] coords the coordinates of the point to check
     * \return true if point does exist within this neighborhood region
     */
    bool point_in_region(const moab::CartVect& coords) const;

    /**
     * \brief Determines if point lies within radius of cylindrical region
     * \param[in] point the calculation point to be tested
     * \return true if point is inside the region, false otherwise
     *
     * Note that this method is only valid for track-based events.  If the
     * event is not a track-based event, then it will always return false.
     * It will also return false if no radius was set due to the neighborhood
     * containing all calculation points.
     */
    bool point_within_max_radius(const TallyEvent& event,
                                 const moab::CartVect& point) const;

  private:
    // Set of calculation points currently in this neighborhood region
    std::set<moab::EntityHandle> points;

    // KD-Tree containing all mesh nodes in the input mesh
    moab::AdaptiveKDTree* kd_tree;
    moab::EntityHandle kd_tree_root;

    // Minimum and maximum corner of a rectangular neighborhood region
    double min_corner[3];
    double max_corner[3];

    // Radius of a cylindrical neighborhood region
    double radius;

    // >>> PRIVATE METHODS

    /**
     * \brief Sets the neighborhood region for a collision event
     * \param[in] collision_point the location of the collision (x, y, z)
     * \param[in] bandwidth the bandwidth vector (hx, hy, hz)
     */
    void set_neighborhood(const moab::CartVect& collision_point,
                          const moab::CartVect& bandwidth);

    /**
     * \brief Sets the neighborhood region for a track-based event
     * \param[in] track_length the total length of the track segment
     * \param[in] start_point the starting location of the particle (xo, yo, zo)
     * \param[in] direction the direction the particle is traveling (uo, vo, wo)
     * \param[in] bandwidth the bandwidth vector (hx, hy, hz)
     */
    void set_neighborhood(double track_length,
                          const moab::CartVect& start_point,
                          const moab::CartVect& direction,
                          const moab::CartVect& bandwidth);

    /**
     * \brief Finds the vertices that exist inside a rectangular region
     *
     * Includes vertices that are within +/- 1e-12 of a box boundary.  This
     * method updates the set of calculation points with all vertices that
     * were located within the current neighborhood region.
     */
    void points_in_box();
};

#endif // DAGMC_KDE_NEIGHBORHOOD_HPP

// end of MCNP5/dagmc/KDENeighborhood.hpp
