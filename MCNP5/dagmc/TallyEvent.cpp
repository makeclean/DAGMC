// MCNP5/dagmc/TallyEvent.cpp

#include <cstdlib>
#include <iostream>
#include <utility>

#include "moab/CartVect.hpp"

#include "TallyEvent.hpp"

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
TallyEvent::TallyEvent()
    : track_length(0.0), position(moab::CartVect(0.0, 0.0, 0.0)),
      direction(moab::CartVect(0.0, 0.0, 0.0)), total_cross_section(0.0),
      particle_energy(0.0), particle_weight(1.0), tally_multiplier(1.0)
{
    // set this tally event to begin with no type
    type = NONE;
}
//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
void TallyEvent::set_track_event(double track_length,
                                 double x, double y, double z,
                                 double u, double v, double w,
                                 double particle_energy,
                                 double particle_weight)
{
    if (type != NONE)
    {
        std::cerr << "Error: Tally event has already been set" << std::endl;
        exit(EXIT_FAILURE);
    }

    // set variables for track-based event
    this->track_length = track_length;
    this->position = moab::CartVect(x, y, z);
    this->direction = moab::CartVect(u, v, w);
    this->particle_energy = particle_energy;
    this->particle_weight = particle_weight;

    // set event type
    type = TRACK;
}
//---------------------------------------------------------------------------//
void TallyEvent::set_collision_event(double total_cross_section,
                                     double x, double y, double z,
                                     double particle_energy,
                                     double particle_weight)
{
    if (type != NONE)
    {
        std::cerr << "Error: Tally event has already been set" << std::endl;
        exit(EXIT_FAILURE);
    }

    // set variables for collision event
    this->total_cross_section = total_cross_section;
    this->position = moab::CartVect(x, y, z);
    this->particle_energy = particle_energy;
    this->particle_weight = particle_weight;

    // set event type
    type = COLLISION;
}
//---------------------------------------------------------------------------//
void TallyEvent::set_tally_multiplier(double value)
{
    tally_multiplier = value;
}
//---------------------------------------------------------------------------//
double TallyEvent::get_tally_multiplier() const
{
    return tally_multiplier;
}
//---------------------------------------------------------------------------//
double TallyEvent::get_weighting_factor() const
{
    return tally_multiplier * particle_weight;
}
//---------------------------------------------------------------------------//

// end of MCNP5/dagmc/TallyEvent.cpp
