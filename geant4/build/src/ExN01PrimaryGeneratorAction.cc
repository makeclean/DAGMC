//
//
// $Id$
//

#include "ExN01PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
//
#include "G4GeneralParticleSource.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

#include "source_tf.h"

ExN01PrimaryGeneratorAction::ExN01PrimaryGeneratorAction()
  :G4VUserPrimaryGeneratorAction(), particleGun(0)
{
  //  fParticleGun = new G4GeneralParticleSource();
  //fParticleGun->SetParticleDefinition(particleTable->FindParticle(particleName="neutron"));
  //fParticleGun->SetParticleEnergy(14.0*MeV);
  
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="neutron"));
  setup_();
}

ExN01PrimaryGeneratorAction::~ExN01PrimaryGeneratorAction()
{
  delete particleGun;
}

void ExN01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //  fParticleGun->GeneratePrimaryVertex(anEvent) ;

  // 
  
  G4double x,y,z;
  G4double energy,weight;
  x = 350.0;
  G4double y_min = -40.0;
  G4double y_max = 40.0;
  G4double z_min = -400.0;
  G4double z_max = 400.0;
  G4double rand1 = G4UniformRand();
  G4double rand2 = G4UniformRand();

  //  sample_linear_(20.,-20.,G4UniformRand(),x);
  linear_sample_(y_max,y_min,rand1,y);
  linear_sample_(z_max,z_min,rand2,z);

  G4double rand3 = G4UniformRand();
  G4double rand4 = G4UniformRand();

  sample_(rand3,rand4,energy,weight);
  particleGun->SetParticleEnergy(energy*MeV);
  particleGun->SetParticlePosition(G4ThreeVector(x*cm,y*cm,z*cm));

  //  particleGun->SetParticlePosition(G4ThreeVector(-10.0*cm,y*cm,z*cm));
  //  G4int i = anEvent->GetEventID() % 3;

  G4double cosTheta = 2*G4UniformRand() - 1.;
  G4double phi = twopi*G4UniformRand();
  G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
  G4double ux = sinTheta*std::cos(phi),
    uy = sinTheta*std::sin(phi),
    uz = cosTheta;
  particleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));
  particleGun->GeneratePrimaryVertex(anEvent);
}


