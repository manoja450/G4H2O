#include <algorithm>

#include "G4d2oStackingAction.hh"
#include "G4d2oDetector.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4Step.hh"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructor
G4d2oStackingAction::G4d2oStackingAction(G4d2oDetector * det) :
  fDetector(det)
{
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Destructor
G4d2oStackingAction::~G4d2oStackingAction()
{
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Prepare New event
void G4d2oStackingAction::PrepareNewEvent()
{
  this->pmt_instance_ids.clear();
  // Fill pmt_instance_ids vector with the instance ids from pmtPhysV_[i]
// in BeginOfRunAction
auto store = G4PhysicalVolumeStore::GetInstance();
for (auto* pv : *store)
  if ( pv->GetName().rfind("pmtPhysV_", 0) == 0 )   // faster prefix test
    pmt_instance_ids.emplace_back( pv->GetInstanceID() );

}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// New Stage
void G4d2oStackingAction::NewStage()
{
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Public Functions
G4ClassificationOfNewTrack G4d2oStackingAction::ClassifyNewTrack(const G4Track * the_track)
{
  G4ClassificationOfNewTrack classification    = G4ClassificationOfNewTrack::fWaiting;
  const G4VProcess * kprocess = the_track->GetCreatorProcess();

  if (the_track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() ) {
    G4double energy_eV = the_track->GetTotalEnergy() / eV;

    // Checks in WCSim StackngAction regardng process type etc...
    // Are those necessary? Or can we just kill based on energy?
    // Checks that the process type !=3 ...which I think just means not created from an optical process
    //  e.g. reflection?
    //  So probably needed so don't throw the dice multiple time for a single photon
// (b) fDetector pointer
if ( kprocess && kprocess->GetProcessType() != 3 && fDetector )  // ← guard
{
  if ( G4UniformRand() > fDetector->EvalPMTQE(energy_eV) )
    classification = G4ClassificationOfNewTrack::fKill;
}

    }
  

  // Remove primary particles generated within the pmt volume
// (a) the_track->GetVolume()
if ( the_track->GetParentID() == 0 )
{
  const auto* pv = the_track->GetVolume();      // ← can be nullptr
  if ( pv &&                                     //   protect here
       std::any_of( pmt_instance_ids.begin(),
                    pmt_instance_ids.end(),
                    [&](int id){ return pv->GetInstanceID() == id; } ) )
  {
    classification = G4ClassificationOfNewTrack::fKill;
  }
}

  return classification;
} // END OF USER STACKING ACTION
