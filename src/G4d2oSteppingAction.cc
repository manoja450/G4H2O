#include "G4d2oSteppingAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>

G4d2oSteppingAction::G4d2oSteppingAction()
    : fReflectionCount(0), fTotalPathLength(0.0), fPhotonCount(0),
      fCurrentPhotonID(-1), fHasPrevDirection(false), fEventID(0) {}

G4d2oSteppingAction::~G4d2oSteppingAction() {}

void G4d2oSteppingAction::ResetEventCounters() {
    fReflectionCount = 0;
    fTotalPathLength = 0.0;
    fPhotonCount = 0;
    fCurrentPhotonID = -1;
    fHasPrevDirection = false;
}

void G4d2oSteppingAction::UserSteppingAction(const G4Step* step) {
    G4Track* track = step->GetTrack();
    G4ParticleDefinition* particle = track->GetDefinition();
    G4String particleName = particle->GetParticleName();

    if (particleName != "opticalphoton") return;

    G4int trackID = track->GetTrackID();
    G4double stepLength = step->GetStepLength();

    if (trackID != fCurrentPhotonID) {
        fCurrentPhotonID = trackID;
        fPhotonCount++;
        fHasPrevDirection = false;
    }

    fTotalPathLength += stepLength;

    G4ThreeVector currentDir = track->GetMomentumDirection();
    if (fHasPrevDirection) {
        G4double dot = fPrevDirection.dot(currentDir);

        // ============================================================
        //  only count reflections if direction changed AND we crossed a real boundary
        // ============================================================
        if (dot < 0.998 && step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
            fReflectionCount++;
            // ReflectionTree is filled inside G4d2oCustomOpBoundary for the
            // Data-Driven model - nothing else to do here.
        }
    }
    fPrevDirection = currentDir;
    fHasPrevDirection = true;
}
