#ifndef G4d2oSteppingAction_h
#define G4d2oSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4d2oSteppingAction : public G4UserSteppingAction {
public:
    G4d2oSteppingAction();
    virtual ~G4d2oSteppingAction();

    virtual void UserSteppingAction(const G4Step* step) override;

    void ResetEventCounters();
    void SetEventID(G4int id) { fEventID = id; }

    // For Sim_Tree instrumentation
    G4int GetReflectionCount() const { return fReflectionCount; }
    G4double GetTotalPathLength() const { return fTotalPathLength; }
    G4int GetPhotonCount() const { return fPhotonCount; }

private:
    // Per-event counters (for Sim_Tree)
    G4int fReflectionCount;
    G4double fTotalPathLength;
    G4int fPhotonCount;

    // Per-photon tracking
    G4int fCurrentPhotonID;
    G4ThreeVector fPrevDirection;
    G4bool fHasPrevDirection;

    // Event ID for ReflectionTree
    G4int fEventID;
};

#endif
