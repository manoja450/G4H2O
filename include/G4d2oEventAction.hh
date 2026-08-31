#ifndef G4d2oEventAction_hh
#define G4d2oEventAction_hh 1

#include "G4d2oDetectorHit.hh"
#include "G4d2oPrimaryGeneratorAction.hh"
#include "inputVariables.hh"
#include "simEvent.h"

#include "TH2D.h"
#include "TVector3.h"

#include "G4UserEventAction.hh"
#include "G4Run.hh"
#include "G4Event.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "TGraph.h"

class G4d2oEventAction : public G4UserEventAction
{
protected:
    G4HCofThisEvent *HCoE;
    G4d2oDetectorHitsCollection *detHC[MAX_SEN_DET];
    G4d2oDetectorHitsCollection *sourceHC[MAX_SEN_DET];
    G4d2oDetectorHitsCollection *targHC[MAX_SEN_DET];

    G4d2oPrimaryGeneratorAction *thePGA;
    G4d2oDetector *theDet;

    G4int numEvents;

    TFile *fout;
    TTree *outTree;
    TH2D *hPMTArray;
    TH2D *hPMTNumVsTime;
    TH1D *hPhotonEnergy;
    TH2D *hSidePMT[2];
    TH1D *hTotalPhotons;

    TGraph *gBialkali_ev = 0;

    inputVariables *input;
    G4int numPMTs, numRows, numInEachRow;

    simEvent *theEventData;

    Int_t eventNumber;

    double qe_ev(double energy_ev);

    void GetHitsCollection(void);
    virtual void ProcessEvent(void);
    virtual void ZeroEventVariables(void);

public:
    G4d2oEventAction();
    ~G4d2oEventAction();

    static constexpr int kNarrowPanels = 4;
    void ResetNarrowVetoEdep()
    {
        for (int i = 0; i < kNarrowPanels; ++i)
            fNarrowVetoEdep[i] = 0.0;
    }
    void AddNarrowVetoEdep(int idx, G4double e)
    {
        if (0 <= idx && idx < kNarrowPanels)
            fNarrowVetoEdep[idx] += e;
    }
    const double *GetNarrowVetoEdep() const { return fNarrowVetoEdep; }

    static constexpr int kWidePanels = 4;
    void ResetWideVetoEdep()
    {
        for (int i = 0; i < kWidePanels; ++i)
            fWideVetoEdep[i] = 0.0;
    }
    void AddWideVetoEdep(int idx, G4double e)
    {
        if (0 <= idx && idx < kWidePanels)
            fWideVetoEdep[idx] += e;
    }
    const double *GetWideVetoEdep() const { return fWideVetoEdep; }

    void AddVetoEdep(G4double e) { fVetoEdep += e; }
    G4double GetVetoEdep() const { return fVetoEdep; }

    void BeginOfEventAction(const G4Event *thisEvent) override;
    void EndOfEventAction(const G4Event *thisEvent) override;

    TFile *OutFilePtr() { return fout; }

private:
    G4double fVetoEdep = 0.0;
    double fNarrowVetoEdep[kNarrowPanels] = {0,0,0,0};
    double fWideVetoEdep[kWidePanels] = {0,0,0,0};

}; // END of class G4d2oEventAction

#endif
