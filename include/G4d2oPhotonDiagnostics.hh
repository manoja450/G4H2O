#ifndef G4d2oPhotonDiagnostics_h
#define G4d2oPhotonDiagnostics_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <map>

// Termination reasons - these tell you where photons actually die
enum TerminationReason {
    kReachedPMT = 0,
    kWaterAbsorbed,
    kTyvekAbsorbed,
    kMaxReflections,
    kBoundaryKill,
    kOther
};

class G4d2oPhotonDiagnostics {
public:
    static G4d2oPhotonDiagnostics* GetInstance();
    static void CreateInstance();
    static void DeleteInstance();

    // Per-photon tracking
    void StartPhoton(G4int eventID, G4int photonID, G4double energy,
                     const G4ThreeVector& position, const G4ThreeVector& direction);
    void RecordReflection(G4int reflectionNumber,
                          const G4ThreeVector& position,
                          const G4ThreeVector& normal,
                          const G4ThreeVector& incomingDir,
                          const G4ThreeVector& outgoingDir,
                          G4double incidentAngleDeg,
                          G4double reflectedAngleDeg,
                          G4double globalPz,
                          G4double dotProduct);
    void RecordAbsorption(const G4ThreeVector& position, const G4String& reason);
    void RecordPMTHit(const G4ThreeVector& position, G4double time);
    void EndPhoton(G4bool reachedPMT, G4int totalReflections);

    // Per-event statistics
    void NewEvent(G4int eventID);
    void RecordEventResult(G4int eventID, G4int nPhotons, G4int nPMTHits,
                           G4int nAbsorbed, G4double totalPE, G4int qualityCutPass,
                           G4double meanReflections, G4int maxReflections);

    // Photon Survival Chain - these are the counters you care about
    void RecordCherenkovPhoton(G4double energy, const G4ThreeVector& position);
    void RecordTyvekReflection();
    void RecordPMTReach();
    void RecordPhotoelectron();
    void RecordWaterAbsorption(const G4ThreeVector& position);
    void RecordEventChain(G4int eventID, G4int nCherenkov,
                          G4int nTyvekReflections, G4int nPMTReach, G4int nPE);

    // Reflection count per photon (for the max reflection limit)
    G4int GetReflectionCountForPhoton(G4int photonID) const;
    void IncrementReflectionCountForPhoton(G4int photonID);
    void ResetReflectionCountForPhoton(G4int photonID);
    void ClearAllPhotonCounts();

    // Termination tracking - THIS IS WHERE THE PHYSICS IS RECORDED
    void RecordTermination(G4int photonID, TerminationReason reason);
    void PrintTerminationSummary() const;

    // Print methods
    void PrintPhotonTrace(G4int eventID, G4int photonID) const;
    void PrintEventSummary(G4int eventID) const;
    void PrintGlobalSummary() const;
    void PrintBadReflections() const;
    void PrintChainSummary() const;  // THIS NOW PRINTS PHYSICALLY CORRECT CHAIN

    void Reset();
    void SetPrintTraces(G4bool print) { fPrintTraces = print; }
    void SetMaxTraces(G4int max) { fMaxTraces = max; }
    void SetTraceEvent(G4int eventID) { fTraceEventID = eventID; }
    void WriteDiagnostics();

private:
    G4d2oPhotonDiagnostics();
    ~G4d2oPhotonDiagnostics();

    struct PhotonTrace {
        G4int eventID;
        G4int photonID;
        G4double initialEnergy;
        G4ThreeVector startPos;
        G4ThreeVector startDir;
        std::vector<G4double> reflectionAngles;
        std::vector<G4double> incidentAngles;
        std::vector<G4ThreeVector> reflectionPositions;
        std::vector<G4ThreeVector> reflectionOutgoing;
        std::vector<G4double> globalPzValues;
        std::vector<G4double> dotProducts;
        G4int nReflections;
        bool reachedPMT;
        G4double finalTime;
        TerminationReason terminationReason;
    };

    struct EventSummary {
        G4int eventID;
        G4int nPhotons;
        G4int nPMTHits;
        G4int nAbsorbed;
        G4double totalPE;
        G4int qualityCutPass;
        G4double meanReflections;
        G4int maxReflections;
        G4double fractionBadReflections;
    };

    struct EventChain {
        G4int eventID;
        G4int nCherenkov;
        G4int nTyvekReflections;
        G4int nPMTReach;
        G4int nPE;
    };

    static G4d2oPhotonDiagnostics* fInstance;

    std::map<G4int, std::map<G4int, PhotonTrace>> fPhotonTraces;
    std::map<G4int, EventSummary> fEventSummaries;
    std::map<G4int, EventChain> fEventChains;
    std::map<G4int, G4int> fPhotonReflectionCounts;
    std::map<G4int, TerminationReason> fPhotonTermination;
    std::map<TerminationReason, G4int> fTerminationCounts;

    G4int fCurrentEventID;
    G4int fCurrentPhotonID;
    PhotonTrace fCurrentTrace;

    G4int fTotalEvents;
    G4int fTotalPhotons;
    G4int fTotalReflections;
    G4int fTotalPMTHits;
    G4int fTotalAbsorbed;
    G4int fTotalBadReflections;
    G4int fTotalQualityPass;
    G4int fTotalCherenkovPhotons;
    G4int fTotalTyvekReflections;
    G4int fTotalPMTReach;
    G4int fTotalPhotoelectrons;

    G4bool fPrintTraces;
    G4int fMaxTraces;
    G4int fTraceEventID;
    G4int fTracesPrinted;
};

#endif
