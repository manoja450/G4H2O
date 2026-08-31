#ifndef G4d2oDataDrivenReflector_h
#define G4d2oDataDrivenReflector_h

#include "G4ThreeVector.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ios.hh"
#include "globals.hh"
#include <vector>
#include <map>
#include <random>

class TFile;
class TTree;

class G4d2oDataDrivenReflector {
public:
    G4d2oDataDrivenReflector(G4double reflectivity);
    ~G4d2oDataDrivenReflector();

    G4ThreeVector GetReflectionDirection(const G4ThreeVector& incomingDirection,
                                          const G4ThreeVector& normal,
                                          G4double incidentAngleRad) const;

    G4double SampleOutgoingAngle(G4double incidentAngleDeg) const;

    G4double GetReflectivity() const { return fReflectivity; }

    void LoadPDF(G4double incidentAngleDeg, const G4String& filename);
    void LoadAllPDFs(const G4String& directory);
    void PrintLoadedAngles() const;

    // Validation method
    void ValidatePDF(G4double incidentAngleDeg, G4int nSamples = 100000) const;

    static G4d2oDataDrivenReflector* GetInstance();
    static void CreateInstance(G4double reflectivity, const G4String& dataDirectory);
    static void DeleteInstance();

    static void SetRootFile(TFile* file);
    static void CloseRootTree();
    static void RecordReflection(G4double incidentAngleDeg, G4double reflectedAngleDeg,
                                 G4double px = 0, G4double py = 0, G4double pz = 0,
                                 G4double nx = 0, G4double ny = 0, G4double nz = 0);
    static void SetCurrentEventNumber(G4int eventNum);

private:
    struct PDF {
        std::vector<G4double> theta;        // reflection angle (degrees) - fine grid
        std::vector<G4double> intensity;    // intensity from thesis (interpolated)
        std::vector<G4double> cdf;          // cumulative distribution
        G4double maxIntensity;
    };

    std::map<G4int, PDF> fPDFs;
    G4double fReflectivity;
    mutable std::mt19937 fRNG;
    mutable std::uniform_real_distribution<G4double> fRandDist;

    // Helper functions
    void NormalizeAndBuildCDF(PDF& pdf) const;
    G4double SampleFromPDF(const PDF& pdf, G4double rand) const;
    G4double GetThetaOut(G4double incidentAngleDeg) const;
    G4double InterpolateThetaOut(G4double incidentAngleDeg, G4double rand) const;
    void FindBoundingPDFs(G4double incidentAngleDeg, const PDF*& pdfLow, const PDF*& pdfHigh,
                          G4double& weightLow, G4double& weightHigh) const;

    // Function B: Continuous interpolation helpers
    std::vector<G4double> CreateContinuousPDF(const std::vector<G4double>& theta,
                                               const std::vector<G4double>& intensity,
                                               G4int nFineBins = 360) const;
    G4double LinearInterpolate(const std::vector<G4double>& x,
                                const std::vector<G4double>& y,
                                G4double x0) const;

    static G4d2oDataDrivenReflector* fInstance;
    static TTree* fReflectionTree;
    static G4double fRecordedIncident;
    static G4double fRecordedReflected;
    static G4double fRecordedPx;
    static G4double fRecordedPy;
    static G4double fRecordedPz;
    static G4double fRecordedNx;
    static G4double fRecordedNy;
    static G4double fRecordedNz;
    static G4int fRecordedEventID;
    static G4int fReflectionCounter;
    static G4int fCurrentEventNumber;
};

#endif
