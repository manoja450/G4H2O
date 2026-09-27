#ifndef G4d2oCustomOpBoundary_h
#define G4d2oCustomOpBoundary_h

#include "G4OpBoundaryProcess.hh"
#include "G4d2oDataDrivenReflector.hh"
#include <random>

class G4d2oCustomOpBoundary : public G4OpBoundaryProcess {
public:
    // Azimuthal model selection
    // The thesis only measures the in-plane (1D) distribution.
    // The out-of-plane (azimuthal) distribution is NOT measured.
    // These models represent different assumptions for extending 1D->3D.
    enum AzimuthalModel {
        kUniform,           
        kGaussian15,        
        kGaussian30,        
        kGaussian35,        
        kGaussian45,        
        kGaussian60,        
        kLambertianAzimuth  
    };

    G4d2oCustomOpBoundary(const G4String& processName = "G4d2oCustomOpBoundary");
    virtual ~G4d2oCustomOpBoundary();

    virtual G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step) override;


    void SetAzimuthalModel(AzimuthalModel model) { fAzimuthalModel = model; }
    AzimuthalModel GetAzimuthalModel() const { return fAzimuthalModel; }

 
    const char* GetAzimuthalModelName() const;

private:
  
    G4double SampleAzimuthalAngle() const;

    
    static G4ThreeVector GetLocalFrameSurfaceNormal(const G4StepPoint* point);

    mutable std::mt19937 fRNG;
    mutable std::uniform_real_distribution<G4double> fRandDist;
    mutable std::normal_distribution<G4double> fGaussDist;
 
    AzimuthalModel fAzimuthalModel = kUniform;  // Default: uniform
};

#endif

