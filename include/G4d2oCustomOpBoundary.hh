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
        kUniform,           // phi = Uniform(0, 2pi) [thesis default assumption]
        kGaussian15,        // phi ~ Gaussian(0, 15deg) [narrow spread hypothesis]
        kGaussian30,        // phi ~ Gaussian(0, 30deg) [thesis Super-K comparison value]
        kGaussian35,        // phi ~ Gaussian(0, 35deg) [test slightly wider than 30deg]
        kGaussian45,        // phi ~ Gaussian(0, 45deg) [moderate spread hypothesis]
        kGaussian60,        // phi ~ Gaussian(0, 60deg) [broad spread hypothesis]
        kLambertianAzimuth  // phi weighted by cos(phi) [Lambertian-like hypothesis]
    };

    G4d2oCustomOpBoundary(const G4String& processName = "G4d2oCustomOpBoundary");
    virtual ~G4d2oCustomOpBoundary();

    virtual G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step) override;

    // Set azimuthal model (call before simulation starts)
    void SetAzimuthalModel(AzimuthalModel model) { fAzimuthalModel = model; }
    AzimuthalModel GetAzimuthalModel() const { return fAzimuthalModel; }

    // Get model name as string (for printing)
    const char* GetAzimuthalModelName() const;

private:
    // Helper: Sample azimuthal angle based on selected model
    G4double SampleAzimuthalAngle() const;

    // Helper: surface normal at a step point, computed correctly for
    // translated/rotated volumes (e.g. tyvekCapBotPhysV) by transforming
    // into the solid's local frame before calling SurfaceNormal(), then
    // rotating the result back to global. See .cc for the full rationale.
    static G4ThreeVector GetLocalFrameSurfaceNormal(const G4StepPoint* point);

    // Random number generators
    mutable std::mt19937 fRNG;
    mutable std::uniform_real_distribution<G4double> fRandDist;
    mutable std::normal_distribution<G4double> fGaussDist;

    // Current azimuthal model
    AzimuthalModel fAzimuthalModel = kUniform;  // Default: uniform
};

#endif

