#include "G4d2oCustomOpBoundary.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpticalPhoton.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ParticleChange.hh"
#include "G4OpticalSurface.hh"
#include "G4TouchableHandle.hh"
#include "G4AffineTransform.hh"
#include "G4NavigationHistory.hh"
#include <cstdio>

static G4int gPrintCount = 0;
static G4int gMaxPrint = 0;
static G4bool gPrintLimitReached = false;
static G4int gReflectionCount = 0;

static G4long gTyvekCrossings = 0;
static G4long gTyvekReflections = 0;
static G4long gOtherCrossings = 0;
static G4long gOtherReflections = 0;

static G4bool IsTyvekVolume(const G4VPhysicalVolume* vol) {
    if (!vol) return false;
    const G4String& name = vol->GetName();
    return name.find("tyvek") != std::string::npos || name.find("Tyvek") != std::string::npos;
}

// ============================================================
// Constructor / Destructor
// ============================================================

G4d2oCustomOpBoundary::G4d2oCustomOpBoundary(const G4String& processName)
    : G4OpBoundaryProcess(processName),
      fRNG(std::random_device{}()),
      fRandDist(0.0, 1.0),
      fGaussDist(0.0, 1.0) {
}

G4d2oCustomOpBoundary::~G4d2oCustomOpBoundary() {
}

const char* G4d2oCustomOpBoundary::GetAzimuthalModelName() const {
    switch(fAzimuthalModel) {
        case kUniform:           return "UNIFORM (0 to 2pi)";
        case kGaussian15:        return "GAUSSIAN sigma = 15deg";
        case kGaussian30:        return "GAUSSIAN sigma = 30deg";
        case kGaussian35:        return "GAUSSIAN sigma = 35deg";
        case kGaussian45:        return "GAUSSIAN sigma = 45deg";
        case kGaussian60:        return "GAUSSIAN sigma = 60deg";
        case kLambertianAzimuth: return "LAMBERTIAN (cos-weighted)";
        default:                 return "UNKNOWN";
    }
}

G4double G4d2oCustomOpBoundary::SampleAzimuthalAngle() const {
    G4double phi = 0.0;

    switch(fAzimuthalModel) {
        case kUniform:
    
            return 2.0 * M_PI * fRandDist(fRNG);

        case kGaussian15:
            phi = fGaussDist(fRNG) * 15.0 * deg;
            break;

        case kGaussian30:
            phi = fGaussDist(fRNG) * 30.0 * deg;
            break;

        case kGaussian35:
            phi = fGaussDist(fRNG) * 35.0 * deg;
            break;

        case kGaussian45:
            phi = fGaussDist(fRNG) * 45.0 * deg;
            break;

        case kGaussian60:
            phi = fGaussDist(fRNG) * 60.0 * deg;
            break;

        case kLambertianAzimuth:
            {
                G4double phi_abs = 0.0;
                G4double maxVal = 1.0;
                int nAttempts = 0;
                while (nAttempts < 1000) {
                    G4double trial = fRandDist(fRNG) * 90.0 * deg;
                    G4double prob = std::cos(trial);
                    if (fRandDist(fRNG) * maxVal <= prob) {
                        phi_abs = trial;
                        break;
                    }
                    nAttempts++;
                }
                phi = phi_abs;
                if (fRandDist(fRNG) < 0.5) phi = -phi;
            }
            break;

        default:
            return 2.0 * M_PI * fRandDist(fRNG);
    }

    if (phi > 90.0 * deg) phi = 90.0 * deg;
    if (phi < -90.0 * deg) phi = -90.0 * deg;

    return phi;
}

G4ThreeVector G4d2oCustomOpBoundary::GetLocalFrameSurfaceNormal(const G4StepPoint* point) {
    if (!point) return G4ThreeVector(0, 0, 1);

    G4VPhysicalVolume* volume = point->GetPhysicalVolume();
    if (!volume || !volume->GetLogicalVolume() || !volume->GetLogicalVolume()->GetSolid()) {
        return G4ThreeVector(0, 0, 1);
    }

    G4TouchableHandle touchable = point->GetTouchableHandle();
    if (!touchable || !touchable->GetHistory()) {
        
        return volume->GetLogicalVolume()->GetSolid()->SurfaceNormal(point->GetPosition()).unit();
    }

    G4AffineTransform globalToLocal = touchable->GetHistory()->GetTopTransform();
    G4ThreeVector localPos = globalToLocal.TransformPoint(point->GetPosition());

    G4ThreeVector localNormal =
        volume->GetLogicalVolume()->GetSolid()->SurfaceNormal(localPos);

    G4AffineTransform localToGlobal = globalToLocal.Inverse();
    G4ThreeVector globalNormal = localToGlobal.TransformAxis(localNormal);

    return globalNormal.unit();
}

G4VParticleChange* G4d2oCustomOpBoundary::PostStepDoIt(const G4Track& track,
                                                         const G4Step& step) {

    // --- 1. Store original direction ---
    G4ThreeVector originalDir = track.GetMomentumDirection();

    // --- 2. Let Geant4 handle the boundary (REFLECTIVITY/TRANSMITTANCE
    // decision happens here) ---
    G4VParticleChange* vpChange = G4OpBoundaryProcess::PostStepDoIt(track, step);

    // --- 3. Cast to G4ParticleChange ---
    G4ParticleChange* particleChange = dynamic_cast<G4ParticleChange*>(vpChange);
    if (!particleChange) return vpChange;

    const G4StepPoint* postStepPoint = step.GetPostStepPoint();
    const G4StepPoint* preStepPoint = step.GetPreStepPoint();
    G4bool isTyvekBoundary = IsTyvekVolume(postStepPoint->GetPhysicalVolume()) ||
                              IsTyvekVolume(preStepPoint->GetPhysicalVolume());

    if (isTyvekBoundary) gTyvekCrossings++; else gOtherCrossings++;

    // --- 4. Check if reflection occurred ---
    if (*particleChange->GetMomentumDirection() == originalDir) {
        return vpChange;
    }

    if (isTyvekBoundary) gTyvekReflections++; else gOtherReflections++;

    // Only apply the thesis-sampled Tyvek angular model at Tyvek
    // boundaries. Anywhere else, this process reflects; leave
    // Geant4's own decision alone - the thesis data describes Tyvek,
    // not those surfaces.
    if (!isTyvekBoundary) {
        return vpChange;
    }

    // --- 5. Get reflector ---
    G4d2oDataDrivenReflector* reflector = G4d2oDataDrivenReflector::GetInstance();
    if (reflector == nullptr) return vpChange;

    // --- 6. Get incident angle ---
    G4ThreeVector incomingDir = track.GetMomentumDirection();

    G4ThreeVector normal = GetLocalFrameSurfaceNormal(postStepPoint);
    if (incomingDir.dot(normal) > 0) normal = -normal;

    G4double incidentRad = incomingDir.angle(-normal);
    G4double incidentDeg = incidentRad / deg;
    if (incidentDeg > 90.0) incidentDeg = 90.0;
    if (incidentDeg < 0) incidentDeg = 0;

    // --- 7. Sample outgoing angle (in-plane) ---
    G4double thetaOutDeg = reflector->SampleOutgoingAngle(incidentDeg);
    G4double thetaOutRad = thetaOutDeg * deg;

    if (!gPrintLimitReached && gPrintCount < gMaxPrint) {
        gPrintCount++;
        gReflectionCount++;
        printf("  [REFLECTION %d] incident=%.2f deg, reflected=%.2f deg [AZIMUTH=%s]\n",
               gReflectionCount, incidentDeg, thetaOutDeg,
               GetAzimuthalModelName());
        fflush(stdout);
        if (gPrintCount >= gMaxPrint) {
            gPrintLimitReached = true;
            printf("  ... (further reflections suppressed)\n");
            fflush(stdout);
        }
    }

    G4ThreeVector norm = normal.unit();
    G4ThreeVector perp = incomingDir - (incomingDir.dot(norm)) * norm;
    if (perp.mag() < 1e-10) {
        perp = G4ThreeVector(1, 0, 0);
        if (std::abs(norm.dot(perp)) > 0.9999) perp = G4ThreeVector(0, 1, 0);
    }
    perp = perp.unit();
    G4ThreeVector perp2 = norm.cross(perp).unit();

    
    G4double phi = SampleAzimuthalAngle();

    G4ThreeVector tangentialDir = std::cos(phi) * perp + std::sin(phi) * perp2;
    G4ThreeVector finalDir = std::cos(thetaOutRad) * norm + std::sin(thetaOutRad) * tangentialDir;
    finalDir = finalDir.unit();

    G4double dotNormal = finalDir.dot(norm);
    if (dotNormal < 0) {
        finalDir = finalDir - 2.0 * dotNormal * norm;
        finalDir = finalDir.unit();
    }
    if (std::abs(finalDir.dot(norm)) < 1e-6) {
        finalDir = finalDir + 1e-6 * norm;
        finalDir = finalDir.unit();
    }

    reflector->RecordReflection(incidentDeg, thetaOutDeg,
                                finalDir.x(), finalDir.y(), finalDir.z(),
                                norm.x(), norm.y(), norm.z());

    particleChange->ProposeMomentumDirection(finalDir);

    return particleChange;
}

