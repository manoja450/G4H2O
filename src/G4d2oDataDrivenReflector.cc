#include "G4d2oDataDrivenReflector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <random>

// ============================================================
// Static member initialization
// ============================================================

G4d2oDataDrivenReflector* G4d2oDataDrivenReflector::fInstance = nullptr;
TTree* G4d2oDataDrivenReflector::fReflectionTree = nullptr;
G4double G4d2oDataDrivenReflector::fRecordedIncident = 0;
G4double G4d2oDataDrivenReflector::fRecordedReflected = 0;
G4double G4d2oDataDrivenReflector::fRecordedPx = 0;
G4double G4d2oDataDrivenReflector::fRecordedPy = 0;
G4double G4d2oDataDrivenReflector::fRecordedPz = 0;
G4double G4d2oDataDrivenReflector::fRecordedNx = 0;
G4double G4d2oDataDrivenReflector::fRecordedNy = 0;
G4double G4d2oDataDrivenReflector::fRecordedNz = 0;
G4int G4d2oDataDrivenReflector::fRecordedEventID = 0;
G4int G4d2oDataDrivenReflector::fReflectionCounter = 0;
G4int G4d2oDataDrivenReflector::fCurrentEventNumber = 0;

// ============================================================
// Constructor / Destructor
// ============================================================

G4d2oDataDrivenReflector::G4d2oDataDrivenReflector(G4double reflectivity)
    : fReflectivity(reflectivity), fRNG(std::random_device{}()), fRandDist(0.0, 1.0) {
    G4cout << "\n=========================================================" << G4endl;
    G4cout << "G4d2oDataDrivenReflector: Created" << G4endl;
    G4cout << "  Reflectivity = " << reflectivity << G4endl;
    G4cout << "  Function A: Incident angle interpolation: ENABLED" << G4endl;
    G4cout << "  Function B: Continuous PDF interpolation: ENABLED" << G4endl;
    G4cout << "  Function C: CDF sampling: ENABLED" << G4endl;
    G4cout << "=========================================================\n" << G4endl;
}

G4d2oDataDrivenReflector::~G4d2oDataDrivenReflector() {
    G4cout << "G4d2oDataDrivenReflector: Destroyed" << G4endl;
}

// ============================================================
// Singleton Management
// ============================================================

G4d2oDataDrivenReflector* G4d2oDataDrivenReflector::GetInstance() {
    return fInstance;
}

void G4d2oDataDrivenReflector::CreateInstance(G4double reflectivity, const G4String& dataDirectory) {
    if (fInstance == nullptr) {
        fInstance = new G4d2oDataDrivenReflector(reflectivity);
        fInstance->LoadAllPDFs(dataDirectory);

        // Run validation
        G4cout << "\n>>> RUNNING VALIDATION <<<" << G4endl;
        for (auto& pair : fInstance->fPDFs) {
            fInstance->ValidatePDF(pair.first, 50000);
        }
    }
}

void G4d2oDataDrivenReflector::DeleteInstance() {
    if (fInstance) {
        delete fInstance;
        fInstance = nullptr;
    }
}

// ============================================================
// FUNCTION B: Continuous Interpolation Helpers
// ============================================================

G4double G4d2oDataDrivenReflector::LinearInterpolate(const std::vector<G4double>& x,
                                                      const std::vector<G4double>& y,
                                                      G4double x0) const {
    if (x.empty() || y.empty()) return 0.0;
    if (x0 <= x.front()) return y.front();
    if (x0 >= x.back()) return y.back();

    for (size_t i = 0; i < x.size() - 1; ++i) {
        if (x0 >= x[i] && x0 <= x[i+1]) {
            G4double t = (x0 - x[i]) / (x[i+1] - x[i]);
            return y[i] * (1.0 - t) + y[i+1] * t;
        }
    }
    return y.back();
}

std::vector<G4double> G4d2oDataDrivenReflector::CreateContinuousPDF(
    const std::vector<G4double>& theta,
    const std::vector<G4double>& intensity,
    G4int nFineBins) const {

    std::vector<G4double> fineIntensity;
    G4double thetaMin = -90.0;
    G4double thetaMax = 90.0;
    G4double step = (thetaMax - thetaMin) / nFineBins;

    for (G4int i = 0; i <= nFineBins; ++i) {
        G4double t = thetaMin + i * step;
        G4double interpIntensity = LinearInterpolate(theta, intensity, t);
        fineIntensity.push_back(interpIntensity);
    }

    return fineIntensity;
}

// ============================================================
// PDF Loading with Function B (Continuous Interpolation)
// ============================================================

void G4d2oDataDrivenReflector::NormalizeAndBuildCDF(PDF& pdf) const {
    if (pdf.intensity.empty()) return;

    // Normalize
    G4double sum = 0.0;
    for (G4double val : pdf.intensity) sum += val;
    if (sum > 0) {
        for (G4double& val : pdf.intensity) val /= sum;
    }

    // Build CDF
    pdf.cdf.resize(pdf.intensity.size());
    G4double cumulative = 0.0;
    for (size_t i = 0; i < pdf.intensity.size(); ++i) {
        cumulative += pdf.intensity[i];
        pdf.cdf[i] = cumulative;
    }

    pdf.maxIntensity = *std::max_element(pdf.intensity.begin(), pdf.intensity.end());
}

void G4d2oDataDrivenReflector::LoadPDF(G4double incidentAngleDeg, const G4String& filename) {
    PDF pdf;

    // Read raw data from thesis (coarse grid, ~5deg resolution)
    std::vector<G4double> coarseTheta, coarseIntensity;
    std::ifstream file(filename);
    if (!file.is_open()) {
        G4cerr << "ERROR: Cannot open file " << filename << G4endl;
        return;
    }

    G4double theta, intensity;
    while (file >> theta >> intensity) {
        coarseTheta.push_back(theta);
        coarseIntensity.push_back(intensity);
    }
    file.close();

    if (coarseTheta.empty()) {
        G4cerr << "ERROR: No data in file " << filename << G4endl;
        return;
    }

    G4cout << "  Read " << coarseTheta.size() << " coarse points from " << filename << G4endl;

    // ============================================================
    // FUNCTION B: Interpolate to continuous fine grid (0.5deg resolution)
    // ============================================================
    const G4int nFineBins = 360;  // -90 to 90 in 0.5deg steps
    pdf.theta.reserve(nFineBins + 1);
    pdf.intensity.reserve(nFineBins + 1);

    G4double thetaMin = -90.0;
    G4double thetaMax = 90.0;
    G4double step = (thetaMax - thetaMin) / nFineBins;

    for (G4int i = 0; i <= nFineBins; ++i) {
        G4double t = thetaMin + i * step;
        pdf.theta.push_back(t);

        G4double interpIntensity = LinearInterpolate(coarseTheta, coarseIntensity, t);
        pdf.intensity.push_back(interpIntensity);
    }

    NormalizeAndBuildCDF(pdf);
    fPDFs[static_cast<G4int>(incidentAngleDeg)] = pdf;

    G4cout << "  Loaded PDF for incident " << incidentAngleDeg << "deg -> interpolated to "
           << pdf.theta.size() << " continuous points" << G4endl;
}

void G4d2oDataDrivenReflector::LoadAllPDFs(const G4String& directory) {
    std::vector<G4int> angles = {0, 10, 20, 30, 40, 50, 60, 70, 80};

    G4int loaded = 0;
    for (G4int angle : angles) {
        std::stringstream ss;
        ss << directory << "/incident_" << angle << "deg.txt";
        G4String filename = ss.str();

        std::ifstream test(filename);
        if (test.good()) {
            test.close();
            LoadPDF(angle, filename);
            loaded++;
        } else {
            G4cerr << "WARNING: File not found: " << filename << G4endl;
        }
    }

    G4cout << "\nLoaded " << loaded << " PDF files from " << directory << G4endl;
    PrintLoadedAngles();
}

void G4d2oDataDrivenReflector::PrintLoadedAngles() const {
    G4cout << "\n========== PDF LIBRARY SUMMARY ==========" << G4endl;
    G4cout << "Function A: Incident angles available:" << G4endl;
    for (auto& pair : fPDFs) {
        G4cout << "  " << pair.first << "deg: " << pair.second.theta.size()
               << " continuous points (0.5deg resolution)" << G4endl;
    }
    G4cout << "Function B: Continuous interpolation applied" << G4endl;
    G4cout << "Function C: CDF sampling ready" << G4endl;
    G4cout << "==========================================\n" << G4endl;
}

// ============================================================
// FUNCTION A: Find Bounding PDFs for Interpolation
// ============================================================

void G4d2oDataDrivenReflector::FindBoundingPDFs(G4double incidentAngleDeg,
                                                  const PDF*& pdfLow,
                                                  const PDF*& pdfHigh,
                                                  G4double& weightLow,
                                                  G4double& weightHigh) const {
    pdfLow = nullptr;
    pdfHigh = nullptr;

    if (fPDFs.empty()) return;

    // Find the two closest incident angles
    G4int lowerAngle = -1, upperAngle = -1;

    for (auto& pair : fPDFs) {
        G4int angle = pair.first;
        if (angle <= incidentAngleDeg) {
            lowerAngle = angle;
        }
        if (angle >= incidentAngleDeg && upperAngle == -1) {
            upperAngle = angle;
        }
    }

    if (lowerAngle == -1) {
        // Use smallest angle
        lowerAngle = upperAngle = fPDFs.begin()->first;
        weightLow = weightHigh = 0.5;
    } else if (upperAngle == -1) {
        // Use largest angle
        upperAngle = lowerAngle = fPDFs.rbegin()->first;
        weightLow = weightHigh = 0.5;
    } else if (lowerAngle == upperAngle) {
        weightLow = weightHigh = 0.5;
    } else {
        weightHigh = (incidentAngleDeg - lowerAngle) / (upperAngle - lowerAngle);
        weightLow = 1.0 - weightHigh;
    }

    auto itLow = fPDFs.find(lowerAngle);
    auto itHigh = fPDFs.find(upperAngle);

    if (itLow != fPDFs.end()) pdfLow = &(itLow->second);
    if (itHigh != fPDFs.end()) pdfHigh = &(itHigh->second);
}

// ============================================================
// FUNCTION C: Sample from PDF
// ============================================================

G4double G4d2oDataDrivenReflector::SampleFromPDF(const PDF& pdf, G4double rand) const {
    if (pdf.cdf.empty()) return 0.0;

    // Find bin
    size_t idx = 0;
    for (size_t i = 0; i < pdf.cdf.size(); ++i) {
        if (rand <= pdf.cdf[i]) {
            idx = i;
            break;
        }
    }

    // Linear interpolation within bin for smooth sampling
    if (idx == 0) return pdf.theta[0];

    G4double cdf_prev = (idx > 0) ? pdf.cdf[idx-1] : 0.0;
    G4double frac = (rand - cdf_prev) / (pdf.cdf[idx] - cdf_prev);

    return pdf.theta[idx-1] + frac * (pdf.theta[idx] - pdf.theta[idx-1]);
}

G4double G4d2oDataDrivenReflector::InterpolateThetaOut(G4double incidentAngleDeg, G4double rand) const {
    const PDF* pdfLow = nullptr;
    const PDF* pdfHigh = nullptr;
    G4double weightLow = 0.5, weightHigh = 0.5;

    FindBoundingPDFs(incidentAngleDeg, pdfLow, pdfHigh, weightLow, weightHigh);

    if (pdfLow == nullptr && pdfHigh == nullptr) {
        // Fallback to Lambertian
        G4double u = fRandDist(fRNG);
        return std::acos(std::sqrt(u)) * 180.0 / M_PI;
    }

    if (pdfLow == pdfHigh || pdfHigh == nullptr) {
        return SampleFromPDF(*pdfLow, rand);
    }

    if (pdfLow == nullptr) {
        return SampleFromPDF(*pdfHigh, rand);
    }

    // ============================================================
    // FUNCTION A: Weighted average of two PDFs
    // ============================================================
    G4double thetaLow = SampleFromPDF(*pdfLow, rand);
    G4double thetaHigh = SampleFromPDF(*pdfHigh, rand);

    return weightLow * thetaLow + weightHigh * thetaHigh;
}

G4double G4d2oDataDrivenReflector::GetThetaOut(G4double incidentAngleDeg) const {
    G4double rand = fRandDist(fRNG);
    return InterpolateThetaOut(incidentAngleDeg, rand);
}

// ============================================================
// FUNCTION C: Sample Outgoing Angle
// ============================================================
G4double G4d2oDataDrivenReflector::SampleOutgoingAngle(G4double incidentAngleDeg) const {
    return GetThetaOut(incidentAngleDeg);
}

// ============================================================
// Reflection Direction (Main interface for Geant4)
// NOTE: not currently called from the active PostStepDoIt path (which
// builds finalDir itself in G4d2oCustomOpBoundary.cc) - kept here for
// API compatibility, with the matching fixes applied.
// ============================================================

G4ThreeVector G4d2oDataDrivenReflector::GetReflectionDirection(const G4ThreeVector& incomingDirection,
                                                                 const G4ThreeVector& normal,
                                                                 G4double incidentAngleRad) const {
    G4double incidentDeg = incidentAngleRad * 180.0 / M_PI;

    // Sample outgoing angle from thesis data
    G4double thetaOutDeg = SampleOutgoingAngle(incidentDeg);
    G4double thetaOutRad = thetaOutDeg * M_PI / 180.0;

    // Construct orthonormal basis
    G4ThreeVector norm = normal.unit();
    G4ThreeVector incoming = incomingDirection.unit();

    // Tangent vector (perpendicular to normal, in plane of incidence)
    G4ThreeVector tangent = incoming - (incoming.dot(norm)) * norm;
    if (tangent.mag() < 1e-10) {
        tangent = G4ThreeVector(1, 0, 0);
        if (std::abs(norm.dot(tangent)) > 0.9999) tangent = G4ThreeVector(0, 1, 0);
    }
    tangent = tangent.unit();

    // Bitangent
    G4ThreeVector bitangent = norm.cross(tangent);

    // Random azimuth (isotropic out-of-plane), built as a true rotation
    // about norm so the normal component stays cos(thetaOut) for every
    // azimuth (see G4d2oCustomOpBoundary.cc step 10 for the same fix).
    G4double azimuth = 2.0 * M_PI * fRandDist(fRNG);
    G4ThreeVector tangentialDir = std::cos(azimuth) * tangent + std::sin(azimuth) * bitangent;
    G4ThreeVector reflectedDir = std::cos(thetaOutRad) * norm + std::sin(thetaOutRad) * tangentialDir;

    return reflectedDir.unit();
}

// ============================================================
// VALIDATION METHOD
// ============================================================

void G4d2oDataDrivenReflector::ValidatePDF(G4double incidentAngleDeg, G4int nSamples) const {
    G4cout << "\n========== VALIDATING incident " << incidentAngleDeg << "deg ==========" << G4endl;

    // Find the PDF for this incident angle
    auto it = fPDFs.find(static_cast<G4int>(incidentAngleDeg));
    if (it == fPDFs.end()) {
        G4cerr << "ERROR: No PDF loaded for " << incidentAngleDeg << "deg" << G4endl;
        return;
    }

    const PDF& refPDF = it->second;
    G4cout << "Reference PDF has " << refPDF.theta.size() << " continuous bins" << G4endl;

    // Sample nSamples angles
    std::vector<G4double> samples;
    samples.reserve(nSamples);
    for (G4int i = 0; i < nSamples; ++i) {
        samples.push_back(SampleOutgoingAngle(incidentAngleDeg));
    }

    // Create histogram (2.5deg bins for comparison)
    const G4int nBins = 72;
    G4double binMin = -90.0, binMax = 90.0;
    G4double binWidth = (binMax - binMin) / nBins;
    std::vector<G4double> histCounts(nBins, 0.0);

    for (G4double theta : samples) {
        int bin = static_cast<int>((theta - binMin) / binWidth);
        if (bin >= 0 && bin < nBins) histCounts[bin]++;
    }

    // Normalize to probability density
    std::vector<G4double> histPDF(nBins, 0.0);
    for (int i = 0; i < nBins; ++i) {
        histPDF[i] = histCounts[i] / (nSamples * binWidth);
    }

    // Compute shape correlation with reference PDF (binned to same bins)
    std::vector<G4double> refBinned(nBins, 0.0);
    for (int i = 0; i < nBins; ++i) {
        G4double binCenter = binMin + (i + 0.5) * binWidth;
        refBinned[i] = LinearInterpolate(refPDF.theta, refPDF.intensity, binCenter);
    }

    // Calculate correlation
    G4double mean1 = 0.0, mean2 = 0.0;
    for (int i = 0; i < nBins; ++i) {
        mean1 += refBinned[i];
        mean2 += histPDF[i];
    }
    mean1 /= nBins;
    mean2 /= nBins;

    G4double cov = 0.0, var1 = 0.0, var2 = 0.0;
    for (int i = 0; i < nBins; ++i) {
        cov += (refBinned[i] - mean1) * (histPDF[i] - mean2);
        var1 += (refBinned[i] - mean1) * (refBinned[i] - mean1);
        var2 += (histPDF[i] - mean2) * (histPDF[i] - mean2);
    }

    G4double shapeCorr = (var1 * var2 > 0) ? cov / std::sqrt(var1 * var2) : 0.0;

    G4cout << "\nVALIDATION RESULTS:" << G4endl;
    G4cout << "  Shape correlation = " << shapeCorr << G4endl;

    if (shapeCorr > 0.99) {
        G4cout << "  EXCELLENT - Functions A, B, C working correctly" << G4endl;
    } else if (shapeCorr > 0.95) {
        G4cout << "  GOOD - Acceptable, minor deviations" << G4endl;
    } else {
        G4cout << "  POOR - Check PDF loading or digitization" << G4endl;
    }
    G4cout << "=========================================================\n" << G4endl;
}

// ============================================================
// ROOT Recording Methods
// ============================================================

void G4d2oDataDrivenReflector::SetRootFile(TFile* file) {
    if (fReflectionTree == nullptr && file != nullptr) {
        file->cd();
        fReflectionTree = new TTree("ReflectionTree", "Reflection Data");
        fReflectionTree->Branch("incident_deg", &fRecordedIncident, "incident_deg/D");
        fReflectionTree->Branch("reflected_deg", &fRecordedReflected, "reflected_deg/D");
        // Global direction branches (diagnostic)
        fReflectionTree->Branch("px_global", &fRecordedPx, "px_global/D");
        fReflectionTree->Branch("py_global", &fRecordedPy, "py_global/D");
        fReflectionTree->Branch("pz_global", &fRecordedPz, "pz_global/D");
        // Surface normal branches (diagnostic)
        fReflectionTree->Branch("nx_global", &fRecordedNx, "nx_global/D");
        fReflectionTree->Branch("ny_global", &fRecordedNy, "ny_global/D");
        fReflectionTree->Branch("nz_global", &fRecordedNz, "nz_global/D");
        fReflectionTree->Branch("event_id", &fRecordedEventID, "event_id/I");
        G4cout << "Created ReflectionTree with global direction branches" << G4endl;
    }
}

void G4d2oDataDrivenReflector::CloseRootTree() {
    if (fReflectionTree) {
        G4cout << "ReflectionTree: Recorded " << fReflectionCounter << " reflections" << G4endl;
        fReflectionTree->AutoSave();
        fReflectionTree = nullptr;
        fReflectionCounter = 0;
    }
}

void G4d2oDataDrivenReflector::RecordReflection(G4double incidentAngleDeg, G4double reflectedAngleDeg,
                                                 G4double px, G4double py, G4double pz,
                                                 G4double nx, G4double ny, G4double nz) {
    fRecordedIncident = incidentAngleDeg;
    fRecordedReflected = reflectedAngleDeg;
    fRecordedPx = px;
    fRecordedPy = py;
    fRecordedPz = pz;
    fRecordedNx = nx;
    fRecordedNy = ny;
    fRecordedNz = nz;
    fRecordedEventID = fCurrentEventNumber;

    if (fReflectionTree) {
        fReflectionTree->Fill();
        fReflectionCounter++;
    }
}

void G4d2oDataDrivenReflector::SetCurrentEventNumber(G4int eventNum) {
    fCurrentEventNumber = eventNum;
}
