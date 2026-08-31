#include "G4d2oPhotonDiagnostics.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <iomanip>
#include <fstream>

G4d2oPhotonDiagnostics* G4d2oPhotonDiagnostics::fInstance = nullptr;

G4d2oPhotonDiagnostics* G4d2oPhotonDiagnostics::GetInstance() {
    if (fInstance == nullptr) CreateInstance();
    return fInstance;
}

void G4d2oPhotonDiagnostics::CreateInstance() {
    if (fInstance == nullptr) fInstance = new G4d2oPhotonDiagnostics();
}

void G4d2oPhotonDiagnostics::DeleteInstance() {
    if (fInstance) { delete fInstance; fInstance = nullptr; }
}

G4d2oPhotonDiagnostics::G4d2oPhotonDiagnostics()
    : fCurrentEventID(-1), fCurrentPhotonID(-1),
      fTotalEvents(0), fTotalPhotons(0), fTotalReflections(0),
      fTotalPMTHits(0), fTotalAbsorbed(0), fTotalBadReflections(0),
      fTotalQualityPass(0), fPrintTraces(true), fMaxTraces(20),
      fTraceEventID(-1), fTracesPrinted(0) {
    fTotalCherenkovPhotons = 0;
    fTotalTyvekReflections = 0;
    fTotalPMTReach = 0;
    fTotalPhotoelectrons = 0;
    for (int i = 0; i <= kOther; ++i) fTerminationCounts[(TerminationReason)i] = 0;

    G4cout << "\n=========================================================" << G4endl;
    G4cout << "G4d2oPhotonDiagnostics: Created" << G4endl;
    G4cout << "=========================================================\n" << G4endl;
}

G4d2oPhotonDiagnostics::~G4d2oPhotonDiagnostics() {
    PrintGlobalSummary();
    PrintChainSummary();
    PrintTerminationSummary();
    WriteDiagnostics();
}

void G4d2oPhotonDiagnostics::Reset() {
    fPhotonTraces.clear();
    fEventSummaries.clear();
    fEventChains.clear();
    fPhotonReflectionCounts.clear();
    fPhotonTermination.clear();
    fTerminationCounts.clear();
    for (int i = 0; i <= kOther; ++i) fTerminationCounts[(TerminationReason)i] = 0;
    fTotalEvents = 0;
    fTotalPhotons = 0;
    fTotalReflections = 0;
    fTotalPMTHits = 0;
    fTotalAbsorbed = 0;
    fTotalBadReflections = 0;
    fTotalQualityPass = 0;
    fTotalCherenkovPhotons = 0;
    fTotalTyvekReflections = 0;
    fTotalPMTReach = 0;
    fTotalPhotoelectrons = 0;
    fTracesPrinted = 0;
}

void G4d2oPhotonDiagnostics::NewEvent(G4int eventID) {
    fCurrentEventID = eventID;
    fTotalEvents++;
    ClearAllPhotonCounts();
    EventSummary summary; summary.eventID = eventID; summary.nPhotons = 0; summary.nPMTHits = 0; summary.nAbsorbed = 0; summary.totalPE = 0; summary.qualityCutPass = 0; summary.meanReflections = 0; summary.maxReflections = 0; summary.fractionBadReflections = 0;
    fEventSummaries[eventID] = summary;
    EventChain chain; chain.eventID = eventID; chain.nCherenkov = 0; chain.nTyvekReflections = 0; chain.nPMTReach = 0; chain.nPE = 0;
    fEventChains[eventID] = chain;
}

void G4d2oPhotonDiagnostics::RecordCherenkovPhoton(G4double energy, const G4ThreeVector& position) {
    fTotalCherenkovPhotons++;
    if (fCurrentEventID >= 0) fEventChains[fCurrentEventID].nCherenkov++;
}

void G4d2oPhotonDiagnostics::RecordTyvekReflection() {
    fTotalTyvekReflections++;
    if (fCurrentEventID >= 0) fEventChains[fCurrentEventID].nTyvekReflections++;
}

void G4d2oPhotonDiagnostics::RecordPMTReach() {
    fTotalPMTReach++;
    if (fCurrentEventID >= 0) fEventChains[fCurrentEventID].nPMTReach++;
}

void G4d2oPhotonDiagnostics::RecordPhotoelectron() {
    fTotalPhotoelectrons++;
    if (fCurrentEventID >= 0) fEventChains[fCurrentEventID].nPE++;
}

void G4d2oPhotonDiagnostics::RecordWaterAbsorption(const G4ThreeVector& position) {
    fTotalAbsorbed++;
}

void G4d2oPhotonDiagnostics::RecordEventChain(G4int eventID, G4int nCherenkov,
                                              G4int nTyvekReflections, G4int nPMTReach, G4int nPE) {
    if (eventID >= 0) {
        fEventChains[eventID].nCherenkov = nCherenkov;
        fEventChains[eventID].nTyvekReflections = nTyvekReflections;
        fEventChains[eventID].nPMTReach = nPMTReach;
        fEventChains[eventID].nPE = nPE;
    }
}

G4int G4d2oPhotonDiagnostics::GetReflectionCountForPhoton(G4int photonID) const {
    auto it = fPhotonReflectionCounts.find(photonID);
    return (it != fPhotonReflectionCounts.end()) ? it->second : 0;
}

void G4d2oPhotonDiagnostics::IncrementReflectionCountForPhoton(G4int photonID) {
    fPhotonReflectionCounts[photonID]++;
}

void G4d2oPhotonDiagnostics::ResetReflectionCountForPhoton(G4int photonID) {
    fPhotonReflectionCounts[photonID] = 0;
}

void G4d2oPhotonDiagnostics::ClearAllPhotonCounts() {
    fPhotonReflectionCounts.clear();
}

void G4d2oPhotonDiagnostics::RecordTermination(G4int photonID, TerminationReason reason) {
    fPhotonTermination[photonID] = reason;
    fTerminationCounts[reason]++;
    if (photonID == fCurrentPhotonID) {
        fCurrentTrace.terminationReason = reason;
    }
}

void G4d2oPhotonDiagnostics::PrintTerminationSummary() const {
    G4cout << "\n" << std::string(70, '=') << G4endl;
    G4cout << "PHOTON TERMINATION REASONS" << G4endl;
    G4cout << std::string(70, '=') << G4endl;
    G4int total = 0;
    for (auto& kv : fTerminationCounts) total += kv.second;
    if (total == 0) {
        G4cout << "  No termination data recorded." << G4endl;
        return;
    }
    const char* names[] = {"Reached PMT", "Water Absorbed", "Tyvek Absorbed", "Max Reflections", "Boundary Kill", "Other"};
    for (int i = 0; i <= kOther; ++i) {
        G4int count = fTerminationCounts.at((TerminationReason)i);
        G4cout << "  " << names[i] << ": " << count << " (" << (G4double)count/total*100 << "%)" << G4endl;
    }
    G4cout << std::string(70, '=') << "\n" << G4endl;
}

void G4d2oPhotonDiagnostics::StartPhoton(G4int eventID, G4int photonID, G4double energy,
                                          const G4ThreeVector& position, const G4ThreeVector& direction) {
    fCurrentPhotonID = photonID;
    fTotalPhotons++;
    ResetReflectionCountForPhoton(photonID);
    fCurrentTrace.eventID = eventID;
    fCurrentTrace.photonID = photonID;
    fCurrentTrace.initialEnergy = energy;
    fCurrentTrace.startPos = position;
    fCurrentTrace.startDir = direction;
    fCurrentTrace.nReflections = 0;
    fCurrentTrace.reachedPMT = false;
    fCurrentTrace.finalTime = 0;
    fCurrentTrace.terminationReason = kOther;
    fCurrentTrace.reflectionAngles.clear();
    fCurrentTrace.incidentAngles.clear();
    fCurrentTrace.reflectionPositions.clear();
    fCurrentTrace.reflectionOutgoing.clear();
    fCurrentTrace.globalPzValues.clear();
    fCurrentTrace.dotProducts.clear();
    RecordCherenkovPhoton(energy, position);
}

void G4d2oPhotonDiagnostics::RecordReflection(G4int reflectionNumber,
                                               const G4ThreeVector& position,
                                               const G4ThreeVector& normal,
                                               const G4ThreeVector& incomingDir,
                                               const G4ThreeVector& outgoingDir,
                                               G4double incidentAngleDeg,
                                               G4double reflectedAngleDeg,
                                               G4double globalPz,
                                               G4double dotProduct) {
    fTotalReflections++;
    fCurrentTrace.nReflections++;
    fCurrentTrace.reflectionAngles.push_back(reflectedAngleDeg);
    fCurrentTrace.incidentAngles.push_back(incidentAngleDeg);
    fCurrentTrace.reflectionPositions.push_back(position);
    fCurrentTrace.reflectionOutgoing.push_back(outgoingDir);
    fCurrentTrace.globalPzValues.push_back(globalPz);
    fCurrentTrace.dotProducts.push_back(dotProduct);
    if (dotProduct > 0) fTotalBadReflections++;
    RecordTyvekReflection();
    IncrementReflectionCountForPhoton(fCurrentPhotonID);

    bool printThis = false;
    if (fPrintTraces && fTracesPrinted < fMaxTraces) printThis = true;
    if (fTraceEventID >= 0 && fCurrentEventID == fTraceEventID) printThis = true;
    if (printThis) {
        G4cout << "\n  Photon " << fCurrentPhotonID << ", Reflection " << reflectionNumber << G4endl;
        G4cout << "    Position:   (" << std::setw(8) << position.x()/mm
               << ", " << std::setw(8) << position.y()/mm
               << ", " << std::setw(8) << position.z()/mm << ") mm" << G4endl;
        G4cout << "    Normal:     (" << std::setw(8) << normal.x()
               << ", " << std::setw(8) << normal.y()
               << ", " << std::setw(8) << normal.z() << ")" << G4endl;
        G4cout << "    Incoming:   (" << std::setw(8) << incomingDir.x()
               << ", " << std::setw(8) << incomingDir.y()
               << ", " << std::setw(8) << incomingDir.z() << ")" << G4endl;
        G4cout << "    Outgoing:   (" << std::setw(8) << outgoingDir.x()
               << ", " << std::setw(8) << outgoingDir.y()
               << ", " << std::setw(8) << outgoingDir.z() << ")" << G4endl;
        G4cout << "    Incident:   " << incidentAngleDeg << "°" << G4endl;
        G4cout << "    Reflected:  " << reflectedAngleDeg << "°" << G4endl;
        G4cout << "    dot(out,n): " << dotProduct << " ";
        if (dotProduct < 0) G4cout << "✅ into WATER" << G4endl;
        else G4cout << "❌ into TYVEK (BAD!)" << G4endl;
        if (fTraceEventID < 0 || fCurrentEventID != fTraceEventID) fTracesPrinted++;
    }
}

void G4d2oPhotonDiagnostics::RecordAbsorption(const G4ThreeVector& position, const G4String& reason) {
    fTotalAbsorbed++;
    RecordWaterAbsorption(position);
    RecordTermination(fCurrentPhotonID, kWaterAbsorbed);
}

void G4d2oPhotonDiagnostics::RecordPMTHit(const G4ThreeVector& position, G4double time) {
    fTotalPMTHits++;
    fCurrentTrace.reachedPMT = true;
    fCurrentTrace.finalTime = time;
    RecordPMTReach();
    RecordTermination(fCurrentPhotonID, kReachedPMT);
}

void G4d2oPhotonDiagnostics::EndPhoton(G4bool reachedPMT, G4int totalReflections) {
    fCurrentTrace.reachedPMT = reachedPMT;
    fCurrentTrace.nReflections = totalReflections;
    fPhotonTraces[fCurrentEventID][fCurrentPhotonID] = fCurrentTrace;
    if (!reachedPMT && fPhotonTermination.find(fCurrentPhotonID) == fPhotonTermination.end()) {
        RecordTermination(fCurrentPhotonID, kOther);
    }
}

void G4d2oPhotonDiagnostics::RecordEventResult(G4int eventID, G4int nPhotons, G4int nPMTHits,
                                                G4int nAbsorbed, G4double totalPE,
                                                G4int qualityCutPass,
                                                G4double meanReflections, G4int maxReflections) {
    if (eventID >= 0) {
        fEventSummaries[eventID].nPhotons = nPhotons;
        fEventSummaries[eventID].nPMTHits = nPMTHits;
        fEventSummaries[eventID].nAbsorbed = nAbsorbed;
        fEventSummaries[eventID].totalPE = totalPE;
        fEventSummaries[eventID].qualityCutPass = qualityCutPass;
        fEventSummaries[eventID].meanReflections = meanReflections;
        fEventSummaries[eventID].maxReflections = maxReflections;
        int totalRefs = 0, badRefs = 0;
        if (fPhotonTraces.find(eventID) != fPhotonTraces.end()) {
            for (auto& pair : fPhotonTraces[eventID]) {
                for (double dot : pair.second.dotProducts) {
                    totalRefs++;
                    if (dot > 0) badRefs++;
                }
            }
        }
        if (totalRefs > 0) fEventSummaries[eventID].fractionBadReflections = (G4double)badRefs / totalRefs;
        if (qualityCutPass) fTotalQualityPass++;
    }

    if (fTotalEvents % 100 == 0 && fTotalEvents > 0) {
        G4cout << "\n" << std::string(70, '=') << G4endl;
        G4cout << "PHOTON SURVIVAL CHAIN - After " << fTotalEvents << " events" << G4endl;
        G4cout << std::string(70, '=') << G4endl;
        PrintChainSummary();
        G4cout << std::string(70, '=') << "\n" << G4endl;
    }
}

void G4d2oPhotonDiagnostics::PrintPhotonTrace(G4int eventID, G4int photonID) const {
    auto itEvent = fPhotonTraces.find(eventID);
    if (itEvent == fPhotonTraces.end()) { G4cout << "Event not found" << G4endl; return; }
    auto itPhoton = itEvent->second.find(photonID);
    if (itPhoton == itEvent->second.end()) { G4cout << "Photon not found" << G4endl; return; }
    const PhotonTrace& trace = itPhoton->second;
    G4cout << "\n=========================================================" << G4endl;
    G4cout << "PHOTON TRACE: Event " << eventID << ", Photon " << photonID << G4endl;
    G4cout << "=========================================================" << G4endl;
    G4cout << "Initial energy: " << trace.initialEnergy/eV << " eV" << G4endl;
    G4cout << "Start position: (" << trace.startPos.x()/mm << ", "
           << trace.startPos.y()/mm << ", " << trace.startPos.z()/mm << ") mm" << G4endl;
    G4cout << "Total reflections: " << trace.nReflections << G4endl;
    G4cout << "Reached PMT: " << (trace.reachedPMT ? "YES" : "NO") << G4endl;
    const char* names[] = {"Reached PMT", "Water Absorbed", "Tyvek Absorbed", "Max Reflections", "Boundary Kill", "Other"};
    G4cout << "Termination reason: " << names[trace.terminationReason] << G4endl;
    if (trace.nReflections > 0) {
        G4cout << "\nReflection details:" << G4endl;
        G4cout << "  # | Incident | Reflected | dot(out,n) | Status" << G4endl;
        G4cout << "----|----------|-----------|------------|--------" << G4endl;
        for (size_t i = 0; i < trace.reflectionAngles.size(); i++) {
            G4cout << "  " << std::setw(2) << i+1 << " | "
                   << std::setw(8) << std::fixed << std::setprecision(1)
                   << trace.incidentAngles[i] << " | "
                   << std::setw(9) << trace.reflectionAngles[i] << " | "
                   << std::setw(10) << std::fixed << std::setprecision(4)
                   << trace.dotProducts[i] << " | ";
            if (trace.dotProducts[i] < 0) G4cout << "water ✅" << G4endl;
            else G4cout << "TYVEK ❌" << G4endl;
        }
    }
    G4cout << "=========================================================\n" << G4endl;
}

void G4d2oPhotonDiagnostics::PrintEventSummary(G4int eventID) const {
    auto it = fEventSummaries.find(eventID);
    if (it == fEventSummaries.end()) { G4cout << "Event not found" << G4endl; return; }
    const EventSummary& summary = it->second;
    G4cout << "\n=========================================================" << G4endl;
    G4cout << "EVENT SUMMARY: Event " << eventID << G4endl;
    G4cout << "=========================================================" << G4endl;
    G4cout << "  Total photons:        " << summary.nPhotons << G4endl;
    G4cout << "  PMT hits:             " << summary.nPMTHits << G4endl;
    G4cout << "  Absorbed:             " << summary.nAbsorbed << G4endl;
    G4cout << "  Total PE:             " << summary.totalPE << G4endl;
    G4cout << "  Quality cut pass:     " << (summary.qualityCutPass ? "YES" : "NO") << G4endl;
    G4cout << "  Mean reflections:     " << summary.meanReflections << G4endl;
    G4cout << "  Max reflections:      " << summary.maxReflections << G4endl;
    G4cout << "  Bad reflections (%):  " << summary.fractionBadReflections * 100 << "%" << G4endl;
    G4cout << "=========================================================\n" << G4endl;
}

void G4d2oPhotonDiagnostics::PrintBadReflections() const {
    G4cout << "\n=========================================================" << G4endl;
    G4cout << "BAD REFLECTIONS SUMMARY" << G4endl;
    G4cout << "=========================================================" << G4endl;
    if (fTotalBadReflections == 0) {
        G4cout << "✅ NO bad reflections found!\n" << G4endl;
        return;
    }
    G4cout << "❌ Total bad reflections: " << fTotalBadReflections << G4endl;
    for (auto& eventPair : fPhotonTraces) {
        int badInEvent = 0, totalInEvent = 0;
        for (auto& photonPair : eventPair.second) {
            for (double dot : photonPair.second.dotProducts) {
                totalInEvent++;
                if (dot > 0) badInEvent++;
            }
        }
        if (badInEvent > 0) G4cout << "  Event " << eventPair.first << ": " << badInEvent << "/" << totalInEvent << " bad" << G4endl;
    }
    G4cout << "=========================================================\n" << G4endl;
}

void G4d2oPhotonDiagnostics::PrintGlobalSummary() const {
    G4cout << "\n" << std::string(70, '=') << G4endl;
    G4cout << "GLOBAL PHOTON DIAGNOSTICS SUMMARY" << G4endl;
    G4cout << std::string(70, '=') << G4endl;
    G4cout << "  Total events:           " << fTotalEvents << G4endl;
    G4cout << "  Total photons:          " << fTotalPhotons << G4endl;
    G4cout << "  Total reflections:      " << fTotalReflections << G4endl;
    G4cout << "  Total PMT hits:         " << fTotalPMTHits << G4endl;
    G4cout << "  Total absorbed:         " << fTotalAbsorbed << G4endl;
    G4cout << "  Total bad reflections:  " << fTotalBadReflections << G4endl;
    G4cout << "  Events passing quality: " << fTotalQualityPass << G4endl;
    if (fTotalPhotons > 0) {
        G4cout << "\n  Average reflections per photon: " << (G4double)fTotalReflections / fTotalPhotons << G4endl;
        G4cout << "  PMT hit efficiency:            " << (G4double)fTotalPMTHits / fTotalPhotons * 100 << "%" << G4endl;
    }
    if (fTotalBadReflections > 0) {
        G4cout << "\n  ❌ WARNING: Bad reflections found!" << G4endl;
    } else {
        G4cout << "\n  ✅ GOOD: No bad reflections found." << G4endl;
    }
    G4cout << std::string(70, '=') << "\n" << G4endl;
}

// ============================================================
// FIXED: PrintChainSummary() - PHYSICALLY CORRECT CHAIN
// ============================================================

void G4d2oPhotonDiagnostics::PrintChainSummary() const {
    G4cout << "\n" << std::string(70, '=') << G4endl;
    G4cout << "PHOTON SURVIVAL CHAIN - GLOBAL SUMMARY" << G4endl;
    G4cout << std::string(70, '=') << G4endl;

    G4int totalC = fTotalCherenkovPhotons;

    // PHYSICS: Get termination counts (these tell us where photons actually ended)
    G4int totalReachedPMT = fTerminationCounts.at(kReachedPMT);
    G4int totalWaterAbsorbed = fTerminationCounts.at(kWaterAbsorbed);
    G4int totalTyvekAbsorbed = fTerminationCounts.at(kTyvekAbsorbed);
    G4int totalMaxReflections = fTerminationCounts.at(kMaxReflections);
    G4int totalOther = fTerminationCounts.at(kOther);

    // ============================================================
    // PHYSICS RELATIONSHIPS (must hold):
    // 1. Water survive = Cherenkov - WaterAbsorbed
    // 2. PMT reach ≤ Water survive (physical law!)
    // 3. PMT reach = termination count of kReachedPMT
    // ============================================================

    G4int waterSurvive = totalC - totalWaterAbsorbed;
    G4int pmtReach = totalReachedPMT;
    G4int totalReflections = fTotalTyvekReflections;
    G4int totalPE = fTotalPhotoelectrons;

    G4cout << "  Step 1: Cherenkov photons created:    " << totalC << G4endl;

    if (totalC > 0) {
        // Step 2: Water survival (fraction of Cherenkov that survived water)
        G4cout << "  Step 2: Survive water attenuation:    " << waterSurvive
               << " (" << (G4double)waterSurvive / totalC * 100 << "%)" << G4endl;

        // Step 3: Total Tyvek reflections (NOT a percentage - it's a count!)
        G4cout << "  Step 3: Total Tyvek reflections:      " << totalReflections
               << " (avg " << (G4double)totalReflections / totalC
               << " reflections per photon)" << G4endl;

        // Step 4: PMT reach (fraction of water survivors that reached PMTs)
        G4cout << "  Step 4: Reach PMTs:                    " << pmtReach
               << " (" << (G4double)pmtReach / waterSurvive * 100
               << "% of water survivors)" << G4endl;

        // Step 5: Photoelectrons (fraction of PMT hits that became PE)
        G4cout << "  Step 5: Become photoelectrons:         " << totalPE
               << " (" << (G4double)totalPE / pmtReach * 100
               << "% of PMT reach)" << G4endl;
    }

    // Averages per event
    if (fTotalEvents > 0) {
        G4cout << "\n  Averages per event:" << G4endl;
        G4cout << "    Cherenkov photons:  " << (G4double)totalC / fTotalEvents << G4endl;
        G4cout << "    Water survive:      " << (G4double)waterSurvive / fTotalEvents << G4endl;
        G4cout << "    Tyvek reflections:  " << (G4double)totalReflections / fTotalEvents << G4endl;
        G4cout << "    PMT reach:          " << (G4double)pmtReach / fTotalEvents << G4endl;
        G4cout << "    Photoelectrons:     " << (G4double)totalPE / fTotalEvents << G4endl;
    }

    // Sanity check - warn if the chain is physically impossible
    if (pmtReach > waterSurvive) {
        G4cout << "\n  ⚠️ WARNING: PMT reach (" << pmtReach << ") > Water survive ("
               << waterSurvive << ") - PHYSICALLY IMPOSSIBLE!" << G4endl;
        G4cout << "  Check that RecordPMTReach() is only called for photons that survive water."
               << G4endl;
    }

    G4cout << std::string(70, '=') << "\n" << G4endl;
}

void G4d2oPhotonDiagnostics::WriteDiagnostics() {
    std::ofstream file("photon_diagnostics.txt");
    if (!file.is_open()) return;

    file << "PHOTON DIAGNOSTICS OUTPUT" << std::endl;
    file << "=========================" << std::endl;
    file << "Total events: " << fTotalEvents << std::endl;
    file << "Total photons: " << fTotalPhotons << std::endl;
    file << "Total reflections: " << fTotalReflections << std::endl;
    file << "Total PMT hits: " << fTotalPMTHits << std::endl;
    file << "Total absorbed: " << fTotalAbsorbed << std::endl;
    file << "Total bad reflections: " << fTotalBadReflections << std::endl;

    file << "\nPHOTON TERMINATION REASONS:" << std::endl;
    const char* names[] = {"Reached PMT", "Water Absorbed", "Tyvek Absorbed", "Max Reflections", "Boundary Kill", "Other"};
    for (int i = 0; i <= kOther; ++i) {
        G4int count = fTerminationCounts.at((TerminationReason)i);
        file << "  " << names[i] << ": " << count << std::endl;
    }

    file << "\nPHOTON SURVIVAL CHAIN" << std::endl;
    file << "Cherenkov photons: " << fTotalCherenkovPhotons << std::endl;
    file << "Survive water: " << (fTotalCherenkovPhotons - fTerminationCounts.at(kWaterAbsorbed)) << std::endl;
    file << "Tyvek reflections: " << fTotalTyvekReflections << std::endl;
    file << "Reach PMTs: " << fTerminationCounts.at(kReachedPMT) << std::endl;
    file << "Photoelectrons: " << fTotalPhotoelectrons << std::endl;

    file.close();
    G4cout << "Diagnostics written to photon_diagnostics.txt" << G4endl;
}
