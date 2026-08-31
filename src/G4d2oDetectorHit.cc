#include "G4d2oDetectorHit.hh"

#include "G4UnitsTable.hh"

G4Allocator<G4d2oDetectorHit> G4d2oDetectorHitAllocator;

G4d2oDetectorHit::G4d2oDetectorHit()
{
    hitTrackID            = 0;
    hitPixelNumber        = 0;
    hitDetectorNumber     = 0;
    hitTrackParentID      = 0;
    hitTotalEnergyDeposit = 0;
    hitGlobalTime         = 0;
    hitKineticEnergy      = 0;
    hitPDGCharge          = 0;
    hitLeptonNumber       = 0;
    hitBaryonNumber       = 0;
    hitDeltaTime          = 0;
    snprintf(hitParticleName, sizeof(hitParticleName), "%s", "");
    snprintf(hitProcessName, sizeof(hitProcessName), "%s", "");
    snprintf(hitProcessTypeName, sizeof(hitProcessTypeName), "%s", "");
    snprintf(hitMaterialName, sizeof(hitMaterialName), "%s", "");
}//END of constructor

G4d2oDetectorHit::~G4d2oDetectorHit()
{
    //  G4cout << "Instance of G4d2oDetectorHit Destructed!" << G4endl;
    
}//END of destructor

void G4d2oDetectorHit::SetPixelNumber(G4int thePixelNumber){ hitPixelNumber = thePixelNumber; }
G4int G4d2oDetectorHit::GetPixelNumber(void) { return hitPixelNumber; }
//END of Set/Get PixelNumber()

void G4d2oDetectorHit::SetDetectorNumber(G4int theDetectorNumber){ hitDetectorNumber = theDetectorNumber; }
G4int G4d2oDetectorHit::GetDetectorNumber(void) { return hitDetectorNumber; }
//END of Set/Get DetectorNumber()

void G4d2oDetectorHit::SetTrackID(G4int theTrackID){ hitTrackID = theTrackID; }
G4int G4d2oDetectorHit::GetTrackID(void) { return hitTrackID; }
//END of Set/Get TrackID()

void G4d2oDetectorHit::SetTrackParentID(G4int theTrackParentID){ hitTrackParentID = theTrackParentID; }
G4int G4d2oDetectorHit::GetTrackParentID(void) { return hitTrackParentID; }
//END of Set/Get TrackParentID()

void G4d2oDetectorHit::SetScatA(G4int theScatA){ hitScatA = theScatA; }
G4int G4d2oDetectorHit::GetScatA(void) { return hitScatA; }
//END of Set/Get ScatA()

void G4d2oDetectorHit::SetScatZ(G4int theScatZ){ hitScatZ = theScatZ; }
G4int G4d2oDetectorHit::GetScatZ(void) { return hitScatZ; }
//END of Set/Get ScatZ()

void G4d2oDetectorHit::SetTotalEnergyDeposit(G4double theEnergyDeposit) {  hitTotalEnergyDeposit = theEnergyDeposit; }
G4double G4d2oDetectorHit::GetTotalEnergyDeposit(void) { return hitTotalEnergyDeposit; }
//END of Set/Get TotalEnergyDeposit()

void G4d2oDetectorHit::SetGlobalTime(G4double theGlobalTime) { hitGlobalTime = theGlobalTime; }
G4double G4d2oDetectorHit::GetGlobalTime(void) { return hitGlobalTime; }
//END of Set/Get GlobalTime()

void G4d2oDetectorHit::SetKineticEnergy(G4double theKineticEnergy) { hitKineticEnergy = theKineticEnergy; }
G4double G4d2oDetectorHit::GetKineticEnergy(void) { return hitKineticEnergy; }
//END of Set/Get KineticEnergy()

void G4d2oDetectorHit::SetPDGCharge(G4double thePDGCharge) { hitPDGCharge = thePDGCharge; }
G4double G4d2oDetectorHit::GetPDGCharge(void) { return hitPDGCharge; }
//END of Set/Get PDGCharge()

void G4d2oDetectorHit::SetBaryonNumber(G4double theBaryonNumber) { hitBaryonNumber = theBaryonNumber; }
G4double G4d2oDetectorHit::GetBaryonNumber(void) { return hitBaryonNumber; }
//END of Set/Get BaryonNumber()

void G4d2oDetectorHit::SetLeptonNumber(G4double theLeptonNumber) { hitLeptonNumber = theLeptonNumber; }
G4double G4d2oDetectorHit::GetLeptonNumber(void) { return hitLeptonNumber; }
//END of Set/Get LeptonNumber()

void G4d2oDetectorHit::SetDeltaTime(G4double theDeltaTime) { hitDeltaTime = theDeltaTime; }
G4double G4d2oDetectorHit::GetDeltaTime(void) { return hitDeltaTime; }
//END of Set/Get GlobalTime()

void G4d2oDetectorHit::SetCalcA(G4double theCalcA) {  hitCalcA = theCalcA; }
G4double G4d2oDetectorHit::GetCalcA(void) { return hitCalcA; }
//END of Set/Get CalcA()

void G4d2oDetectorHit::SetParticleName(char *theParticleName){ snprintf(hitParticleName, sizeof(hitParticleName), "%s", theParticleName); }
char* G4d2oDetectorHit::GetParticleName(void) { return hitParticleName; }
//END of Set/Get ParticleName()

void G4d2oDetectorHit::SetProcessName(char *theProcessName){ snprintf(hitProcessName, sizeof(hitProcessName), "%s", theProcessName); }
char* G4d2oDetectorHit::GetProcessName(void) { return hitProcessName; }
//END of Set/Get ProcessName()

void G4d2oDetectorHit::SetProcessTypeName(char *theProcessTypeName){ snprintf(hitProcessTypeName, sizeof(hitProcessTypeName), "%s", theProcessTypeName); }
char* G4d2oDetectorHit::GetProcessTypeName(void) { return hitProcessTypeName; }
//END of Set/Get ProcessTypeName()

void G4d2oDetectorHit::SetMaterialName(char *theMaterialName){ snprintf(hitMaterialName, sizeof(hitMaterialName), "%s", theMaterialName); }
char* G4d2oDetectorHit::GetMaterialName(void) { return hitMaterialName; }
//END of Set/Get MaterialName()

void G4d2oDetectorHit::SetCreatorProcessName(char *theCreatorProcessName){ snprintf(hitCreatorProcessName, sizeof(hitCreatorProcessName), "%s", theCreatorProcessName); }
char* G4d2oDetectorHit::GetCreatorProcessName(void) { return hitCreatorProcessName; }
//END of Set/Get CreatorProcessName()

void G4d2oDetectorHit::SetCreatorProcessType(char *theCreatorProcessType){ snprintf(hitCreatorProcessType, sizeof(hitCreatorProcessType), "%s", theCreatorProcessType); }
char* G4d2oDetectorHit::GetCreatorProcessType(void) { return hitCreatorProcessType; }
//END of Set/Get CreatorProcessType()

void G4d2oDetectorHit::SetMomentumDirection(G4ThreeVector theMomentumDirection){ hitMomentumDirection = theMomentumDirection; }
G4ThreeVector G4d2oDetectorHit::GetMomentumDirection(void) { return hitMomentumDirection; }
//END of Set/Get MomentumDirection()

void G4d2oDetectorHit::SetPostMomentumDirection(G4ThreeVector thePostMomentumDirection){ hitPostMomentumDirection = thePostMomentumDirection; }
G4ThreeVector G4d2oDetectorHit::GetPostMomentumDirection(void) { return hitPostMomentumDirection; }
//END of Set/Get PostMomentumDirection()

void G4d2oDetectorHit::SetPosition(G4ThreeVector thePosition){ hitPosition = thePosition; }
G4ThreeVector G4d2oDetectorHit::GetPosition(void) { return hitPosition; }
//END of Set/Get Position()

void G4d2oDetectorHit::SetPrePosition(G4ThreeVector thePrePosition){ hitPrePosition = thePrePosition; }
G4ThreeVector G4d2oDetectorHit::GetPrePosition(void) { return hitPrePosition; }
//END of Set/Get Position()

void G4d2oDetectorHit::SetVertex(G4ThreeVector theVertex){ hitVertex = theVertex; }
G4ThreeVector G4d2oDetectorHit::GetVertex(void) { return hitVertex; }
//END of Set/Get Vertex()

void G4d2oDetectorHit::Draw(void)
{
}//END of Draw()

void G4d2oDetectorHit::Print(void)
{
}//END of Print()
