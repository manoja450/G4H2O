#include "G4d2oMuonGun.hh"
#include "inputVariables.hh"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TSystem.h"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"

G4d2oMuonGun::G4d2oMuonGun()
{
    
    G4cout << "\tConstructing G4d2oMuonGun..." ;
    
    //Set random seed variables
    input = inputVariables::GetIVPointer();
    G4int irand = input->GetRandomStatus();
    if(irand==0) gRandom->SetSeed(1);
    if(irand==1) gRandom->SetSeed(0);
    totalEvents = input->GetNumberOfEvents();

    //Initialize a few things
    G4int n_particle = 1;
    particleGun = new G4ParticleGun(n_particle);
//    sourceEnergy = 4000.0*MeV;
    sourceEnergy = 10.0*MeV;
//    sourceEnergy = 140.0*keV;
    
    //Set the particle name
    G4String particleName = "mu-"; // return to THIS
    //G4String particleName = "mu-";

    //set up the gun
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    particleGun->SetParticleDefinition(particleTable->FindParticle(particleName));
    
    SourcePosition.setX(0);
    SourcePosition.setY(0);
    SourcePosition.setZ(0);
    
    G4double thetaAngle1 = 0.0; //in degrees
    G4double thetaAngle2 = 180.0; //in degrees
    
    cthetaRange1 = cos(thetaAngle1*deg);
    cthetaRange2 = cos(thetaAngle2*deg);
    
    G4d2oNeutrinoAlley *detCon = (G4d2oNeutrinoAlley*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    G4d2oDetector *theDet = (G4d2oDetector*)detCon->GetDetectorPtr();
    // Inner Tank of D20
    // tankSize.set(theDet->TankX(),theDet->TankY(),theDet->TankZ());
    // Outer Tank of H20
    tankSize.set(theDet->OuterTankX(),theDet->OuterTankY(),theDet->OuterTankZ());

    etfTimer = new ReplayTools();
    etfTimer->PrepareETFTimer(5000, input->GetNumberOfEvents()); //time in ms
    
    G4cout << "done." << G4endl;

}//END of constructor

G4d2oMuonGun::~G4d2oMuonGun()
{
    
	G4cout<<"Deleting G4d2oMuonGun...";
	
//    delete particleGun;
	
	G4cout<<"done."<<G4endl;
	
}//END of destructor

void G4d2oMuonGun::GeneratePrimaries(G4Event* anEvent)
{
    
    theEventNum = anEvent->GetEventID();
    etfTimer->SetCurrentEvent(theEventNum);
    gSystem->ProcessEvents();
    
    if(anEvent->GetEventID()==0 && input->GetPrintStatus()!=3) etfTimer->StartUpdateTimer();
        
    G4double px, py, pz;
    G4double ctheta, stheta, phi;

    sourceEnergy = (1000.0)*MeV;

    particleGun->SetParticleEnergy( sourceEnergy ); // return to THIS
  
    G4Navigator* Navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

    int ntries = 0;

    G4double randvarn1; // rho, radius
    G4double randvarn2; // theta, xy angle
    G4double randvarn3; // length of the tank
    while(true){//ntries < anEvent->GetEventID()){
      randvarn1 = 44.45*sqrt(G4UniformRand()); // rho, radius
      randvarn2 = 2.0*TMath::Pi()*G4UniformRand(); // theta, xy angle
      randvarn3 = 89.9275*(-1.0+2.0*G4UniformRand()); // length of the tank
      //SourcePosition.set( (0.0+randvarn1*cos(randvarn2))*cm, (110.73384+randvarn1*sin(randvarn2))*cm, (-80.74526+randvarn3)*cm  );
      SourcePosition.set( (0.0+0.0*randvarn1*cos(randvarn2))*cm, (110.73384+0.0*randvarn1*sin(randvarn2))*cm, (-80.74526+89.9275+100)*cm  );
      G4VPhysicalVolume* volume = Navigator->LocateGlobalPointAndSetup(SourcePosition);

       if(volume->GetName() == "h2oPhysV"
          || volume->GetName() == "d2oPhysV"
          || volume->GetName() == "acrylicPhysV"
          || 1==1
          ){
//      if(volume->GetName() == "d2oPhysV"){
        break;}
      ntries++;

    //  std::cout << "Sampling gun "  << ntries << std::endl;

    }
  
//    particleGun->SetParticlePosition( {(0.0+randvarn1*cos(randvarn2))*cm, (110.73384+randvarn1*sin(randvarn2))*cm, (-80.74526+randvarn3)*cm} );
    particleGun->SetParticlePosition( SourcePosition );
 //   particleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, -1.0));
    ctheta = G4UniformRand()*(cthetaRange1-cthetaRange2) + cthetaRange2;
    phi = 2.0*TMath::Pi()*G4UniformRand();
    
    stheta = sqrt( 1 - pow(ctheta,2) );
    
    pz = ctheta;
    py = stheta*cos(phi);
    px = stheta*sin(phi);
    
    
    //initDir.set(px,py,pz);
    G4double gamble = G4UniformRand();

    if (gamble >= 0.6)
    {initDir.set(-(1/sqrt(8)),-(0/sqrt(8)),-(sqrt(7.0/8.0)));}
    if (gamble < 0.6 && gamble >= 0.3)
    {initDir.set(-(0/sqrt(8)),+(1/sqrt(8)),-(sqrt(7.0/8.0)));}
    if (gamble < 0.3 && gamble >= 0.1)
    {initDir.set(+(0/sqrt(8)),-(1/sqrt(8)),-(sqrt(7.0/8.0)));}
    if (gamble < 0.1)
    {initDir.set(+(1/sqrt(8)),+(0/sqrt(8)),-(sqrt(7.0/8.0)));}
    
    
    particleGun->SetParticleMomentumDirection(initDir);
    
    particleGun->GeneratePrimaryVertex(anEvent);
    
}//END of GeneratePrimaries()

