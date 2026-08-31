
#include "G4d2oPhotonGun.hh"
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
#include "globals.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"

G4d2oPhotonGun::G4d2oPhotonGun()
{
    
    G4cout << "\tConstructing G4d2oPhotonGun..." ;
    
    //Set random seed variables
    input = inputVariables::GetIVPointer();
    G4int irand = input->GetRandomStatus();
    if(irand==0) gRandom->SetSeed(1);
    if(irand==1) gRandom->SetSeed(0);
    totalEvents = input->GetNumberOfEvents();

    //Initialize a few things
    G4int n_particle = 1;
    particleGun = new G4ParticleGun(n_particle);
    sourceEnergy = 3.0*eV;
    
    //Set the particle name
    G4String particleName = "opticalphoton";

    //set up the gun
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    particleGun->SetParticleDefinition(particleTable->FindParticle(particleName));
    
    SourcePosition.setX(0);
    SourcePosition.setY(0);
    SourcePosition.setZ(0);
    
    //uncomment next 2 lines for isotropic
    G4double thetaAngle1 = 0.0; //in degrees
    G4double thetaAngle2 = 180.0; //in degrees
//    //uncomment next 2 lines to specify a different polar angle range
//    G4double thetaAngle1 = 30.0; //in degrees
//    G4double thetaAngle2 = 35.0; //in degrees
    
//    photonsPerEvent = 10000;
    photonsPerEvent = 1;

    cthetaRange1 = cos(thetaAngle1*deg);
    cthetaRange2 = cos(thetaAngle2*deg);
    
    etfTimer = new ReplayTools();
    etfTimer->PrepareETFTimer(5000, input->GetNumberOfEvents()); //time in ms
    
    G4cout << "done." << G4endl;

}//END of constructor

G4d2oPhotonGun::~G4d2oPhotonGun()
{
    
	G4cout<<"Deleting G4d2oPhotonGun...";
	
//    delete particleGun;
	
	G4cout<<"done."<<G4endl;
	
}//END of destructor

void G4d2oPhotonGun::GeneratePrimaries(G4Event* anEvent)
{
    
    theEventNum = anEvent->GetEventID();
    etfTimer->SetCurrentEvent(theEventNum);
    gSystem->ProcessEvents();
    
    if(anEvent->GetEventID()==0 && input->GetPrintStatus()!=3) etfTimer->StartUpdateTimer();
        
    G4double px, py, pz;
    G4double ctheta, stheta, phi;
    
    particleGun->SetParticleEnergy( sourceEnergy );
    //    particleGun->SetParticlePosition( SourcePosition );
    //    particleGun->SetParticlePosition({0*m,1*m,-1*m});
    particleGun->SetParticlePosition({0*m,0*1*m,0*(-1)*m});

    for(G4int iPhoton=0; iPhoton<photonsPerEvent; iPhoton++){
        ctheta = G4UniformRand()*(cthetaRange1-cthetaRange2) + cthetaRange2;
        phi = 2.0*TMath::Pi()*G4UniformRand();
        
        stheta = sqrt( 1 - pow(ctheta,2) );
        
        pz = ctheta;
        py = stheta*cos(phi);
        px = stheta*sin(phi);
        
        initDir.set(px,py,pz);
    
        particleGun->SetParticleMomentumDirection(initDir);
    
        //optical photons need a polarization. randomizing it here
        if(particleGun->GetParticleDefinition()->GetParticleName()=="opticalphoton"){
            G4ThreeVector polVector = initDir.orthogonal();
            polVector.rotate( G4UniformRand()*TMath::TwoPi(), initDir );
            particleGun->SetParticlePolarization(polVector);
        }
    
        particleGun->GeneratePrimaryVertex(anEvent);
    }
    
}//END of GeneratePrimaries()

