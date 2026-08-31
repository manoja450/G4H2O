#include "G4d2oElectronGunvD.hh"
#include "inputVariables.hh"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH2D.h"
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

G4d2oElectronGunvD::G4d2oElectronGunvD()
{
    G4cout << "\tConstructing G4d2oElectronGunvD..." ;
    
    // Set random seed variables
    input = inputVariables::GetIVPointer();
    G4int irand = input->GetRandomStatus();
    if(irand==0) gRandom->SetSeed(1);
    if(irand==1) gRandom->SetSeed(0);
    totalEvents = input->GetNumberOfEvents();

    // Initialize a few things
    G4int n_particle = 1;
    particleGun = new G4ParticleGun(n_particle);
    
    // Set the particle name
    G4String particleName = "e-"; // return to THIS
    // G4String particleName = "mu-";

    // set up the gun
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    particleGun->SetParticleDefinition(particleTable->FindParticle(particleName));
    
    SourcePosition.setX(0);
    SourcePosition.setY(0);
    SourcePosition.setZ(0);
    
    G4double thetaAngle1 = 0.0;   // in degrees
    G4double thetaAngle2 = 180.0; // in degrees
    
    cthetaRange1 = cos(thetaAngle1*deg);
    cthetaRange2 = cos(thetaAngle2*deg);
    
    G4d2oNeutrinoAlley *detCon = (G4d2oNeutrinoAlley*)
        G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    G4d2oDetector *theDet = (G4d2oDetector*)detCon->GetDetectorPtr();
    // Inner Tank of D20
    // tankSize.set(theDet->TankX(),theDet->TankY(),theDet->TankZ());
    // Outer Tank of H20
    tankSize.set(theDet->OuterTankX(),theDet->OuterTankY(),theDet->OuterTankZ());

    etfTimer = new ReplayTools();
    etfTimer->PrepareETFTimer(5000, input->GetNumberOfEvents()); // time in ms

    // -------- NEW: load fluxWeight.root and fluxW histogram --------
    angleOffset = -std::atan(8.077 / 16.76);  // same offset as in readFW.C

    // Adjust path if needed
    const char* fluxPath = "./xscnData/fluxWeight.root";

    fluxFile = TFile::Open(fluxPath, "READ");
    if (!fluxFile || fluxFile->IsZombie()) {
        G4Exception("G4d2oElectronGunvD::G4d2oElectronGunvD",
                    "FluxFileError", FatalException,
                    ("Cannot open flux weight file: " + G4String(fluxPath)).c_str());
    }

    fluxW = dynamic_cast<TH2D*>(fluxFile->Get("fluxW"));
    if (!fluxW) {
        G4Exception("G4d2oElectronGunvD::G4d2oElectronGunvD",
                    "FluxHistError", FatalException,
                    "Histogram 'fluxW' not found in flux weight file.");
    }

    // Optional: normalize like in analysis macro
    double integral = fluxW->Integral();
    if (integral > 0.0) {
        fluxW->Scale(1.0 / integral);
    } else {
        G4Exception("G4d2oElectronGunvD::G4d2oElectronGunvD",
                    "FluxHistError", FatalException,
                    "Histogram 'fluxW' has non-positive integral.");
    }

    G4cout << "\n[ElectronGun] Loaded fluxW from " << fluxPath
           << " (integral = " << fluxW->Integral() << ")" << G4endl;
    // ---------------------------------------------------------------

    G4cout << "done." << G4endl;

}// END of constructor


G4d2oElectronGunvD::~G4d2oElectronGunvD()
{
    
	G4cout<<"Deleting G4d2oElectronGunvD...";
	
//    delete particleGun;
	
	G4cout<<"done."<<G4endl;
	
}//END of destructor

void G4d2oElectronGunvD::GeneratePrimaries(G4Event* anEvent)
{
    theEventNum = anEvent->GetEventID();
    etfTimer->SetCurrentEvent(theEventNum);
    gSystem->ProcessEvents();
    
    if(anEvent->GetEventID()==0 && input->GetPrintStatus()!=3)
        etfTimer->StartUpdateTimer();
        
    G4double px, py, pz, v0x, v0y, v0z, v1x, v1y, v1z, vsize, w1x, w1y, w1z, wsize;
    v0x = (0.1182589088);
    v0y = (-0.37975131476);
    v0z = (-0.91749864818);
    G4double ctheta, stheta, phi;

    // -------- NEW: sample (E, new_cos_theta) from fluxW --------
    if (!fluxW) {
        G4Exception("G4d2oElectronGunvD::GeneratePrimaries",
                    "NoFluxHist", FatalException,
                    "Flux weight histogram not loaded.");
    }

    double E_sample = 0.0;      // in MeV (numeric, as in the histogram)
    double new_cos_theta = 1.0; // dimensionless

    // x-axis = energy, y-axis = new_cos_theta (as in readFW.C)
    fluxW->GetRandom2(E_sample, new_cos_theta);
//     G4cout << "\tnew_cos_theta = " << new_cos_theta;
    // Convert new_cos_theta back to angle (this is new_angle in macro)
    if(new_cos_theta < -1.0){new_cos_theta = -1.0;}
    if(new_cos_theta > +1.0){new_cos_theta = +1.0;}
    double theta_prime = std::acos(new_cos_theta);
//     G4cout << "\tthetaprime = " << theta_prime;
    while(true){

    

    while(true){
        v1x = G4UniformRand();
        v1y = G4UniformRand();
        v1z = G4UniformRand();
        vsize = std::sqrt(v1x * v1x + v1y * v1y + v1z * v1z);

        v1x = v1x / vsize;
        v1y = v1y / vsize;
        v1z = v1z / vsize;

 //             G4cout << "\tvsquare = " << v1x*v1x+v1y*v1y+v1z*v1z;

        if (v0x != v1x || v0y != v1y || v0z != v1z){
            break;        
    }


    }

    w1x = v1y * v0z - v1z * v0y;
    w1y = v1z * v0x - v1x * v0z;
    w1z = v1x * v0y - v1y * v0x;

    wsize = std::sqrt(w1x * w1x + w1y * w1y + w1z * w1z);

    w1x = w1x / wsize;
    w1y = w1y / wsize;
    w1z = w1z / wsize;

   //      G4cout << "\twsquare = " << w1x*w1x+w1y*w1y+w1z*w1z;

    px = v0x * std::cos(theta_prime) + w1x * std::sin(theta_prime);
    py = v0y * std::cos(theta_prime) + w1y * std::sin(theta_prime);
    pz = v0z * std::cos(theta_prime) + w1z * std::sin(theta_prime);

 //    G4cout << "\tpx = " << px;
 //    G4cout << "\tpy = " << py;
   //       G4cout << "\tpz = " << pz;
    // In readFW.C: new_angle = theta + offset => theta = new_angle - offset
//    double theta = theta_prime - angleOffset;
//    double theta = std::asin(1.0) - theta_prime;
    
    // Optional: clamp to [0, pi]
//    if (theta < 0.0) theta = 0.0;
//    if (theta > CLHEP::pi) theta = CLHEP::pi;

    // Random azimuth
  //  phi = 2.0 * TMath::Pi() * G4UniformRand();
  //  auto theta = 1.0 * TMath::Pi() * G4UniformRand();// comment this out

    
  //  stheta = std::sin(theta);
  //  ctheta = std::cos(theta);

  //  pz = ctheta;
  //  py = stheta * std::cos(phi);
  //  px = stheta * std::sin(phi);
   //   G4cout << "\tpsquare = " << px*px+py*py+pz*pz;
  if (px*px+py*py+pz*pz > .999){
    break;
  }

    }
    initDir.set(px, py, pz);

    // Energy in Geant4 units (assuming histogram x-axis is MeV)
    sourceEnergy = E_sample * MeV;
    particleGun->SetParticleEnergy(sourceEnergy);
    // -----------------------------------------------------------

    G4Navigator* Navigator =
        G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

    int ntries = 0;

    G4double randvarn1; // rho, radius
    G4double randvarn2; // theta, xy angle
    G4double randvarn3; // length of the tank
    while(true){
        randvarn1 = 44.45*std::sqrt(G4UniformRand());                // rho, radius
        randvarn2 = 2.0*TMath::Pi()*G4UniformRand();                 // theta, xy angle
        randvarn3 = 89.9275*(-1.0+2.0*G4UniformRand());              // length of the tank
        SourcePosition.set( (0.0+randvarn1*std::cos(randvarn2))*cm,
                            (110.73384+randvarn1*std::sin(randvarn2))*cm,
                            (-80.74526+randvarn3)*cm );
        G4VPhysicalVolume* volume = Navigator->LocateGlobalPointAndSetup(SourcePosition);
        if(volume->GetName() == "d2oPhysV"
//           || volume->GetName() == "h2oPhysV"
//           || volume->GetName() == "acrylicPhysV"
          ){
            break;
        }
        ntries++;
    }

    particleGun->SetParticlePosition(SourcePosition);
    particleGun->SetParticleMomentumDirection(initDir);
    
    particleGun->GeneratePrimaryVertex(anEvent);
    
}// END of GeneratePrimaries()


