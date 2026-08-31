
#include "G4d2oMichelElectronGun.hh"
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

G4d2oMichelElectronGun::G4d2oMichelElectronGun()
{

  G4cout << "\tConstructing G4d2oMichelElectronGun..." ;

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
  ////sourceEnergy = 10.0*MeV;
  //    sourceEnergy = 140.0*keV;
  sourceEnergy = (0.0+G4UniformRand()*53.0)*MeV;

  //Set the particle name
  G4String particleName = "e-";
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
  //tankSize.set(theDet->TankX(),theDet->TankY(),theDet->TankZ());
  // Outer Tank of H20
  tankSize.set(theDet->OuterTankX(),theDet->OuterTankY(),theDet->OuterTankZ());

  etfTimer = new ReplayTools();
  etfTimer->PrepareETFTimer(5000, input->GetNumberOfEvents()); //time in ms

  // Some parameters needed for electron energy sampling
  EMMU = 105.66 * MeV; // muon mass
  EMASS = 0.511 * MeV;
  michel_rho   = 0.75;
  michel_delta = 0.75;
  michel_xsi   = 1.00;
  michel_eta   = 0.00;

  G4cout << "done." << G4endl;

}//END of constructor

G4d2oMichelElectronGun::~G4d2oMichelElectronGun()
{

  G4cout<<"Deleting G4d2oMichelElectronGun...";

  //    delete particleGun;

  G4cout<<"done."<<G4endl;

}//END of destructor

// Taken from G4MuonDecayChannelWithSpin::DecayIt
// Assumes V-A coupling with 1st order radiative corrections,
//  the SM Michel parameters
void G4d2oMichelElectronGun::GeneratePrimaries(G4Event* anEvent)
{

  theEventNum = anEvent->GetEventID();
  etfTimer->SetCurrentEvent(theEventNum);
  gSystem->ProcessEvents();

  if(anEvent->GetEventID()==0 && input->GetPrintStatus()!=3) etfTimer->StartUpdateTimer();

  G4double px, py, pz;
  G4double ctheta, stheta, phi;

  // FIrst let's generate the energy
  G4double rndm, x;
  G4double FG;
  G4double FG_max = 2.00;

  G4double W_mue = (EMMU*EMMU + EMASS*EMASS) / (2.0*EMMU);
  G4double x0    = EMASS/W_mue;
  G4double x0_squared = x0 * x0;

  // ***************************************************
  //     x0 <= x <= 1.   and   -1 <= y <= 1
  //
  //     F(x,y) = f(x)*g(x,y);   g(x,y) = 1.+g(x)*y
  // ***************************************************

  // ***** sampling F(x,y) directly (brute force) *****
  do {
    // Sample the positron energy by sampling from F
    rndm = G4UniformRand();
    x = x0 + rndm*(1.0 - x0);
    G4double x_squared = x*x;

    G4double F_IS = 1./6.*(-2.*x_squared+3.*x-x0_squared);

    G4double F_AS = 1./6.*std::sqrt(x_squared-x0_squared)*(2.*x-2.+std::sqrt(1.-x0_squared));

    G4double G_IS = 2./9.*(michel_rho-0.75)*(4.*x_squared-3.*x-x0_squared);
    G_IS = G_IS + michel_eta*(1.-x)*x0;

    G4double G_AS = 3.*(michel_xsi-1.)*(1.-x);
    G_AS = G_AS+2.*(michel_xsi*michel_delta-0.75)*(4.*x-4.+std::sqrt(1.-x0_squared));
    G_AS = 1./9.*std::sqrt(x_squared-x0_squared)*G_AS;

    F_IS = F_IS + G_IS;
    F_AS = F_AS + G_AS;

    // *** Radiative Corrections ***
    G4double R_IS = F_c(x, x0);
    G4double F = 6.0*F_IS + R_IS/std::sqrt(x_squared - x0_squared);

    // *** Radiative Corrections cont. ***
    G4double R_AS = F_theta(x, x0);

    rndm = G4UniformRand();
    ctheta = G4UniformRand()*(cthetaRange1-cthetaRange2) + cthetaRange2;

    G4double G = 6.0 * F_AS - R_AS/std::sqrt(x_squared - x0_squared);
    FG = std::sqrt(x_squared - x0_squared) * F * (1.0 + (G/F) * ctheta);

    if (FG > FG_max) {
      G4cout<<"***Problem in generating Michel spectrum***: FG > FG_max\n";
      FG_max = FG;
    }

    rndm = G4UniformRand();
  } while (FG < rndm * FG_max);

  sourceEnergy = x * W_mue;
  if (sourceEnergy < EMASS) sourceEnergy = EMASS;
  particleGun->SetParticleEnergy( sourceEnergy );

  G4Navigator* Navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

  int ntries = 0;
  while(false){
    SourcePosition.set(tankSize.x()*(G4UniformRand()-0.5),
        tankSize.y()*(G4UniformRand()-0.5),
        tankSize.z()*(G4UniformRand()-0.5));
    G4VPhysicalVolume* volume = Navigator->LocateGlobalPointAndSetup(SourcePosition);
    if(volume->GetName() == "h2oPhysV"
        || volume->GetName() == "d2oPhysV"
        || volume->GetName() == "acrylicPhysV"
      )
      //        if(volume->GetName() == "d2oPhysV")
      break;
    ntries++;
  }

  //  particleGun->SetParticlePosition( SourcePosition );
  //  particleGun->SetParticlePosition( {0*m, 1*m, -1*m}  );
  G4double randvarn1 = 44.45*sqrt(G4UniformRand()); // rho, radius
  G4double randvarn2 = 2.0*TMath::Pi()*G4UniformRand(); // theta, xy angle
  G4double randvarn3 = 89.9275*(-1.0+2.0*G4UniformRand()); // length of the tank

  particleGun->SetParticlePosition( {(0.0+randvarn1*cos(randvarn2))*cm, (110.73384+randvarn1*sin(randvarn2))*cm, (-80.74526+randvarn3)*cm}  );
  //particleGun->SetParticlePosition( {(0.0+randvarn1*cos(randvarn2))*cm, (110.73384+randvarn1*sin(randvarn2))*cm, (-80.74526-(3.0/4.0)*(89.9275)+randvarn3/4.00)*cm}  );

  ctheta = G4UniformRand()*(cthetaRange1-cthetaRange2) + cthetaRange2;
  phi = 2.0*TMath::Pi()*G4UniformRand();

  stheta = sqrt( 1 - pow(ctheta,2) );

  pz = ctheta;
  py = stheta*cos(phi);
  px = stheta*sin(phi);


  ////    ///// temp for fixed initial direction /////
  //    px = 0.0;
  //    py = -1.0;
  //    pz = 0.0;
  //particleGun->SetParticlePosition(G4ThreeVector(-5.0*cm,0.0,92.5*cm+3.8*mm));
  //particleGun->SetParticlePosition(G4ThreeVector(0.0*cm,200.0*cm,0.0*cm));
  ////    ///// end temp /////

  initDir.set(px,py,pz);

  particleGun->SetParticleMomentumDirection(initDir);

  particleGun->GeneratePrimaryVertex(anEvent);

}//END of GeneratePrimaries()

G4double G4d2oMichelElectronGun::F_c(G4double x, G4double x0)
{
  G4double f_c = 0.0;

  G4double omega = std::log(EMMU / EMASS);

  f_c = (5.+17.*x-34.*x*x)*(omega+std::log(x))-22.*x+34.*x*x;
  f_c = (1.-x)/(3.*x*x)*f_c;
  f_c = (6.-4.*x)*R_c(x)+(6.-6.*x)*std::log(x) + f_c;
  f_c = (CLHEP::fine_structure_const/CLHEP::twopi) * (x*x-x0*x0) * f_c;

  return f_c;
}

G4double G4d2oMichelElectronGun::F_theta(G4double x, G4double x0)
{
  G4double f_theta = 0.0;
  G4double omega = std::log(EMMU / EMASS);

  f_theta = (1.+x+34*x*x)*(omega+std::log(x))+3.-7.*x-32.*x*x;
  f_theta = f_theta + ((4.*(1.-x)*(1.-x))/x)*std::log(1.-x);
  f_theta = (1.-x)/(3.*x*x) * f_theta;
  f_theta = (2.-4.*x)*R_c(x)+(2.-6.*x)*std::log(x)-f_theta;
  f_theta = (CLHEP::fine_structure_const/CLHEP::twopi) * (x*x-x0*x0) * f_theta;

  return f_theta;
}

G4double G4d2oMichelElectronGun::R_c(G4double x) {
  G4int n_max = (int)(100.*x);

  if(n_max<10)n_max=10;

  G4double L2 = 0.0;

  for(G4int n=1; n<=n_max; n++){
    L2 += std::pow(x,n)/(n*n);
  }

  G4double omega = std::log(EMMU/EMASS);

  G4double r_c;

  r_c = 2.*L2-(pi*pi/3.)-2.;
  r_c = r_c + omega * (1.5+2.*std::log((1.-x)/x));
  r_c = r_c - std::log(x)*(2.*std::log(x)-1.);
  r_c = r_c + (3.*std::log(x)-1.-1./x)*std::log(1.-x);

  return r_c;
}
