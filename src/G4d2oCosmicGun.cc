
#include "G4d2oCosmicGun.hh"
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

G4d2oCosmicGun::G4d2oCosmicGun()
{
    
    G4cout << "\tConstructing G4d2oCosmicGun..." ;
    
    //Set random seed variables
    input = inputVariables::GetIVPointer();
    G4int irand = input->GetRandomStatus();
    if(irand==0) gRandom->SetSeed(1);
    if(irand==1) gRandom->SetSeed(0);
    totalEvents = input->GetNumberOfEvents();

    //Initialize a few things
    G4int n_particle = 1;
    particleGun = new G4ParticleGun(n_particle);
    
    // Read the cry input file
    std::ifstream inputFile;
    std::string inputfile("crysetup.txt");
    inputFile.open(inputfile,std::ios::in);
    char buffer[1000];
    
    if (inputFile.fail()) {
        if( inputfile !="")  //....only complain if a filename was given
            G4cout << "PrimaryGeneratorAction: Failed to open CRY input file= " << inputfile << G4endl;
        InputState=-1;
    }else{
#ifdef USE_CRY
        std::string setupString("");
        while ( !inputFile.getline(buffer,1000).eof()) {
            setupString.append(buffer);
            setupString.append(" ");
        }
        G4String crydata = gSystem->Getenv("CRYHOME");
        CRYSetup *setup=new CRYSetup(setupString,Form("%s/data",crydata.data()));
        
        gen = new CRYGenerator(setup);
        
        CLHEP::HepRandom::setTheSeed(input->GetUserSEED());
        // set random number generator
    RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(),&CLHEP::HepRandomEngine::flat);
        setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
        InputState=0;
#else
        G4cout << "PrimaryGeneratorAction: CRY not included in build " << G4endl;
        InputState=-1;
#endif
    }
    
    const G4d2oNeutrinoAlley *neutrinoAlley = dynamic_cast<const G4d2oNeutrinoAlley*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    zPlaneOfCosmics = neutrinoAlley->GetTopOfConcreteCeiling();
    std::cout<<"Getting Neutrino Alley TopOfCeiling " <<zPlaneOfCosmics<<std::endl;
    
    etfTimer = new ReplayTools();
    etfTimer->PrepareETFTimer(5000, input->GetNumberOfEvents()); //time in ms
    
    G4cout << "done." << G4endl;

}//END of constructor

G4d2oCosmicGun::~G4d2oCosmicGun()
{
    
	G4cout<<"Deleting G4d2oCosmicGun...";
	
//    delete particleGun;
	
	G4cout<<"done."<<G4endl;
	
}//END of destructor

void G4d2oCosmicGun::GeneratePrimaries(G4Event* anEvent)
{
    
    theEventNum = anEvent->GetEventID();
    etfTimer->SetCurrentEvent(theEventNum);
    gSystem->ProcessEvents();
    
    if(anEvent->GetEventID()==0 && input->GetPrintStatus()!=3) etfTimer->StartUpdateTimer();


    if (InputState != 0) {
        G4String* str = new G4String("CRY library was not successfully initialized");
        //G4Exception(*str);
        G4Exception("PrimaryGeneratorAction", "1",
                    RunMustBeAborted, *str);
    }
    
#ifdef USE_CRY

    G4String particleName;
    vect.clear();
    gen->genEvent(&vect);
    
    //....debug output

    if (anEvent->GetEventID() % 1000 == 0) {
    G4cout << "\nEvent=" << anEvent->GetEventID() << " "
    << "CRY generated nparticles=" << vect.size()
    << "Time=" << gen->timeSimulated()
    << G4endl;
    }

G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

for (unsigned j = 0; j < vect.size(); ++j) {
    auto* cp = vect[j];

    // Definition
    particleGun->SetParticleDefinition(
        particleTable->FindParticle(cp->PDGid())
    );

    // Kinetic energy (CRY in MeV)
    particleGun->SetParticleEnergy(cp->ke() * MeV);

    // Direction (normalize for safety)
    G4ThreeVector dir(cp->u(), cp->v(), cp->w());
    particleGun->SetParticleMomentumDirection(dir.unit());

    // Position:
    // Option A (recommended): use CRY’s plane coordinates directly (cm -> Geant4)
//    G4ThreeVector pos(cp->x() * cm, cp->y() * cm, cp->z() * cm);

    // Option B (only if you deliberately use a constant plane):
    G4ThreeVector pos(cp->x() * m, cp->y() * m, zPlaneOfCosmics /* in G4 units, e.g., 20*m */);

    particleGun->SetParticlePosition(pos);

    // Time (CRY in ns)
    particleGun->SetParticleTime(cp->t() * ns);

    // Add to event
    particleGun->GeneratePrimaryVertex(anEvent);

    delete cp;
}




 /***
    for ( unsigned j=0; j<vect.size(); j++) {
        particleName=CRYUtils::partName(vect[j]->id());
        
        //....debug output
        if (anEvent->GetEventID() % 1000 == 0) {
        cout << "  "          << particleName << " "
        << "charge="      << vect[j]->charge() << " "
        << setprecision(4)
        << "energy (MeV)=" << vect[j]->ke()*MeV << " "
        << "pos (m)"
        //<< G4ThreeVector(vect[j]->x(), vect[j]->y(), vect[j]->z())
        << G4ThreeVector(vect[j]->x(), vect[j]->y(), zPlaneOfCosmics)
        << " " << "direction cosines "
        << G4ThreeVector(vect[j]->u(), vect[j]->v(), vect[j]->w())
        << " " << endl;
        }
        
        particleGun->SetParticleDefinition(particleTable->FindParticle(vect[j]->PDGid()));
        particleGun->SetParticleEnergy(vect[j]->ke()*MeV);
        particleGun->SetParticlePosition(G4ThreeVector(vect[j]->x()*m, vect[j]->y()*m, zPlaneOfCosmics));
        particleGun->SetParticleMomentumDirection(G4ThreeVector(vect[j]->u(), vect[j]->v(), vect[j]->w()));
        particleGun->SetParticleTime(vect[j]->t());
        particleGun->GeneratePrimaryVertex(anEvent);
        delete vect[j];

    } ***/

#endif

}//END of GeneratePrimaries()

