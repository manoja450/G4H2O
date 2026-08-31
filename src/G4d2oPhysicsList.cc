#include "G4d2oPhysicsList.hh"
#include "G4d2oCustomOpBoundary.hh"
#include "G4d2oDataDrivenReflector.hh"

#include <fstream>

G4d2oPhysicsList::G4d2oPhysicsList(G4bool bOn, G4bool bNHP)
: G4VModularPhysicsList()
{
    defaultCutValue = 0.01*mm;

    bNull = !(bOn);
    bNeutronHP = bNHP;

    if(bNull) G4cerr<<"\tG4d2oPhysicsList constructor: Physics list is empty. \n\t\tPass the argument 'bOn=true' to turn it on."<<G4endl;
    else{
        if(!bNeutronHP) G4cerr<<"\tG4d2oPhysicsList constructor: Using non-HP neutron physics. \n\t\tPass the argument 'bNHP=true' to turn on HP neutron physics."<<G4endl;
    }

    verbosityLevel = 1;
    SetVerboseLevel(verbosityLevel);
}

G4d2oPhysicsList::~G4d2oPhysicsList() {
    G4cerr<<"Deleting G4d2oPhysicsList...";
    G4cerr<<"done"<<G4endl;
}

//---------------------------------------------------------
//Construct Particles
//---------------------------------------------------------

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

void G4d2oPhysicsList::ConstructParticle()
{
    G4cerr<<"\tConstructing particles...\n";

    G4cerr<<"\t\tbosons"<<G4endl;
    G4BosonConstructor bConstructor;
    bConstructor.ConstructParticle();

    G4cerr<<"\t\tleptons"<<G4endl;
    G4LeptonConstructor lConstructor;
    lConstructor.ConstructParticle();

    G4cerr<<"\t\tmesons"<<G4endl;
    G4MesonConstructor mConstructor;
    mConstructor.ConstructParticle();

    G4cerr<<"\t\tbaryons"<<G4endl;
    G4BaryonConstructor  hConstructor;
    hConstructor.ConstructParticle();

    G4cerr<<"\t\tlight ions"<<G4endl;
    G4IonConstructor iConstructor;
    iConstructor.ConstructParticle();

    G4cerr<<"\t\tshort lived"<<G4endl;
    G4ShortLivedConstructor sConstructor;
    sConstructor.ConstructParticle();

    G4cerr<<"\tDone constructing particles."<<G4endl;
}

//---------------------------------------------------------
//Construct Processes
//---------------------------------------------------------

#include "G4ProcessManager.hh"

void G4d2oPhysicsList::ConstructProcess()
{
    AddTransportation();

    if(!bNull){
        G4cerr << "G4d2oPhysicsList::ConstructProcess() -> Registering EM processes" << G4endl;
        ConstructEM();

        G4bool bTurnOnRAD = false;
        if(bTurnOnRAD){
            G4cerr << "G4d2oPhysicsList::ConstructProcess() -> Registering radioactive decay processes" << G4endl;
            ConstructRadioactiveDecay();
        }
        else{
            G4cerr << "G4d2oPhysicsList::ConstructProcess() -> Radioactive decay processes DISABLED" << G4endl;
        }

        G4cerr << "G4d2oPhysicsList::ConstructProcess() -> Registering hadronic processes" << G4endl;
        ConstructHadronic();

        G4cerr << "G4d2oPhysicsList::ConstructProcess() -> Registering optical processes" << G4endl;
        ConstructOptical();
    }

    G4cerr << "\tG4d2oPhysicsList::ConstructProcess() -> Done registering physics processes" << G4endl;
}

#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"

void G4d2oPhysicsList::ConstructEM()
{
    G4VPhysicsConstructor *emPhysicsList = new G4EmStandardPhysics(verbosityLevel);
    emPhysicsList->ConstructProcess();

    G4VPhysicsConstructor *emExtraPhysicsList = new G4EmExtraPhysics(verbosityLevel);
    emExtraPhysicsList->ConstructProcess();
}

#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

void G4d2oPhysicsList::ConstructRadioactiveDecay()
{
    G4VPhysicsConstructor *decPhysicsList = new G4DecayPhysics();
    decPhysicsList->ConstructProcess();

    G4VPhysicsConstructor *radDecPhysicsList = new G4RadioactiveDecayPhysics();
    radDecPhysicsList->ConstructProcess();
}

#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4IonPhysics.hh"

void G4d2oPhysicsList::ConstructHadronic()
{
    G4bool bBIC = true;

    G4VPhysicsConstructor *hadStopping = 0;
    G4VPhysicsConstructor *hadElastic = 0;
    G4VPhysicsConstructor *hadPhysics = 0;
    G4VPhysicsConstructor *ionPhysics = 0;

    if(bNeutronHP) hadElastic = new G4HadronElasticPhysicsHP(verbosityLevel);
    else hadElastic = new G4HadronElasticPhysics(verbosityLevel);

    if(bBIC){
        ionPhysics = new G4IonBinaryCascadePhysics(verbosityLevel);
        if(bNeutronHP) hadPhysics = new G4HadronPhysicsQGSP_BIC_HP(verbosityLevel);
        else hadPhysics = new G4HadronPhysicsQGSP_BIC(verbosityLevel);
    }
    else{
        ionPhysics = new G4IonPhysics(verbosityLevel);
        if(bNeutronHP) hadPhysics = new G4HadronPhysicsQGSP_BERT_HP(verbosityLevel);
        else hadPhysics = new G4HadronPhysicsQGSP_BERT(verbosityLevel);
    }

    hadStopping = new G4StoppingPhysics(verbosityLevel);

    if(hadElastic) hadElastic->ConstructProcess();
    if(hadPhysics) hadPhysics->ConstructProcess();
    if(hadStopping) hadStopping->ConstructProcess();
    if(ionPhysics) ionPhysics->ConstructProcess();
}

#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4Cerenkov.hh"

void G4d2oPhysicsList::ConstructOptical()
{
    G4Scintillation* theScintillationProcess = new G4Scintillation("Scintillation");
    theScintillationProcess->SetTrackSecondariesFirst(true);

    G4int MaxNumPhotons = 300;
    G4Cerenkov *theCerenkovProcess = new G4Cerenkov("Cerenkov");
    theCerenkovProcess->SetTrackSecondariesFirst(true);
    theCerenkovProcess->SetMaxBetaChangePerStep(10.0);
    theCerenkovProcess->SetMaxNumPhotonsPerStep(MaxNumPhotons);

    G4OpAbsorption* theAbsorptionProcess = new G4OpAbsorption();
    G4OpRayleigh* theRayleighScatteringProcess = new G4OpRayleigh();

    // ============================================================
    // BOUNDARY PROCESS: Data-Driven only (Unified model removed)
    // ============================================================
    G4d2oCustomOpBoundary* theBoundaryProcess = new G4d2oCustomOpBoundary();
    G4cout << "  Using DATA-DRIVEN boundary process" << G4endl;
    theBoundaryProcess->SetAzimuthalModel(G4d2oCustomOpBoundary::kGaussian35);
    G4cout << "  Azimuthal model: " << theBoundaryProcess->GetAzimuthalModelName() << G4endl;

    GetParticleIterator()->reset();

    while( (*GetParticleIterator())() ){
        G4ParticleDefinition* particle = GetParticleIterator()->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();

        if (theScintillationProcess->IsApplicable(*particle)) {
            pmanager->AddProcess(theScintillationProcess);
            pmanager->SetProcessOrderingToLast(theScintillationProcess, idxAtRest);
            pmanager->SetProcessOrderingToLast(theScintillationProcess, idxPostStep);
        }

        if (theCerenkovProcess->IsApplicable(*particle)) {
            pmanager->AddProcess(theCerenkovProcess);
            pmanager->SetProcessOrdering(theCerenkovProcess, idxPostStep);
        }

        if (particleName == "opticalphoton") {
            G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
            pmanager->AddDiscreteProcess(theAbsorptionProcess);
            pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
            pmanager->AddDiscreteProcess(theBoundaryProcess);
        }
    }
}

#include "G4UnitsTable.hh"

void G4d2oPhysicsList::SetCuts()
{
    if (verboseLevel >0)
    {
        G4cerr << "\tG4d2oPhysicsList::SetCuts:";
        G4cerr << "\tCutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
    }

    SetCutValue(defaultCutValue, "gamma");
    SetCutValue(defaultCutValue, "e-");
    SetCutValue(defaultCutValue, "e+");

    if (verboseLevel>0) DumpCutValuesTable();
    G4cerr<<"\tDone initializing G4d2oPhysicsList"<<G4endl;
}

//end of file G4d2oPhysicsList.cc
