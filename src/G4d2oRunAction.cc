//#include "TreeMaker.hh"

#include "G4d2oRunAction.hh"
#include "G4d2oEventAction.hh"
#include "G4d2oNeutrinoAlley.hh"
#include "G4d2oMaterialsDefinition.hh"
#include "G4d2oPrimaryGeneratorAction.hh"

#include "G4RunManager.hh"

G4d2oMaterialsDefinition * G4d2oRunAction::materialsPtr = NULL;

G4d2oRunAction::G4d2oRunAction()
{
    G4cerr << "\tConstructing G4d2oRunAction..." ;

    input = inputVariables::GetIVPointer();
    runno = input->GetRunNumber();
    numPrimaryEvents = input->GetNumberOfEvents();

    TDatime *dateTime = new TDatime();
    dateTime->Set();
    runDate.month      = dateTime->GetMonth();
    runDate.day        = dateTime->GetDay();
    runDate.year       = dateTime->GetYear();

    TString theDirPath = input->GetOutputDirectory();

    if(theDirPath == TString("")) sprintf(path, "%s","./data/");
    else sprintf(path, "%s/", theDirPath.Data());

    sprintf(name, "%s%03d", "Sim_D2ODetector",runno);

    outFileName = Form("%s/%s.root",path,name);

    system(Form("tar -czf %s/fileHistories/run%03d.tgz src/*.cc include/*.hh *.cc G4d2o.cc mac simEvent -C %s/beamOnFiles beamOn%03d.dat",path,runno,path,runno));
    mBeamOn = TMacro(Form("%s/beamOnFiles/beamOn%03d.dat",path,runno));

    theRT = (ReplayTools*)((G4d2oPrimaryGeneratorAction*)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction())->GetReplayTool();

    G4cerr << "done." << G4endl;
}

G4d2oRunAction::~G4d2oRunAction()
{
    G4cerr<<"Deleting G4d2oRunAction...";
    G4cerr << "done." << G4endl;
}

void G4d2oRunAction::BeginOfRunAction(const G4Run *aRun)
{
    ((G4Run *)(aRun))->SetRunID(runno);
    G4cout<<"\nRun ID: "<<aRun->GetRunID()<<G4endl<<G4endl;

    G4d2oEventAction *evAct = (G4d2oEventAction*)G4RunManager::GetRunManager()->GetUserEventAction();
    evAct->OutFilePtr()->cd();
    mBeamOn.Write("beamOn.dat");

    setupTree = new TTree("Setup_Tree","Setup tree");
    setupTree->Branch( "numPrimaryEvents", &numPrimaryEvents, "numPrimaryEvents/I");
    setupTree->Branch( "runDate", &runDate, "month/i:day:year");
    setupTree->Branch( "totPMTs", &totPMTs, "totPMTs/I");
    setupTree->Branch( "numPMTRows", &pmtRows, "numPMTRows/I");
    setupTree->Branch( "tankSize", "TVector3", &tankSize);

    // Attach reflection tree to output file
    TFile* outFile = evAct->OutFilePtr();
    G4d2oMaterialsDefinition::AttachReflectionTree(outFile);

    G4cerr << "\t*********************" << G4endl;
    G4cerr << "\t***** BEGIN Run *****" << G4endl;
    G4cerr << "\t*********************" << G4endl;
    G4cerr << G4endl;

    theRT->DisplayProgress(-1);
}

void G4d2oRunAction::EndOfRunAction(const G4Run*)
{
    G4d2oNeutrinoAlley *detCon = (G4d2oNeutrinoAlley*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    G4d2oDetector *theDet = (G4d2oDetector*)detCon->GetDetectorPtr();

    totPMTs = theDet->GetTotalPMTS();
    pmtRows = theDet->GetPMTRows();
    tankSize.SetXYZ(theDet->TankX(),
                    theDet->TankY(),
                    theDet->TankZ());

    theRT->DisplayProgress(-2);

    G4cerr << G4endl;
    G4cerr << "\t*********************" << G4endl;
    G4cerr << "\t******* END Run *****" << G4endl;
    G4cerr << "\t*********************" << G4endl;

    G4d2oEventAction *evAct = (G4d2oEventAction*)G4RunManager::GetRunManager()->GetUserEventAction();
    evAct->OutFilePtr()->cd();
    setupTree->Fill();
    setupTree->AutoSave();

    // Close reflection tree
    G4d2oMaterialsDefinition::CloseReflectionTree();
}

G4d2oMaterialsDefinition* G4d2oRunAction::GetMaterialsPointer( void )
{
    if( materialsPtr==NULL )
        materialsPtr = new G4d2oMaterialsDefinition();

    return materialsPtr;
}
