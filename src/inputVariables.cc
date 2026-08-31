#include <fstream>
#include <iostream>
#include <ctype.h>

#include "TMath.h"

#include "inputVariables.hh"

inputVariables * inputVariables::ptrIV = NULL;
ofstream * inputVariables::ptrOF = NULL;

inputVariables::inputVariables(int argc, char** argv)
{
    G4cerr << "Instance of inputVariables Constructed!" << G4endl;

    currentVals = new GenericInputHolder();
    newVals = new GenericInputHolder();
    thePTracFiles = new TObjArray();

    ReadCurrentFile();

    if(argc==2 && (strcmp(argv[1],"--help")==0 || strcmp(argv[1], "help")==0)){
        PrintHelp(argv[0]);
        exit(0);
    }
    else if(argc==2 && strcmp(argv[1],"ShowCurrent")==0){
        printf("\nCurrent beamOn.dat\n\n");
        WriteValues(currentVals,false);
        printf("\n\n");
        exit(0);
    }
    else if(argc>1 && strcmp(argv[1],"Setup")==0){
        printf("\nSetting up beamOn.dat\n\n");
        Bool_t bUpdate = UpdateValues(argc, argv);
        isave = 0;
        if(bUpdate){
            SetValues(newVals);
            SetupDirectory();
            WriteValues(newVals);
            WriteValues(newVals,false);
        }
        printf("\n\n");
        exit(0);
    }
    else if(argc>1 && strcmp(argv[1],"CommandLine")==0){
        PrintCommandLine(argv,currentVals);
        exit(0);
    }

    //Setting values here (no change to beamOn.dat)
    if(argc==1) SetValues(currentVals); //use values from the file
    else{
        Bool_t bUpdate = UpdateValues(argc, argv);
        if(!bUpdate){
            exit(0);
        }
        SetValues(newVals);
    }
    SetupDirectory();

    //Writing a new beamOn.dat file
    WriteValues(newVals,true,Form("%s/beamOnFiles/beamOn%03d.dat",outputDir.Data(),runno));

    dateTime = new TDatime();

    numHitsColl = 0;
    numDetHC = 0;
    numSrcHC = 0;
    numTargHC = 0;
}

inputVariables::~inputVariables()
{
    G4cerr<<"Deleting inputVariables...";
    delete currentVals;
    delete newVals;
    G4cerr<<"done"<<G4endl;
}

void inputVariables::SetupDirectory(){
    system(Form("mkdir -p %s",outputDir.Data()));
    system(Form("mkdir -p %s/fileHistories",outputDir.Data()));
    system(Form("mkdir -p %s/beamOnFiles",outputDir.Data()));
}

inputVariables* inputVariables::GetIVPointer( int argc, char** argv )
{
    if(!ptrIV) ptrIV = new inputVariables(argc, argv);
    return ptrIV;
}

ofstream * inputVariables::GetOutputFile( void )
{
    if( ptrOF==NULL ){
        dateTime->Set();
        ptrOF = new ofstream( "tracking_hits.dat", ios::out );
        *ptrOF <<dateTime->GetMonth()<<"/"<<dateTime->GetDay()<<"/"<<dateTime->GetYear()<<"   "<<dateTime->GetHour();
        if(dateTime->GetMinute() < 10) *ptrOF<<":0";
        if(dateTime->GetMinute() >= 10) *ptrOF<<":";
        *ptrOF<<dateTime->GetMinute()<<G4endl;
        *ptrOF <<" Run Number: "<<runno<<G4endl;
    }
    return ptrOF;
}

Bool_t inputVariables::ReadCurrentFile(){
    ifstream infile("beamOn.dat");

    for(Int_t i=0; i<NUMINT; i++)
        infile >> currentVals->iPar[i] >> currentVals->iDesc[i];
    for(Int_t i=0; i<NUMLONG; i++)
        infile >> currentVals->lPar[i] >> currentVals->lDesc[i];
    for(Int_t i=0; i<NUMDOUBLE; i++)
        infile >> currentVals->dPar[i] >> currentVals->dDesc[i];
    for(Int_t i=0; i<NUMCHAR; i++)
        infile >> currentVals->sPar[i] >> currentVals->sDesc[i];

    infile.close();

    ifstream infile0("beamOn0.dat");
    if(!infile0){
        printf("\tError in inputVariables::ReadCurrentFile() - can't find the template input file: beamOn0.dat.\n\n");
        exit(0);
    }

    TString sDum;
    for(Int_t i=0; i<NUMINT; i++)
        infile0 >> sDum >> currentVals->iDesc[i];
    for(Int_t i=0; i<NUMLONG; i++)
        infile0 >> sDum >> currentVals->lDesc[i];
    for(Int_t i=0; i<NUMDOUBLE; i++)
        infile0 >> sDum >> currentVals->dDesc[i];
    for(Int_t i=0; i<NUMCHAR; i++)
        infile0 >> sDum >> currentVals->sDesc[i];

    infile0.close();

    for(Int_t i=0; i<NUMINT; i++){
        newVals->iPar[i] = currentVals->iPar[i];
        sprintf(newVals->iDesc[i],"%s",currentVals->iDesc[i]);
    }
    for(Int_t i=0; i<NUMLONG; i++){
        newVals->lPar[i] = currentVals->lPar[i];
        sprintf(newVals->lDesc[i],"%s",currentVals->lDesc[i]);
    }
    for(Int_t i=0; i<NUMDOUBLE; i++){
        newVals->dPar[i] = currentVals->dPar[i];
        sprintf(newVals->dDesc[i],"%s",currentVals->dDesc[i]);
    }
    for(Int_t i=0; i<NUMCHAR; i++){
        newVals->sPar[i] = currentVals->sPar[i];
        sprintf(newVals->sDesc[i],"%s",currentVals->sDesc[i]);
    }

    return true;
}

void inputVariables::WriteValues(GenericInputHolder *theVals, Bool_t bFile, TString sBOFileName){
    ofstream newFile;
    if(bFile) newFile.open(sBOFileName);

    for(Int_t i=0; i<NUMINT; i++){
        if(!bFile) cout<<theVals->iPar[i]<<"\t\t"<<theVals->iDesc[i]<<endl;
        else{
            if(i==0 && isave==1) newFile<<theVals->iPar[i]+1<<"\t\t"<<theVals->iDesc[i]<<endl;
            else newFile<<theVals->iPar[i]<<"\t\t"<<theVals->iDesc[i]<<endl;
        }
    }

    for(Int_t i=0; i<NUMLONG; i++){
        if(!bFile) cout<<theVals->lPar[i]<<"\t\t"<<theVals->lDesc[i]<<endl;
        else newFile<<theVals->lPar[i]<<"\t\t"<<theVals->lDesc[i]<<endl;
    }

    for(Int_t i=0; i<NUMDOUBLE; i++){
        if(!bFile) cout<<theVals->dPar[i]<<"\t\t"<<theVals->dDesc[i]<<endl;
        else newFile<<theVals->dPar[i]<<"\t\t"<<theVals->dDesc[i]<<endl;
    }

    for(Int_t i=0; i<NUMCHAR; i++){
        if(!bFile) cout<<theVals->sPar[i]<<"\t\t"<<theVals->sDesc[i]<<endl;
        else newFile<<theVals->sPar[i]<<"\t\t"<<theVals->sDesc[i]<<endl;
    }
}

void inputVariables::SetValues(GenericInputHolder *theVals){
    runno               = theVals->iPar[indexRunNum];
    nevents             = theVals->iPar[indexNumEvents];
    ivis                = theVals->iPar[indexVis];
    iprint              = theVals->iPar[indexDebug];
    irand               = theVals->iPar[indexSeed];
    iPhysics            = theVals->iPar[indexPhysics];
    iNeutronHP          = theVals->iPar[indexNeutronHP];
    iPGA                = theVals->iPar[indexPGA];
    iSideLining         = theVals->iPar[indexSideLining];
    iUseBottomPMTs      = theVals->iPar[indexUseBottomPMTs];
    iPMTQE              = theVals->iPar[indexPMTQE];
    iUseBottomVeto      = theVals->iPar[indexUseBottomVeto];
    iUseBottomShielding = theVals->iPar[indexUseBottomShielding];

    pmtDiameter         = theVals->dPar[indexPMTDiam];
    tailThick           = theVals->dPar[indexTailThick];
    shieldThick         = theVals->dPar[indexShieldThick];
    h2oRefl             = theVals->dPar[indexH2oRefl];
    reflectivity        = theVals->dPar[indexReflectivity];
    specialEnergy       = theVals->dPar[indexSpecialEnergy];

    userSEED            = theVals->lPar[indexUserSEED];
    outputDir           = theVals->sPar[indexOutDir];
}

void inputVariables::PrintCommandLine(char** argv, GenericInputHolder *theVals){
    printf("\n");
    printf("Shell command:\n");
    printf("%s %d ",argv[0],theVals->iPar[indexRunNum]);

    printf("Events %d ",theVals->iPar[indexNumEvents]);

    if(theVals->iPar[indexVis]==0) printf("NoVis ");
    else if(theVals->iPar[indexVis]==1) printf("OpenGL ");
    else if(theVals->iPar[indexVis]==2) printf("Wired ");
    else if(theVals->iPar[indexVis]==3) printf("Dawn ");

    if(theVals->iPar[indexDebug]==0) printf("NoTracks ");
    else if(theVals->iPar[indexDebug]==1) printf("TrackScreen ");
    else if(theVals->iPar[indexDebug]==2) printf("TrackFile ");

    if(theVals->iPar[indexSeed]==0) printf("FixSeed ");
    else if(theVals->iPar[indexSeed]==1) printf("ClockSeed ");

    if(theVals->iPar[indexPhysics]==0) printf("PhysicsOff ");
    else if(theVals->iPar[indexPhysics]==1) printf("PhysicsOn ");

    if(theVals->iPar[indexNeutronHP]==0) printf("NoNeutronHP ");
    else if(theVals->iPar[indexNeutronHP]==1) printf("NeutronHP ");

    if(theVals->iPar[indexPGA]==0) printf("PhotonPGA ");
    else if(theVals->iPar[indexPGA]==1) printf("ElectronPGA ");
    else if(theVals->iPar[indexPGA]==2) printf("CosmicPGA ");
    else if(theVals->iPar[indexPGA]==3) printf("MichelElectronPGA ");

    if(theVals->iPar[indexSideLining]==0) printf("SideTeflon ");
    else if(theVals->iPar[indexSideLining]==1) printf("SidePMTs ");

    if(theVals->iPar[indexUseBottomPMTs]==0) printf("BottomTeflon ");
    else if(theVals->iPar[indexUseBottomPMTs]==1) printf("BottomPMTs ");

    if(theVals->iPar[indexPMTQE]==0) printf("StandardQE ");
    else if(theVals->iPar[indexPMTQE]==1) printf("HighQE ");

    if(theVals->iPar[indexUseBottomVeto]==0) printf("NoBottomVeto ");
    else if(theVals->iPar[indexUseBottomVeto]==1) printf("UseBottomVeto ");

    if(theVals->iPar[indexUseBottomShielding]==0) printf("NoBottomShielding ");
    else if(theVals->iPar[indexUseBottomShielding]==1) printf("UseBottomShielding ");

    printf("PMTd %g ",theVals->dPar[indexPMTDiam]);
    printf("TailThick %g ",theVals->dPar[indexTailThick]);
    printf("ShieldThick %g ",theVals->dPar[indexShieldThick]);
    printf("H2oRefl %g ",theVals->dPar[indexH2oRefl]);
    printf("Reflectivity %g ",theVals->dPar[indexReflectivity]);
    printf("SpecialEnergy %g ",theVals->dPar[indexSpecialEnergy]);

    printf("UserSEED %lld ",theVals->lPar[indexUserSEED]);
    printf("OutputDir %s ",theVals->sPar[indexOutDir].Data());
    printf("\n\n");
}

Bool_t inputVariables::UpdateValues(int argc, char** argv){
    if(argc<2) return true;

    Int_t iStart = 2;

    if(argc>1){
        Bool_t bDigit = false;
        if(isdigit(argv[1][0])) bDigit = true;
        else if( ((argv[1][0]=='+') || (argv[1][0]=='-')) && isdigit(argv[1][1]) ) bDigit = true;

        Bool_t bSetup = false;
        if( strcmp(argv[1], "Setup")==0) bSetup = true;

        if(!bDigit && !bSetup){
            printf("\nError in inputVariables::UpdateValues() - first argument \"%s\" must be:\n",argv[1]);
            printf("\ta run number, \"ShowCurrent\", \"Setup\", or \"help\". \n\n");
            return false;
        }

        if(!bSetup){
            if(argv[1][0]=='+' || argv[1][0]=='-' || argv[1][0]=='0'){
                newVals->iPar[indexRunNum] += atoi(argv[1]);
            }
            else newVals->iPar[indexRunNum] = atoi(argv[1]);
        }
        else if(argc>2){
            if(argv[2][0]=='+' || argv[2][0]=='-' || argv[2][0]=='0'){
                newVals->iPar[indexRunNum] += atoi(argv[2]);
            }
            else newVals->iPar[indexRunNum] = atoi(argv[2]);
            iStart++;
        }
    }

    if(argc>iStart){
        for(Int_t iArg=iStart; iArg<argc; iArg++){
            if(!strcmp(argv[iArg], "EVENTS") || !strcmp(argv[iArg], "Events") || !strcmp(argv[iArg], "events")){
                if(iArg+1 < argc){
                    newVals->iPar[indexNumEvents] = atoi(argv[iArg+1]);
                    iArg++;
                }
                else{
                    printf("Error in inputVariables::UpdateValues() - number of events must follow command \"%s\"\n",argv[iArg]);
                    return false;
                }
            }
            else if(!strcmp(argv[iArg], "NOVIS") || !strcmp(argv[iArg], "NoVis") || !strcmp(argv[iArg], "novis"))
                newVals->iPar[indexVis] = 0;
            else if(!strcmp(argv[iArg], "WIRED") || !strcmp(argv[iArg], "Wired") || !strcmp(argv[iArg], "wired"))
                newVals->iPar[indexVis] = 2;
            else if(!strcmp(argv[iArg], "OPENGL") || !strcmp(argv[iArg], "OpenGL") || !strcmp(argv[iArg], "opengl"))
                newVals->iPar[indexVis] = 1;
            else if(!strcmp(argv[iArg], "DAWN") || !strcmp(argv[iArg], "Dawn") || !strcmp(argv[iArg], "dawn"))
                newVals->iPar[indexVis] = 3;
            else if(!strcmp(argv[iArg], "NOTRACKS") || !strcmp(argv[iArg], "NoTracks") || !strcmp(argv[iArg], "notracks"))
                newVals->iPar[indexDebug] = 0;
            else if(!strcmp(argv[iArg], "TRACKFILE") || !strcmp(argv[iArg], "TrackFile") || !strcmp(argv[iArg], "trackfile"))
                newVals->iPar[indexDebug] = 2;
            else if(!strcmp(argv[iArg], "TRACKSCREEN") || !strcmp(argv[iArg], "TrackScreen") || !strcmp(argv[iArg], "trackScreen"))
                newVals->iPar[indexDebug] = 1;
            else if(!strcmp(argv[iArg], "FIXSEED") || !strcmp(argv[iArg], "FixSeed") || !strcmp(argv[iArg], "fixseed"))
                newVals->iPar[indexSeed] = 0;
            else if(!strcmp(argv[iArg], "CLOCKSEED") || !strcmp(argv[iArg], "ClockSeed") || !strcmp(argv[iArg], "clockseed"))
                newVals->iPar[indexSeed] = 1;
            else if(!strcmp(argv[iArg], "PHOTONPGA") || !strcmp(argv[iArg], "PhotonPGA") || !strcmp(argv[iArg], "photonPGA"))
                newVals->iPar[indexPGA] = 0;
            else if(!strcmp(argv[iArg], "ELECTRONPGA") || !strcmp(argv[iArg], "ElectronPGA") || !strcmp(argv[iArg], "electronPGA"))
                newVals->iPar[indexPGA] = 1;
            else if(!strcmp(argv[iArg], "COSMICPGA") || !strcmp(argv[iArg], "CosmicPGA") || !strcmp(argv[iArg], "cosmicPGA"))
                newVals->iPar[indexPGA] = 2;
            else if(!strcmp(argv[iArg], "MICHELELECTRONPGA") || !strcmp(argv[iArg], "MichelElectronPGA") || !strcmp(argv[iArg], "michelelectronPGA"))
                newVals->iPar[indexPGA] = 3;
            else if(!strcmp(argv[iArg], "PHYSICSON") || !strcmp(argv[iArg], "PhysicsOn") || !strcmp(argv[iArg], "physicson"))
                newVals->iPar[indexPhysics] = 1;
            else if(!strcmp(argv[iArg], "PHYSICSOFF") || !strcmp(argv[iArg], "PhysicsOff") || !strcmp(argv[iArg], "physicsoff"))
                newVals->iPar[indexPhysics] = 0;
            else if(!strcmp(argv[iArg], "NEUTRONHP") || !strcmp(argv[iArg], "NeutronHP") || !strcmp(argv[iArg], "neutronHP"))
                newVals->iPar[indexNeutronHP] = 1;
            else if(!strcmp(argv[iArg], "NONEUTRONHP") || !strcmp(argv[iArg], "NoNeutronHP") || !strcmp(argv[iArg], "noneutronHP"))
                newVals->iPar[indexNeutronHP] = 0;
            else if(!strcmp(argv[iArg], "SIDETEFLON") || !strcmp(argv[iArg], "SideTeflon") || !strcmp(argv[iArg], "sideTeflon"))
                newVals->iPar[indexSideLining] = 0;
            else if(!strcmp(argv[iArg], "SIDEPMTS") || !strcmp(argv[iArg], "SidePMTs") || !strcmp(argv[iArg], "sidePMTs") ||
                    !strcmp(argv[iArg], "SIDEPMT") || !strcmp(argv[iArg], "SidePMT") || !strcmp(argv[iArg], "sidePMT") ||
                    !strcmp(argv[iArg], "sidepmt") || !strcmp(argv[iArg], "sidepmts"))
                newVals->iPar[indexSideLining] = 1;
            else if(!strcmp(argv[iArg], "BOTTOMTEFLON") || !strcmp(argv[iArg], "BottomTeflon") || !strcmp(argv[iArg], "bottomTeflon"))
                newVals->iPar[indexUseBottomPMTs] = 0;
            else if(!strcmp(argv[iArg], "BOTTOMPMTS") || !strcmp(argv[iArg], "BottomPMTs") || !strcmp(argv[iArg], "bottomPMTs") ||
                    !strcmp(argv[iArg], "BOTTOMPMT") || !strcmp(argv[iArg], "BottomPMT") || !strcmp(argv[iArg], "bottomPMT") ||
                    !strcmp(argv[iArg], "bottompmt") || !strcmp(argv[iArg], "bottompmts"))
                newVals->iPar[indexUseBottomPMTs] = 1;
            else if(!strcmp(argv[iArg], "STANDARDQE") || !strcmp(argv[iArg], "StandardQE") || !strcmp(argv[iArg], "standardQE"))
                newVals->iPar[indexPMTQE] = 0;
            else if(!strcmp(argv[iArg], "HIGHQE") || !strcmp(argv[iArg], "HighQE") || !strcmp(argv[iArg], "highQE"))
                newVals->iPar[indexPMTQE] = 1;
            else if(!strcmp(argv[iArg], "NOBOTTOMVETO") || !strcmp(argv[iArg], "NoBottomVeto") || !strcmp(argv[iArg], "noBottomVeto"))
                newVals->iPar[indexUseBottomVeto] = 0;
            else if(!strcmp(argv[iArg], "USEBOTTOMVETO") || !strcmp(argv[iArg], "UseBottomVeto") || !strcmp(argv[iArg], "useBottomVeto"))
                newVals->iPar[indexUseBottomVeto] = 1;
            else if(!strcmp(argv[iArg], "NOBOTTOMSHIELDING") || !strcmp(argv[iArg], "NoBottomShielding") || !strcmp(argv[iArg], "noBottomShielding"))
                newVals->iPar[indexUseBottomShielding] = 0;
            else if(!strcmp(argv[iArg], "USEBOTTOMSHIELDING") || !strcmp(argv[iArg], "UseBottomShielding") || !strcmp(argv[iArg], "useBottomShielding"))
                newVals->iPar[indexUseBottomShielding] = 1;
            else if(!strcmp(argv[iArg], "PMTd") || !strcmp(argv[iArg], "pmtd") || !strcmp(argv[iArg], "PMTdiameter") || !strcmp(argv[iArg], "pmtdiameter")){
                if(iArg+1 < argc){
                    newVals->dPar[indexPMTDiam] = strtod(argv[iArg+1],NULL);
                    iArg++;
                }
                else{
                    printf("Error in inputVariables::UpdateValues() - PMT diameter in inches must follow command \"%s\"\n",argv[iArg]);
                    return false;
                }
            }
            else if(!strcmp(argv[iArg], "TailThick") || !strcmp(argv[iArg], "tailthick") || !strcmp(argv[iArg], "Tailthick") || !strcmp(argv[iArg], "tailThick")) {
                if(iArg+1 < argc) {
                    newVals->dPar[indexTailThick] = strtod(argv[iArg+1],NULL);
                    iArg++;
                }
                else{
                    printf("Error in inputVariables::UpdateValues() - Tail Catcher thickness in cm must follow command \"%s\"\n",argv[iArg]);
                    return false;
                }
            }
            else if(!strcmp(argv[iArg], "ShieldThick") || !strcmp(argv[iArg], "shieldthick") || !strcmp(argv[iArg], "Shieldthick") || !strcmp(argv[iArg], "shieldThick")) {
                if(iArg+1 < argc) {
                    newVals->dPar[indexShieldThick] = strtod(argv[iArg+1],NULL);
                    iArg++;
                }
                else{
                    printf("Error in inputVariables::UpdateValues() - Shield thickness in cm must follow command \"%s\"\n",argv[iArg]);
                    return false;
                }
            }
            else if(!strcmp(argv[iArg], "H2oRefl") || !strcmp(argv[iArg], "h2orefl") || !strcmp(argv[iArg], "H2orefl") || !strcmp(argv[iArg], "h2oRefl")) {
                if(iArg+1 < argc) {
                    newVals->dPar[indexH2oRefl] = strtod(argv[iArg+1],NULL);
                    iArg++;
                }
                else{
                    printf("Error in inputVariables::UpdateValues() - Refl ofh2o must follow command \"%s\"\n",argv[iArg]);
                    return false;
                }
            }
            else if(!strcmp(argv[iArg], "Reflectivity") || !strcmp(argv[iArg], "reflectivity") ) {
                if(iArg+1 < argc) {
                    newVals->dPar[indexReflectivity] = strtod(argv[iArg+1],NULL);
                    iArg++;
                }
                else{
                    printf("Error in inputVariables::UpdateValues() - Reflectivity must follow command \"%s\"\n",argv[iArg]);
                    return false;
                }
            }
            else if(!strcmp(argv[iArg], "SpecialEnergy") || !strcmp(argv[iArg], "specialEnergy") ) {
                if(iArg+1 < argc) {
                    newVals->dPar[indexSpecialEnergy] = strtod(argv[iArg+1],NULL);
                    iArg++;
                }
                else{
                    printf("Error in inputVariables::UpdateValues() - SpecialEnergy must follow command \"%s\"\n",argv[iArg]);
                    return false;
                }
            }
            else if(!strcmp(argv[iArg], "UserSEED") || !strcmp(argv[iArg], "UserSeed") || !strcmp(argv[iArg], "userSeed") || !strcmp(argv[iArg], "userseed")){
                if(iArg+1 < argc){
                    newVals->lPar[indexUserSEED] = (G4long)(strtoul(argv[iArg+1],NULL,10)&0x7FFFFFFFFFFFFFFF);
                    iArg++;
                }
                else{
                    printf("Error in inputVariables::UpdateValues() - User SEED must follow command \"%s\"\n",argv[iArg]);
                    return false;
                }
            }
            else if(!strcmp(argv[iArg], "OUTPUTDIR") || !strcmp(argv[iArg], "OutputDir") || !strcmp(argv[iArg], "outputdir") || !strcmp(argv[iArg], "outputDir") ){
                if(iArg+1 < argc){
                    newVals->sPar[indexOutDir] = argv[iArg+1];
                    iArg++;
                }
                else{
                    printf("Error in inputVariables::UpdateValues() - output directory name must follow command \"%s\"\n",argv[iArg]);
                    return false;
                }
            }
            else{
                printf("Error in inputVariables::UpdateValues() - Unrecognized command \"%s\"\n",argv[iArg]);
                return false;
            }
        }
    }
    return true;
}

void inputVariables::PrintHelp(Char_t *exeName){
    printf("\n\tNo arguments executes with no changes to current beamOn.dat:\n");
    printf("\t\t%s \n",exeName);
    printf("\n\tCommand line usage with absolute run number specified:\n");
    printf("\t\t%s [ N ]  [options (below)]       # (N = absolute run number)\n",exeName);
    printf("\n\tCommand line usage with relative run number specified:\n");
    printf("\t\t%s [ +n ] [options (below)]       # (n = current run number incremented by n)\n",exeName);
    printf("\t\t%s [ -n ] [options (below)]       # (n = current run number reduced by n)\n",exeName);
    printf("\n\tTo see the command line argument corresponding to the current setup:\n");
    printf("\t\t%s CommandLine\n",exeName);

    cout<<"\n\tOptions:"<<endl;
    cout<<"\t\t[EVENTS #]                                                 "<<currentVals->iDesc[indexNumEvents]<<endl;
    cout<<"\t\t[NoVis | OPENGL | WIRED | DAWN]                            "<<currentVals->iDesc[indexVis]<<endl;
    cout<<"\t\t[NoTracks | TrackFile | TrackScreen ]                      "<<currentVals->iDesc[indexDebug]<<endl;
    cout<<"\t\t[FixSeed | ClockSeed]                                      "<<currentVals->iDesc[indexSeed]<<endl;
    cout<<"\t\t[PhotonPGA | ElectronPGA | CosmicPGA | MichelElectronPGA ] "<<currentVals->iDesc[indexPGA]<<endl;
    cout<<"\t\t[PhysicsOn | PhysicsOff]                                   "<<currentVals->iDesc[indexPhysics]<<endl;
    cout<<"\t\t[NoNeutronHP | NeutronHP]                                  "<<currentVals->iDesc[indexNeutronHP]<<endl;
    cout<<"\t\t[SideTeflon | SidePMTs]                                    "<<currentVals->iDesc[indexSideLining]<<endl;
    cout<<"\t\t[BottomTeflon | BottomPMTs]                                "<<currentVals->iDesc[indexUseBottomPMTs]<<endl;
    cout<<"\t\t[StandardQE | HighQD]                                      "<<currentVals->iDesc[indexPMTQE]<<endl;
    cout<<"\t\t[NoBottomVeto | UseBottomVeto]                             "<<currentVals->iDesc[indexUseBottomVeto]<<endl;
    cout<<"\t\t[NoBottomShielding | UseBottomShielding]                   "<<currentVals->iDesc[indexUseBottomShielding]<<endl;
    cout<<"\t\t[PMTd # ]                                                  "<<currentVals->dDesc[indexPMTDiam]<<endl;
    cout<<"\t\t[TailThick # ]                                             "<<currentVals->dDesc[indexTailThick]<<endl;
    cout<<"\t\t[ShieldThick # ]                                           "<<currentVals->dDesc[indexShieldThick]<<endl;
    cout<<"\t\t[H2oRefl # ]                                               "<<currentVals->dDesc[indexH2oRefl]<<endl;
    cout<<"\t\t[Reflectivity # ]                                          "<<currentVals->dDesc[indexReflectivity]<<endl;
    cout<<"\t\t[SpecialEnergy # ]                                         "<<currentVals->dDesc[indexSpecialEnergy]<<endl;
    cout<<"\t\t[UserSeed # ]                                              "<<currentVals->lDesc[indexUserSEED]<<endl;
    cout<<"\t\t[OutputDir <dirPath>]                                      "<<currentVals->sDesc[indexOutDir]<<endl;

    printf("\n\n");
}

G4int inputVariables::GetNewCollectionID(G4int iType)
{
    G4int iReturn = GetNumHitsCollections();
    numHitsColl++;
    if(iType==0){
        if(numDetHC<MAX_SEN_DET-1){
            detCollID[numDetHC] = iReturn;
            numDetHC++;
        }
        else{
            printf("Error in inputVariables::GetNewCollectionID() - too many sensitive detectors requested (%d).\n\tMax allowed is %d\n\n",numDetHC+1,MAX_SEN_DET);
            exit(0);
        }
    }
    else if(iType==1){
        if(numSrcHC<MAX_SEN_DET-1){
            srcCollID[numSrcHC] = iReturn;
            numSrcHC++;
        }
        else{
            printf("Error in inputVariables::GetNewCollectionID() - too many sensitive sources requested (%d).\n\tMax allowed is %d\n\n",numSrcHC+1,MAX_SEN_DET);
            exit(0);
        }
    }
    else if(iType==2){
        if(numTargHC<MAX_SEN_DET-1){
            targCollID[numTargHC] = iReturn;
            numTargHC++;
        }
        else{
            printf("Error in inputVariables::GetNewCollectionID() - too many sensitive targets requested (%d).\n\tMax allowed is %d\n\n",numTargHC+1,MAX_SEN_DET);
            exit(0);
        }
    }
    return iReturn;
}

