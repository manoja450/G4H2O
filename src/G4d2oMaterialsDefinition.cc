#include "G4d2oMaterialsDefinition.hh"
#include "G4d2oDataDrivenReflector.hh"
#include "inputVariables.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "TString.h"
#include <fstream>
#include <limits>

G4d2oMaterialsDefinition::G4d2oMaterialsDefinition()
{
    input = inputVariables::GetIVPointer();
    h2oRefl = input->GetH2oRefl();
    manager = G4NistManager::Instance();
    
    matAir = manager->FindOrBuildMaterial("G4_AIR");
    G4double indexAir = 1.0;
    G4double absAir = 10000000.0*m;
    SetUniformOpticalProperties(matAir, indexAir, absAir);
   
    matVinylToluene = manager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    G4double indexVinylToluene = 1.58;
    G4double absVinylToluene = 0.0405*m;
    SetUniformOpticalProperties(matVinylToluene, indexVinylToluene, absVinylToluene);
    
    matDelrin = manager->FindOrBuildMaterial("G4_POLYOXYMETHYLENE");
    
    matAl=0;
    matSteel=0;
    matPb=0; 
    matPoly=0;
    matPMMA=0;
    matMuMetal=0;
    matCopper=0;
    matPlastic=0;
    matH2O=0;
    matD2O=0;
    matVacuum=0;
    matFusedSilica=0;
    matBorosilicate=0;
    matPhotoCathode=0;
    matTeflon=0;
    matTyvek=0;
    matEpoxy=0;
    matConcrete=0;
}

G4d2oMaterialsDefinition::~G4d2oMaterialsDefinition()
{
    G4cout<<"Deleting G4d2oMaterialsDefinition... done"<<G4endl;
}

G4Material * G4d2oMaterialsDefinition::GetMaterial( materialName matName )
{
    G4Material *theMaterial = NULL;
    G4double density, fractionmass, per[12];
    G4int ncomponents, natoms, z[12];
    
    switch( matName )
    {
        case AIR:
            theMaterial = matAir;
            break;
            
        case ALUMINUM:
            if(!matAl){
                density = 2.7*g/cm3;
                z[0] = 13;
                matAl = new G4Material("Aluminum", density, ncomponents=1);
                matAl->AddElement(manager->FindOrBuildElement(z[0]), fractionmass=1.0);
                printf("\tAluminum has been built.\n");
            }
            theMaterial = matAl;
            break;
          
        case VINYLTOLUENE:
            theMaterial = matVinylToluene;
            break;
            
        case LEAD:
            if(!matPb){
                density = 11.35*g/cm3;
                z[0] = 82;
                matPb = new G4Material("Lead", density, ncomponents=1);
                matPb->AddElement(manager->FindOrBuildElement(z[0]), fractionmass=1.0);
                printf("\tLead has been built.\n");
            }
            theMaterial = matPb;
            break;
            
        case POLY:
            if(!matPoly){
                density = 0.95*g/cm3;
                z[0] = 1;
                z[1] = 6;
                matPoly = new G4Material("Poly", density, ncomponents=2);
                matPoly->AddElement(manager->FindOrBuildElement(z[0]), natoms=2);
                matPoly->AddElement(manager->FindOrBuildElement(z[1]), natoms=1);
                printf("\tPolyethylene has been built.\n");
            }
            theMaterial = matPoly;
            break;

        case STEEL:
            if(!matSteel){
                per[0] = 0.08*perCent; z[0] = 6;
                per[1] = 2.0*perCent; z[1] = 25;
                per[2] = 18.0*perCent; z[2] = 24;
                per[3] = 1.0*perCent; z[3] = 42;
                per[4] = 10.0*perCent; z[4] = 28;
                per[5] = 1.0*perCent; z[5] = 14;
                per[6] = 0.045*perCent; z[6] = 15;
                per[7] = 0.03*perCent; z[7] = 16;
                per[8] = 100.0*perCent - (per[0]+per[1]+per[2]+per[3]+per[4]+per[5]+per[6]+per[7]); z[8] = 26;
                density = 8.0*g/cm3;
                matSteel = new G4Material("Steel", density, ncomponents=9);
                for(G4int i=0; i<9; i++)
                    matSteel->AddElement(manager->FindOrBuildElement(z[i]), per[i]);
                printf("\tSteel has been built.\n");
            }
            theMaterial = matSteel;
            break;
            
        case PLASTIC:
            if(!matPlastic){
                density = 1.415*g/cm3;
                z[0] = 1; z[1] = 6; z[2] = 8;
                matPlastic = new G4Material("Plastic", density, ncomponents=3);
                matPlastic->AddElement(manager->FindOrBuildElement(z[0]), natoms=2);
                matPlastic->AddElement(manager->FindOrBuildElement(z[1]), natoms=1);
                matPlastic->AddElement(manager->FindOrBuildElement(z[2]), natoms=1);
                printf("\tPlastic (polyoxymethylene) has been built.\n");
            }
            theMaterial = matPlastic;
            break;

        case PMMA:
            if(!matPMMA){
                density = 1.19*g/cm3;
                z[0] = 6; z[1] = 8; z[2] = 1;
                matPMMA = new G4Material("PMMA", density, ncomponents=3);
                matPMMA->AddElement(manager->FindOrBuildElement(z[0]), natoms=5);
                matPMMA->AddElement(manager->FindOrBuildElement(z[1]), natoms=2);
                matPMMA->AddElement(manager->FindOrBuildElement(z[2]), natoms=8);
                G4double indexPMMA = 1.4986;
                G4double absPMMA = 1000.0*cm;
                SetUniformOpticalProperties(matPMMA, indexPMMA, absPMMA);
                printf("\t%s has been built.\n",matPMMA->GetName().data());
            }
            theMaterial = matPMMA;
            break;

        case MUMETAL:
            if(!matMuMetal){
                density = 8.25*g/cm3;
                z[0] = 28; z[1] = 26;
                matMuMetal = new G4Material("Mu Metal", density, ncomponents=2);
                matMuMetal->AddElement(manager->FindOrBuildElement(z[0]), 0.8);
                matMuMetal->AddElement(manager->FindOrBuildElement(z[1]), 0.2);
                printf("\tMu Metal has been built.\n");
            }
            theMaterial = matMuMetal;
            break;
            
        case COPPER:
            if(!matCopper){
                density = 8.960*g/cm3;
                z[0] = 29;
                matCopper = new G4Material("Copper", density, ncomponents=1);
                matCopper->AddElement(manager->FindOrBuildElement(z[0]), fractionmass=1.0);
                printf("\tCopper has been built.\n");
            }
            theMaterial = matCopper;
            break;

        case VACUUM:
            if(!matVacuum){
                G4double atomicNumber = 1.;
                G4double massOfMole = 1.008*g/mole;
                density = 1.e-25*g/cm3;
                matVacuum = new G4Material("vacuum", atomicNumber, massOfMole, density, kStateGas, 2.73*kelvin, 3.e-18*pascal);
                G4double indexVacuum = 1.0;
                G4double absVacuum = std::numeric_limits<G4double>::max();
                SetUniformOpticalProperties(matVacuum, indexVacuum, absVacuum);
                printf("\t%s has been built.\n",matVacuum->GetName().data());
            }
            theMaterial = matVacuum;
            break;
            
        case H2O:
            if(!matH2O){
                density = 1.0*g/cm3;
                z[0] = 1; z[1] = 8;
                matH2O = new G4Material("Water", density, ncomponents=2, kStateLiquid);
                matH2O->AddElement(manager->FindOrBuildElement(z[0]), natoms=2);
                matH2O->AddElement(manager->FindOrBuildElement(z[1]), natoms=1);
                G4double indexH2O = 1.33;
                SetOpticalProperties(matH2O, indexH2O, h2oRefl, "opticalData/pureWaterAbsLength_formatted.dat");
                printf("\t%s has been built.\n",matH2O->GetName().data());
            }
            theMaterial = matH2O;
            break;

        case D2O:
            if(!matD2O){
                density = 1.11*g/cm3;
                G4Element* elD = new G4Element("Deuterium", "D", 1);
                elD->AddIsotope(new G4Isotope("H-2", 1, 2, 2.01410178*g/mole), 1.0);
                z[0] = 8;
                matD2O = new G4Material("Heavy water", density, ncomponents=2, kStateLiquid);
                matD2O->AddElement(manager->FindOrBuildElement(z[0]), natoms=1);
                matD2O->AddElement(elD, natoms=2);
                G4double indexD2O = 1.33;
                SetOpticalProperties(matD2O, indexD2O, h2oRefl, "opticalData/pureWaterAbsLength_formatted.dat");
                printf("\t%s has been built.\n",matD2O->GetName().data());
            }
            theMaterial = matD2O;
            break;

        case FUSEDSILICA:
            if(!matFusedSilica){
                density = 2.203*g/cm3;
                z[0] = 14; z[1] = 8;
                matFusedSilica = new G4Material("Fused Silica", density, ncomponents=2);
                matFusedSilica->AddElement(manager->FindOrBuildElement(z[0]), natoms=1);
                matFusedSilica->AddElement(manager->FindOrBuildElement(z[1]), natoms=2);
                G4double indexFusedSilica = 1.4585;
                G4double absFusedSilica = 1000.0*cm;
                SetUniformOpticalProperties(matFusedSilica, indexFusedSilica, absFusedSilica);
                printf("\t%s has been built.\n",matFusedSilica->GetName().data());
            }
            theMaterial = matFusedSilica;
            break;

        case BOROSILICATE:
            if(!matBorosilicate){
                density = 2.23*g/cm3;
                per[0] = 0.040064; z[0] = 5;
                per[1] = 0.539562; z[1] = 8;
                per[2] = 0.028191; z[2] = 11;
                per[3] = 0.011644; z[3] = 13;
                per[4] = 0.377220; z[4] = 14;
                per[5] = 0.003321; z[5] = 19;
                matBorosilicate = new G4Material("Borosilicate glass", density, ncomponents=6);
                for(G4int i=0; i<6; i++)
                    matBorosilicate->AddElement(manager->FindOrBuildElement(z[i]), per[i]);
                G4double indexBorosilicate = 1.517;
                G4double absBorosilicate = 1000.0*cm;
                SetUniformOpticalProperties(matBorosilicate, indexBorosilicate, absBorosilicate);
                printf("\t%s has been built.\n",matBorosilicate->GetName().data());
            }
            theMaterial = matBorosilicate;
            break;
            
        case PHOTOCATHODE:
            if(!matPhotoCathode){
                density = 2.23*g/cm3;
                z[0] = 13;
                matPhotoCathode = new G4Material("Photo cathode (arb.)", density, ncomponents=1);
                matPhotoCathode->AddElement(manager->FindOrBuildElement(z[0]), natoms=1);
                G4double indexPhotoCathode = 1.517;
                G4double absPhotoCathode = 0.0000001*mm;
                SetUniformOpticalProperties(matPhotoCathode, indexPhotoCathode, absPhotoCathode);
                printf("\t%s has been built.\n",matPhotoCathode->GetName().data());
            }
            theMaterial = matPhotoCathode;
            break;

        case TEFLON:
            if(!matTeflon){
                density = 2.2*g/cm3;
                z[0] = 6; z[1] = 9;
                matTeflon = new G4Material("Teflon", density, ncomponents=2);
                matTeflon->AddElement(manager->FindOrBuildElement(z[0]), natoms=2);
                matTeflon->AddElement(manager->FindOrBuildElement(z[1]), natoms=4);
                G4double indexTeflon = 1.4;
                G4double absTeflon = 0.01*cm;
                SetUniformOpticalProperties(matTeflon, indexTeflon, absTeflon);
                printf("\t%s has been built.\n",matTeflon->GetName().data());
            }
            theMaterial = matTeflon;
            break;

        case TYVEK:
            if(!matTyvek){ 
                density = 0.406367041199*g/cm3;
                z[0] = 6; z[1] = 1;
                matTyvek = new G4Material("Tyvek", density, ncomponents=2);
                matTyvek->AddElement(manager->FindOrBuildElement(z[0]), natoms=2);
                matTyvek->AddElement(manager->FindOrBuildElement(z[1]), natoms=4);
                G4double indexTyvek = 1.5;
                G4double absTyvek = 0.01*cm;
                SetUniformOpticalProperties(matTyvek, indexTyvek, absTyvek);
                printf("\t%s has been built.\n",matTyvek->GetName().data());
            }
            theMaterial = matTyvek;
            break;

        case DELRIN:
            theMaterial = matDelrin;          
            break;

        case EPOXY:
            if(!matEpoxy){ 
                density = 1.217;
                z[0] = 6; z[1] = 1; z[2] = 8;
                matEpoxy = new G4Material("Epoxy", density, ncomponents=3);
                matEpoxy->AddElement(manager->FindOrBuildElement(z[0]), natoms=15);
                matEpoxy->AddElement(manager->FindOrBuildElement(z[1]), natoms=16);
                matEpoxy->AddElement(manager->FindOrBuildElement(z[2]), natoms=2);
                G4double indexEpoxy = 1.5;
                G4double absEpoxy = 0.01*cm;
                SetUniformOpticalProperties(matEpoxy, indexEpoxy, absEpoxy);
                printf("\t%s has been built.\n",matEpoxy->GetName().data());
            }
            theMaterial = matEpoxy;
            break;

        case CONCRETE:
            density = 2.30*g/cm3;
            z[0] = 1; per[0] = 0.022100;
            z[1] = 6; per[1] = 0.002484;
            z[2] = 8; per[2] = 0.574930;
            z[3] = 11; per[3] = 0.015208;
            z[4] = 12; per[4] = 0.001266;
            z[5] = 13; per[5] = 0.019953;
            z[6] = 14; per[6] = 0.304627;
            z[7] = 19; per[7] = 0.010045;
            z[8] = 20; per[8] = 0.042951;
            z[9] = 26; per[9] = 0.006435;
            matConcrete = new G4Material("Concrete", density, ncomponents=10);
            for(G4int i=0; i<10; i++)
                matConcrete->AddElement(manager->FindOrBuildElement(z[i]), per[i]);
            theMaterial = matConcrete;
            break;

        default:
            G4cout<<"ERROR: Requested Material does not exist.\n";
            exit(0);
            break;
    }
    
    return theMaterial;
}

void G4d2oMaterialsDefinition::SetOpticalProperties(G4Material *theMat, G4double theIndex, G4double h2oReflF, G4String absLengthFile)
{
    G4int numEn;
    G4double *energy = 0;
    G4double *absLength = 0;
    G4double *refrIndex = 0;
    G4bool bFileGood = true;
    std::ifstream inFile(absLengthFile);
    if(inFile){
        inFile >> numEn;
        if(numEn >= 1){
            G4double *lambda = new G4double[numEn];
            energy = new G4double[numEn];
            absLength = new G4double[numEn];
            refrIndex = new G4double[numEn];
            for(G4int iEn=numEn-1; iEn>=0; iEn--){
                inFile >> lambda[iEn] >> absLength[iEn];
                absLength[iEn] = h2oReflF * absLength[iEn];
                if(lambda[iEn]>0) energy[iEn] = h_Planck*c_light/(lambda[iEn]*nm);
                else bFileGood = false;
                absLength[iEn] *= cm;
                refrIndex[iEn] = theIndex;
                G4cout<<"\t"<<iEn<<" "<<energy[iEn]/eV<<" "<<refrIndex[iEn]<<" "<<absLength[iEn]/m<<G4endl;
            }
        }
        else bFileGood = false;
        inFile.close();
    }
    else bFileGood = false;

    if(!bFileGood){
        printf("Error in G4d2oMaterialsDefinition::SetOpticalProperties() - check %s\n",absLengthFile.data());
        exit(0);
    }

    if(numEn>0 && energy && absLength && refrIndex){
        G4MaterialPropertiesTable *matProp = new G4MaterialPropertiesTable();
        matProp->AddProperty("RINDEX",energy,refrIndex,numEn);
        matProp->AddProperty("ABSLENGTH",energy,absLength,numEn);
        theMat->SetMaterialPropertiesTable(matProp);
    }
    else{
        printf("Error in G4d2oMaterialsDefinition::SetOpticalProperties()\n");
        exit(0);
    }
}

void G4d2oMaterialsDefinition::SetUniformOpticalProperties(G4Material *theMat, G4double theIndex, G4double theAbsLength)
{
    const G4int numEnergies = 3;
    G4double energy[numEnergies] = {1.5*eV, 4.25*eV, 7.0*eV};
    G4double refrIndex[numEnergies] = {theIndex,theIndex,theIndex};
    G4double absLength[numEnergies] = {theAbsLength,theAbsLength,theAbsLength};
    
    G4MaterialPropertiesTable *matProp = new G4MaterialPropertiesTable();
    matProp->AddProperty("RINDEX",energy,refrIndex,numEnergies);
    matProp->AddProperty("ABSLENGTH",energy,absLength,numEnergies);
    theMat->SetMaterialPropertiesTable(matProp);
}

void G4d2oMaterialsDefinition::SetReflector(G4VPhysicalVolume *theExitingVolume, 
                                            G4VPhysicalVolume *theEnteringVolume,
                                            G4double theReflectivity, 
                                            G4double theSigmaAlpha, 
                                            G4SurfaceType type)
{
    G4OpticalSurface *OptSurfReflector = new G4OpticalSurface("Reflector");
    OptSurfReflector->SetType(type);
    OptSurfReflector->SetModel(unified);
    if(theSigmaAlpha <= 0.0) {
        OptSurfReflector->SetFinish(polished);
    } else {
        OptSurfReflector->SetFinish(ground);
        OptSurfReflector->SetSigmaAlpha(theSigmaAlpha);
    }
    
    const G4int numEnergies = 3;
    G4double energy[numEnergies] = {1.5*eV, 4.25*eV, 7.0*eV};
    G4double reflectivity[numEnergies] = {theReflectivity,theReflectivity,theReflectivity};
    G4double transmittance[numEnergies] = {1.0 - theReflectivity, 1.0 - theReflectivity, 1.0 - theReflectivity};
    
    G4MaterialPropertiesTable *MPTReflect = new G4MaterialPropertiesTable();
    MPTReflect->AddProperty("REFLECTIVITY",energy,reflectivity,numEnergies);
    MPTReflect->AddProperty("TRANSMITTANCE",energy,transmittance,numEnergies);
    
    OptSurfReflector->SetMaterialPropertiesTable(MPTReflect);
    
    G4String surfaceName(Form("optSurf_%s_%s",theExitingVolume->GetName().data(),theEnteringVolume->GetName().data()));
    new G4LogicalBorderSurface(surfaceName,theExitingVolume,theEnteringVolume,OptSurfReflector);
}

// ============================================================
// PURE DATA-DRIVEN Tyvek reflector (from Chavarria 2007 thesis)
// ============================================================
void G4d2oMaterialsDefinition::SetDataDrivenReflector(G4VPhysicalVolume *theExitingVolume,
                                                       G4VPhysicalVolume *theEnteringVolume,
                                                       G4double theReflectivity,
                                                       const G4String& dataDirectory)
{
    G4d2oDataDrivenReflector::CreateInstance(theReflectivity, dataDirectory);
    
    G4OpticalSurface* optSurf = new G4OpticalSurface("Tyvek_Reflector");
    optSurf->SetType(dielectric_dielectric);
    optSurf->SetModel(unified);
    optSurf->SetFinish(polishedfrontpainted);
    
    const G4int numEnergies = 3;
    G4double energy[numEnergies] = {1.5*eV, 4.25*eV, 7.0*eV};
    G4double reflectivityArr[numEnergies] = {theReflectivity, theReflectivity, theReflectivity};
    G4double transmittance[numEnergies] = {1.0 - theReflectivity, 1.0 - theReflectivity, 1.0 - theReflectivity};
    
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("REFLECTIVITY", energy, reflectivityArr, numEnergies);
    mpt->AddProperty("TRANSMITTANCE", energy, transmittance, numEnergies);
    
    optSurf->SetMaterialPropertiesTable(mpt);
    
    G4String surfaceName(Form("DataDriven_%s_%s", 
                              theExitingVolume->GetName().data(),
                              theEnteringVolume->GetName().data()));
    new G4LogicalBorderSurface(surfaceName, theExitingVolume, theEnteringVolume, optSurf);
    
    G4cout << "SetDataDrivenReflector: Created PURE DATA-DRIVEN reflector for " 
           << theExitingVolume->GetName() << " -> " << theEnteringVolume->GetName() << G4endl;
}

void G4d2oMaterialsDefinition::AttachReflectionTree(TFile* file)
{
    G4d2oDataDrivenReflector::SetRootFile(file);
}

void G4d2oMaterialsDefinition::CloseReflectionTree()
{
    G4d2oDataDrivenReflector::CloseRootTree();
}

void G4d2oMaterialsDefinition::SetCurrentEventNumber(G4int eventNum)
{
    G4d2oDataDrivenReflector::SetCurrentEventNumber(eventNum);
}

