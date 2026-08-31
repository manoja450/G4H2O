/************************************************************
 * G4d2oPhysicsList.hh
 ********************************************************************/


#ifndef G4d2oPhysicsList_h
#define G4d2oPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

class G4d2oPhysicsList: public G4VModularPhysicsList
{
public:
    G4d2oPhysicsList(G4bool bOn=true, G4bool bNHP=false);
   	virtual ~G4d2oPhysicsList();
    
protected:
    
    //Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();
    
    //These methods construct physics processes and register them
    void ConstructEM();
    void ConstructRadioactiveDecay();
	void ConstructHadronic();
    void ConstructOptical();

	//Sets the Cuts
	void SetCuts();
    
private:
    
    G4bool bNull, bNeutronHP;
    G4int verbosityLevel;
    
};

#endif

//End of file G4d2oPhysicsList.hh

