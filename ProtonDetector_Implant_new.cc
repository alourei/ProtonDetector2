/*
 * ProtonDetector_Implant.cc
 *
 *  Created on: Dec 18, 2013
 *      Author: perezlou
 */


#include <ctime>

#include <G4RunManager.hh>
#include "G4UImanager.hh"
#include <G4UIExecutive.hh>
#include "globals.hh"

#include "ProtonDetectorConstruction.hh"
#include "ProtonDetectorRunAction.hh"
#include "ProtonDetectorEventAction.hh"
#include "ProtonDetectorSteppingAction.hh"
#include "MargotPhysicsList.hh"
#include "ProtonDetectorPrimaryGeneratorAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include <ProtonDetectorROOTAnalysis.hh>

int main(int argc,char** argv)
{
	 // Simulation Setup Variables
	  G4String FragName = "Al23";   // Fragment Name
	  //G4String FragName = "proton";   // Fragment Name
	  //G4String FragName = "geantino";   // Fragment Name
	  G4int FragMass = 23;          // Fragment Mass Number

	  //G4double BeamEngPerA = 0.0*MeV;    // Stopped Beam with distrib.
	  G4double BeamEngPerA = 40.0126*MeV;  // Beam Energy in A MeV

	  G4String SourceType;
	  //SourceType = "test";     // Test beam at given beam energy
	  //SourceType = "diffuse";  // New diffuse beam
	  //SourceType = "pencil";   // pencil beam with +-0.5% energy
	  SourceType = "RealBeam";   // the above w/ 5mm beam spot + 1 deg. div.
	  //SourceType = "iso";      // isotropic beam
	  //SourceType = "conic";    // conic beam with dims. of DemonDet.
	  //SourceType = "General";    // GPS source

	  G4int numSourcesPerEvt = 1;
	  G4int numberOfEvent=1;

	  G4int VerboseFlag = 0;
	  G4int VisFlag = 1;

	  //*******************************************************************
	  // Construct the default run manager
	  //
	  G4RunManager* runManager = new G4RunManager;

	  // Start Time and Random Number Engine
	  G4int start_time = time(NULL);  // stores initial time of program start
	  CLHEP::RanecuEngine* theEngine = new CLHEP::RanecuEngine();
	  CLHEP::HepRandom::setTheEngine(theEngine);
	  CLHEP::HepRandom::setTheSeed(start_time);
	  CLHEP::HepRandom::showEngineStatus();

	  // set mandatory initialization classes

	  ProtonDetectorConstruction *detector = new ProtonDetectorConstruction();
	  runManager->SetUserInitialization(detector);

	  MargotPhysicsList* physics = new MargotPhysicsList();
	  runManager->SetUserInitialization(physics);

	#ifdef G4UI_USE

	  G4UIExecutive* session =new G4UIExecutive(argc,argv,"tcsh");

 	#endif

	#ifdef G4VIS_USE
	  if(VisFlag ==1)
		  session =new G4UIExecutive(argc,argv,"qt");
	  G4VisManager* visManager = new G4VisExecutive;
	  visManager->Initialize();
	#endif

	 //Histograming
	 ProtonDetectorROOTAnalysis *theAnalysis = new ProtonDetectorROOTAnalysis();

	 // Initialize Primary Generator Action
	  G4double Beam_Energy = static_cast<G4double>(FragMass)*BeamEngPerA;
	  G4VUserPrimaryGeneratorAction* gen_action = new ProtonDetectorPrimaryGeneratorAction(numSourcesPerEvt,FragName,SourceType,Beam_Energy);
	  runManager->SetUserAction(gen_action);

	  //set optional user action classes

	  ProtonDetectorRunAction *run_action=new ProtonDetectorRunAction;
	  runManager->SetUserAction(run_action);

	  ProtonDetectorEventAction* event_action = new ProtonDetectorEventAction();
	  runManager->SetUserAction(event_action);

	  ProtonDetectorSteppingAction *step_action = new ProtonDetectorSteppingAction();
	  runManager->SetUserAction(step_action);




}


