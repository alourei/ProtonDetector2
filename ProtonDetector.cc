//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: AstroBox2.cc,v 1.6 2011/08/16 17:47:10 roeder Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
// Modified by Brian Roeder, TAMU on 08/16/2011
// email - broeder@comp.tamu.edu
// 
// -------------------------------------------------------------------
// Based on     GEANT 4 - exampleN01, adding elements from other exps.
// -------------------------------------------------------------------
//
// 8/16/2011 - Simulation of factors for the AstroBox2 concept!
// 6/01/2015 - Merging implant and decay
//
//

#include <ctime>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include <G4UIExecutive.hh>
#include <G4UIterminal.hh>
#include "globals.hh"

#include "MargotDetectorConstruction.hh"
#include "ProtonDetectorConstruction.hh"
#include "ProtonDetectorRunAction.hh"
#include "ProtonDetectorEventAction.hh"
#include "ProtonDetectorSteppingAction.hh"
#include "ProtonDetectorPhysicsList.hh"
#include "MargotPhysicsList.hh"
#include "MargotPrimaryGeneratorAction.hh"
#include "ProtonDetectorPrimaryGeneratorAction.hh"
#include "MargotEventAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

//#include "MargotDataRecordTree.hh"
#include <ProtonDetectorROOTAnalysis.hh>

int main(int argc,char** argv)
{
  // Simulation Setup Variables

  // Choose "Fragment" - Possible Beams are Al23
  G4String FragName = "Ar32";	// Fragment Name  
  //G4String FragName = "Si25";   // Fragment Name
  //G4String FragName = "Al23";   // Fragment Name
  //G4String FragName = "Cl31";   // Fragment Name
  //G4String FragName = "Mg20";   // Fragment Name
  //   G4String FragName = "proton";   // Fragment Name
  //G4String FragName = "geantino";   // Fragment Name
  //G4int FragMass = 20;          // Fragment Mass Number
  //G4int FragMass = 31;          // Fragment Mass Number
  //G4int FragMass = 23;          // Fragment Mass Number
  //G4int FragMass = 25;          // Fragment Mass Number
  G4int FragMass = 32;          // Fragment Mass Number
  //G4int FragMass = 1;          // Fragment Mass Number

  G4double BeamEngPerA = 0.0*MeV;    // Stopped Beam with distrib.
  //G4double BeamEngPerA = 72.312 *MeV;  // Beam Energy in A MeV
  //G4double BeamEngPerA = 6.52174 *MeV;  // Beam Energy in A MeV

  G4String SourceType;
  //SourceType = "test";     // Test beam at given beam energy
  //SourceType = "diffuse";  // New diffuse beam
  //SourceType = "pencil";   // pencil beam with +-0.5% energy
  //SourceType = "RealBeam";   // the above w/ 5mm beam spot + 1 deg. div.
  //SourceType = "FragmentBeam";   // the above w/ 5mm beam spot + 1 deg. div.
   //SourceType = "iso";      // isotropic beam
  //SourceType = "conic";    // conic beam with dims. of DemonDet.
  SourceType = "General";    // General
  //SourceType = "GPS";    // GPS source

  G4String DetType;
  DetType = "Al23_GasIC_NoGap";    // No Gap between MicroBulks
  //DetType = "Al23_GasIC_wGap"; // 100um Gap between MicroBulks
  //DetType = "Al23_GasIC_wGap1mm"; // 1mm Gap between MicroBulks

 
  G4int numSourcesPerEvt = 1;

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
  //
  G4VUserDetectorConstruction* detector = new ProtonDetectorConstruction();
  runManager->SetUserInitialization(detector);
  //
  //G4VUserPhysicsList* physics = new MargotPhysicsList;
  G4VUserPhysicsList* physics = new ProtonDetectorPhysicsList;
   runManager->SetUserInitialization(physics);


 G4UIExecutive *session=0;
 G4UImanager* UI=G4UImanager::GetUIpointer();
 G4VisManager* visManager=0; 

 if (argc==1){   // Define UI session for interactive mode.

  // Visualization, if you choose to have it!
#ifdef G4UI_USE
   session =new G4UIExecutive(argc,argv,"tcsh");
#endif

#ifdef G4VIS_USE
   session =new G4UIExecutive(argc,argv,"qt");
   visManager = new G4VisExecutive;
   visManager->Initialize();
   G4cout<<UI<<G4endl;
#endif

 }
 else if(argc==2 && atoi(argv[1])==1){
   session =new G4UIExecutive(argc,argv,"tcsh");
 
 }
 else {;
   //This part should be used for batch
 }
 

  ProtonDetectorROOTAnalysis *theAnalysis = new ProtonDetectorROOTAnalysis();

  //Histograming

  // Initialize Primary Generator Action
  G4double Beam_Energy = static_cast<G4double>(FragMass)*BeamEngPerA;
  //G4VUserPrimaryGeneratorAction* gen_action = new MargotPrimaryGeneratorAction(numSourcesPerEvt,FragName,SourceType,Beam_Energy);
  //G4VUserPrimaryGeneratorAction* gen_action = new ProtonDetectorPrimaryGeneratorAction(numSourcesPerEvt,FragName,SourceType,Beam_Energy);
  ProtonDetectorPrimaryGeneratorAction* gen_action = new ProtonDetectorPrimaryGeneratorAction(numSourcesPerEvt,FragName,SourceType,Beam_Energy);
  //gen_action->SetBeamData("BeamProperties.dat");
  gen_action->SetBeamData("BeamProperties_LISE.dat");
  runManager->SetUserAction(gen_action);

  //set optional user action classes

  //Added Run actions
   ProtonDetectorRunAction *run_action=new ProtonDetectorRunAction;
   runManager->SetUserAction(run_action);

  //Include Event Action Sequence to invoke event by event analysis
  //G4UserEventAction* event_action = new MargotEventAction();
  ProtonDetectorEventAction* event_action = new ProtonDetectorEventAction();
  runManager->SetUserAction(event_action);
  
  ProtonDetectorSteppingAction *step_action = new ProtonDetectorSteppingAction();
  runManager->SetUserAction(step_action);


  // Initialize G4 kernel
  //
  runManager->Initialize();

  // Get the pointer to the UI manager and set verbosities

 
  UI->ApplyCommand("/tracking/storeTrajectory 1");

#ifdef G4UI_USE
  if (argc==1){
#ifdef G4VIS_USE
   UI->ApplyCommand("/vis/scene/create") ;
   UI->ApplyCommand("/vis/open OGL") ;
   UI->ApplyCommand("/vis/drawVolume");
   UI->ApplyCommand("/tracking/storeTrajectory 1");
   UI->ApplyCommand("/vis/scene/endOfEventAction accumulate");
   UI->ApplyCommand("/vis/scene/add/trajectories");
#endif
  session->SessionStart();
  delete session;
 }

  else if(argc==2 && atoi(argv[1])==1){
  session->SessionStart();
  delete session;
  }
  else{
   G4String command = "/control/execute ";
   G4String fileName = argv[1];
   UI->ApplyCommand(command+fileName);
}
#endif
  // Start a run  
  //runManager->BeamOn(numberOfEvent);
  // Job termination
  //
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
  //
  G4cout
  << "############################################################" << G4endl;
  G4cout 
    << " Run Summary" << G4endl;
  G4cout
    << "############################################################" << G4endl;
  //  G4cout << "The Catcher Angle was set to : " << DetAng/deg << " deg." << G4endl; 

  delete theAnalysis;
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete theEngine;
  delete runManager;

  G4cout << "End of Simulation! " << G4endl;

  G4int end_time = time(NULL);
  G4double run_time = static_cast<double>(end_time) - static_cast<double>(start_time);
  if(run_time < 60.)
    {G4cout << "The total run time was " << run_time << " seconds." << G4endl;}
  else
    {
      run_time /= 60.;
      G4cout << "The total run time was " << run_time << " minutes." << G4endl;
    }

  return 0;
}


