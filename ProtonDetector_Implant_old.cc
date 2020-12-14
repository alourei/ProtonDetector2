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
// $Id: BetaP_GasIC.cc,v 1.6 2009/11/29 17:47:10 roeder Exp $
// GEANT4 tag $Name: geant4-09-02.p02 $
//
// Modified by Brian Roeder, TAMU on 11/29/2009
// email - broeder@comp.tamu.edu
// 
// -------------------------------------------------------------------
// Based on     GEANT 4 - exampleN01, adding elements from other exps.
// -------------------------------------------------------------------
//
// 11/29/09 - Implantation of 23Al into a Gas Ion Chamber.
//
//


#include <ctime>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIQt.hh"

#include "globals.hh"

#include "MargotDetectorConstruction.hh"
#include "ProtonDetectorConstruction.hh"
#include "MargotPhysicsList.hh"
#include "MargotPrimaryGeneratorAction.hh"
#include "MargotEventAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "MargotDataRecordTree.hh"

int main(int argc,char** argv)
{
  // Simulation Setup Variables

  // Choose "Fragment" - Possible Beams are Al23,Cl31,P27
  G4String FragName = "Al23";   // Fragment Name
  G4int FragMass = 23;          // Fragment Mass Number
  //FragName = "proton";
  //FragMass = 1;

  G4double BeamEngPerA = 40.0126*MeV;  // Beam Energy in A MeV
  //G4double BeamEngPerA = 2.5*MeV;  // Beam Energy in A MeV
 
  G4String SourceType;
  //SourceType = "test";     // Test beam at given beam energy
  //SourceType = "diffuse";  // New diffuse beam
  //SourceType = "pencil";   // pencil beam with +-0.5% energy
  SourceType = "RealBeam";   // the above w/ 5mm beam spot +1 deg. div. 
  //SourceType = "iso";      // isotropic beam
  //SourceType = "conic";    // conic beam with dims. of DemonDet.

  G4String DetType;
  DetType = "Al23_GasIC";    // Gas ionization chamber

  // Degraders for "Setup3" 2/14/2012
  //G4double DegThickness = 0.9985*mm; // ~good for 1/2 atm 23Al, setup 5
  G4double DegThickness = 0.961*mm; // ~good for 800 torr 23Al, setup 5
  //G4double DegThickness = 0.824*mm; // ~good for 3 atm 23Al, setup 5
  //G4double DegThickness = 0.991*mm; // ~good for 600 torr
  //G4double DegThickness = 0.978*mm; // ~good for 800 torr 23Al, setup 3
  //G4double DegThickness = 0.87*mm;   // ~good for 800 torr 27P
  //G4double DegThickness = 0.788*mm;   // ~good for 800 torr 31Cl
  //G4double DegThickness = 0.963*mm; // ~good for 1000 torr
 // Degraders for "Setup4" 3/01/2012
  //G4double DegThickness = 0.990*mm;
  G4double DegRotation = 0.0*deg;
 
  G4int numberOfEvent;
  //numberOfEvent = 2e4;
  numberOfEvent = 10; 

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
  //
  //G4VUserDetectorConstruction* detector = new MargotDetectorConstruction(DetType,DegThickness,DegRotation);
  G4VUserDetectorConstruction* detector = new ProtonDetectorConstruction();
  runManager->SetUserInitialization(detector);
  //
  G4VUserPhysicsList* physics = new MargotPhysicsList;
  runManager->SetUserInitialization(physics);

  // Visualization, if you choose to have it!

  #ifdef G4UI_USE
  G4UIsession* session =new G4UIQt(argc,argv);
 #endif

  
  #ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // Generate Analysis Pointer Class
  MargotDataRecordTree* MargotPointer = new MargotDataRecordTree;
  cout << MargotPointer << endl;

  // Initialize Primary Generator Action
  G4double Beam_Energy = static_cast<G4double>(FragMass)*BeamEngPerA;
  G4VUserPrimaryGeneratorAction* gen_action = new MargotPrimaryGeneratorAction(FragName,SourceType,Beam_Energy);  
  //G4VUserPrimaryGeneratorAction* gen_action = new ProtonDetectorPrimaryGeneratorAction(numSourcesPerEvt,FragName,SourceType,Beam_Energy);

  runManager->SetUserAction(gen_action);

  //set optional user action classes

  //Include Event Action Sequence to invoke event by event analysis
  G4UserEventAction* event_action = new MargotEventAction();
  runManager->SetUserAction(event_action);
  
  // Initialize G4 kernel
  //
  runManager->Initialize();

  // Get the pointer to the UI manager and set verbosities

  G4UImanager* UI = G4UImanager::GetUIpointer();
 
  if (VerboseFlag == 1)
    {
      UI->ApplyCommand("/run/verbose 0");
      UI->ApplyCommand("/event/verbose 1");  
      UI->ApplyCommand("/tracking/verbose 2");
    }
  else
    {
     // silent running 
    }

  // Open visualization window to see stuff made; particles.
  
  if (VisFlag == 1)
    {
      UI->ApplyCommand("/vis/scene/create");
      UI->ApplyCommand("/vis/open OGL");
      //UI->ApplyCommand("/vis/open VRML2FILE");
      UI->ApplyCommand("/tracking/storeTrajectory 1");
      UI->ApplyCommand("/vis/scene/endOfEventAction accumulate");
      UI->ApplyCommand("/vis/scene/add/trajectories");
    }

  // Start a run  
  session->SessionStart();
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
  MargotPointer->GetParticleTotals();

  delete MargotPointer;  // data files managed by MargotDataRecordTree.hh class
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


