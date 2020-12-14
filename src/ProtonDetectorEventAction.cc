/*
 * ProtonDetectorEventAction.cc
 *
 *  Created on: Dec 11, 2013
 *      Author: perezlou
 */

#include <ProtonDetectorEventAction.hh>


#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "G4ThreeVector.hh"

#include "G4RunManager.hh"

#include "Randomize.hh"
#include <iomanip>

#include <ProtonDetectorROOTAnalysis.hh>


ProtonDetectorEventAction::ProtonDetectorEventAction():drawFlag("all"),printModulo(100) {
	// TODO Auto-generated constructor stub


}

ProtonDetectorEventAction::~ProtonDetectorEventAction() {
	// TODO Auto-generated destructor stub
}

void ProtonDetectorEventAction::BeginOfEventAction(const G4Event* evt){

	  G4int evtNb = evt->GetEventID();

	  const G4int verboseLevel = G4RunManager::GetRunManager()->GetVerboseLevel();
	    if(verboseLevel>0){
	      if (evtNb%printModulo == 0) {
	        G4cout << "##################################################################"
	  	     << G4endl
	  	     << "########    ProtonDetectorEventAction::BeginOfEventAction()   ##########"
	  	     << G4endl
	  	     << "########           Begin of event: " << evtNb << "        ########"<<  G4endl;
	        CLHEP::HepRandom::showEngineStatus();
	        G4cout << "##################################################################"
	  	     << G4endl;
	      }
	    }

	    if (gProtonDetectorROOTAnalysis) gProtonDetectorROOTAnalysis->BeginOfEventAction(evt);


}

void ProtonDetectorEventAction::EndOfEventAction(const G4Event *evt){

	  //
	  //  After the end of the event...
	  //
	  G4int evtNb = evt->GetEventID();

	  if (evtNb%printModulo == 0){
	    G4cout << "##################################################################"
		   << G4endl
		   << "#########    ProtonDetectorEventAction::EndOfEventAction()   #########"
		   << G4endl
		   << "#### End of event: " << evtNb << G4endl;
	    G4cout << "##################################################################"
		   << G4endl;
	  }

	  // Histogramming
	  if(gProtonDetectorROOTAnalysis)
	    gProtonDetectorROOTAnalysis->EndOfEventAction(evt);

}
