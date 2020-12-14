/*
 * ProtonDetectorRunAction.cc
 *
 *  Created on: Dec 11, 2013
 *      Author: perezlou
 */

#include <ProtonDetectorRunAction.hh>
#include "ProtonDetectorROOTAnalysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"


ProtonDetectorRunAction::ProtonDetectorRunAction() {
	// TODO Auto-generated constructor stub

}

ProtonDetectorRunAction::~ProtonDetectorRunAction() {
	// TODO Auto-generated destructor stub
}

void ProtonDetectorRunAction::BeginOfRunAction(const G4Run *aRun){

	  //
	  // Actions to perform at the beginning og the run
	  //

	  const G4int verboseLevel = G4RunManager::GetRunManager()->GetVerboseLevel();
	  if(verboseLevel>2){
	    G4cout << "##################################################################"
		  << G4endl
		  << "###########   ProtonDetectorRunAction::BeginOfRunAction()  ##############"
		  << G4endl
		  << "###    Run " << aRun->GetRunID() << " start." << G4endl;
	    G4cout << "##################################################################"
		 << G4endl;
	  }

	  //inform the runManager to save random number seed
	  G4RunManager::GetRunManager()->SetRandomNumberStore(true);

	  // Histogramming
	  if (gProtonDetectorROOTAnalysis) gProtonDetectorROOTAnalysis->BeginOfRunAction(aRun);

}

void ProtonDetectorRunAction::EndOfRunAction(const G4Run *aRun){

	  //
	  // Actions to perform at the end of the run
	  //
	  if (gProtonDetectorROOTAnalysis) gProtonDetectorROOTAnalysis->EndOfRunAction(aRun);
}

