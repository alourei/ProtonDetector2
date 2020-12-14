/*
 * ProtonDetectorSteppingAction.cc
 *
 *  Created on: Dec 13, 2013
 *      Author: perezlou
 */

#include <ProtonDetectorSteppingAction.hh>

#include <ProtonDetectorROOTAnalysis.hh>

ProtonDetectorSteppingAction::ProtonDetectorSteppingAction() {
	// TODO Auto-generated constructor stub

}

ProtonDetectorSteppingAction::~ProtonDetectorSteppingAction() {
	// TODO Auto-generated destructor stub
}

void ProtonDetectorSteppingAction::UserSteppingAction(const G4Step* aStep){


	G4double Edep = aStep->GetTotalEnergyDeposit();

	G4double KineticEnergy = aStep->GetTrack()->GetKineticEnergy();

	G4String volume =aStep->GetTrack()->GetVolume()->GetName();

	G4String material = aStep->GetTrack()->GetMaterial()->GetChemicalFormula();


	//if(volume == "AlDegrader"){
	//	G4double Z = aStep->GetTrack()->GetMaterial()->GetZ();

	//	G4double A = aStep->GetTrack()->GetMaterial()->GetA();

	//	G4cout<<material<<" "<<Z<<" "<<A<<" Energy Deposit "<<Edep/MeV<<" "<<KineticEnergy/MeV<<G4endl;

	//}

	if (gProtonDetectorROOTAnalysis)
	    gProtonDetectorROOTAnalysis->UserSteppingAction(aStep); // original


}
