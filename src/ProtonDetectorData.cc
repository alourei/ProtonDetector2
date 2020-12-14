/*
 * ProtonDetectorData.cc
 *
 *  Created on: Dec 19, 2013
 *      Author: perezlou
 */

#include <ProtonDetectorData.hh>

ClassImp(ProtonDetectorData);

ProtonDetectorData::ProtonDetectorData() {
	// TODO Auto-generated constructor stub

	 energyOnGas_beam =0;
	 energyOnGas_protons =0;
	 energyOnGas_betas =0;
	 energyOnGas_alphas =0;
	 energyOnGas_total =0;
	 energyOnDegrader =0;
	 energyOnVeto=0;
	 energyOnVeto_protons=0;
	 energyOnVeto_betas=0;
	 stepSumLengthOnGas_beam =0;
	 stepSumLengthOnGas_protons =0;
	 stepSumLengthOnGas_betas =0;
	 stepSumLengthOnGas_alphas =0;
	 beam_last_positionX =-999;
	 beam_last_positionY =-999;
	 beam_last_positionZ =-999;
	 proton_last_positionX =-999;
	 proton_last_positionY =-999;
	 proton_last_positionZ =-999;
	 alpha_last_positionX =-999;
	 alpha_last_positionY =-999;
	 alpha_last_positionZ =-999;
         proton_multiplicity=0;
         beta_multiplicity=0;
         total_multiplicity=0;
         proton_Vetomultiplicity=0;
         beta_Vetomultiplicity=0;
         total_Vetomultiplicity=0;
	 eventID =0;
	 runID =0;

	 BetaKineticEnergy=0;
	 ProtonKineticEnergy=0;
	 RecoilKineticEnergy=0;


	 for(Int_t j=0;j<NPADS;j++){
	   energyOnPad_betas[j]=0.;
	   energyOnPad_protons[j]=0.;
	   energyOnPad_total[j]=0.;
	   
	 }
	 for(Int_t j=0;j<NVETO;j++){
	   energyOnveto_betas[j]=0.;
	   energyOnveto_protons[j]=0.;
	   energyOnveto_total[j]=0.;	   
	 }

}

ProtonDetectorData::~ProtonDetectorData() {
	// TODO Auto-generated destructor stub
}

