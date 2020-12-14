/*
 * ProtonDetectorData.hh
 *
 *  Created on: Dec 19, 2013
 *      Author: perezlou
 */

#ifndef PROTONDETECTORDATA_HH_
#define PROTONDETECTORDATA_HH_

#define NPADS 5
#define NVETO 8

#include "TROOT.h"

class ProtonDetectorData {

private:

	Double_t energyOnGas_beam;
	Double_t energyOnGas_protons;
	Double_t energyOnGas_betas;
	Double_t energyOnGas_alphas;
	Double_t energyOnGas_total;
	Double_t energyOnPad_betas[NPADS];
	Double_t energyOnPad_protons[NPADS];
	Double_t energyOnPad_total[NPADS];
	Double_t energyOnveto_betas[NVETO];
	Double_t energyOnveto_protons[NVETO];
	Double_t energyOnveto_total[NVETO];
	Double_t energyOnDegrader;
	Double_t energyOnVeto;
	Double_t energyOnVeto_protons;
	Double_t energyOnVeto_betas;
	Double_t stepSumLengthOnGas_beam;
	Double_t stepSumLengthOnGas_protons;
	Double_t stepSumLengthOnGas_betas;
	Double_t stepSumLengthOnGas_alphas;
	Double_t beam_last_positionX;
	Double_t beam_last_positionY;
	Double_t beam_last_positionZ;
	Double_t proton_last_positionX;
	Double_t proton_last_positionY;
	Double_t proton_last_positionZ;
	Double_t alpha_last_positionX;
	Double_t alpha_last_positionY;
	Double_t alpha_last_positionZ;
        Int_t proton_multiplicity;
        Int_t beta_multiplicity;
        Int_t total_multiplicity;
        Int_t proton_Vetomultiplicity;
        Int_t beta_Vetomultiplicity;
        Int_t total_Vetomultiplicity;
	Int_t eventID;
	Int_t runID;

        //New initial Energies

        Double_t BetaKineticEnergy;
        Double_t ProtonKineticEnergy;
        Double_t RecoilKineticEnergy;

public:
	ProtonDetectorData();
	~ProtonDetectorData();

	Double_t GetEnergyOnGas_beam(){return energyOnGas_beam;}
	Double_t GetEnergyOnGas_protons(){return energyOnGas_protons;}
	Double_t GetEnergyOnGas_betas(){return energyOnGas_betas;}
	Double_t GetEnergyOnGas_alphas(){return energyOnGas_alphas;}
	Double_t GetEnergyOnGas_total(){return energyOnGas_total;}
	Double_t GetEnergyOnDegrader(){return energyOnDegrader;}
	Double_t GetEnergyOnVeto(){return energyOnVeto;}
	Double_t GetEnergyOnVeto_protons(){return energyOnVeto_protons;}
	Double_t GetEnergyOnVeto_betas(){return energyOnVeto_betas;}
	Double_t GetStepSumLengthOnGas_beam(){return stepSumLengthOnGas_beam;}
	Double_t GetStepSumLengthOnGas_protons(){return stepSumLengthOnGas_protons;}
	Double_t GetStepSumLengthOnGas_betas(){return stepSumLengthOnGas_betas;}
	Double_t GetStepSumLengthOnGas_alphas(){return stepSumLengthOnGas_alphas;}
	Double_t GetBeamLastPositionX(){return beam_last_positionX;}
	Double_t GetBeamLastPositionY(){return beam_last_positionY;}
	Double_t GetBeamLastPositionZ(){return beam_last_positionZ;}
	Double_t GetProtonLastPositionX(){return proton_last_positionX;}
	Double_t GetProtonLastPositionY(){return proton_last_positionY;}
	Double_t GetProtonLastPositionZ(){return proton_last_positionZ;}
	Double_t GetAlphaLastPositionX(){return alpha_last_positionX;}
	Double_t GetAlphaLastPositionY(){return alpha_last_positionY;}
	Double_t GetAlphaLastPositionZ(){return alpha_last_positionZ;}
        Int_t GetProtonMultiplicity(){return proton_multiplicity;}
        Int_t GetBetaMultiplicity(){return beta_multiplicity;}
        Int_t GetTotalMultiplicity(){return total_multiplicity;}
        Int_t GetProtonVetoMultiplicity(){return proton_Vetomultiplicity;}
        Int_t GetBetaVetoMultiplicity(){return beta_Vetomultiplicity;}
        Int_t GetTotalVetoMultiplicity(){return total_Vetomultiplicity;}
	Int_t GetEventID(){return eventID;}
	Int_t GetRunID(){return runID;}
	Double_t GetEnergyOnPad_protons(Int_t i){return energyOnPad_protons[i];}
	Double_t GetEnergyOnPad_betas(Int_t i){return energyOnPad_betas[i];}
        Double_t GetEnergyOnPad_total(Int_t i){return energyOnPad_total[i];}
	Double_t GetEnergyOnVeto_protons(Int_t i){return energyOnveto_protons[i];}
	Double_t GetEnergyOnVeto_betas(Int_t i){return energyOnveto_betas[i];}
        Double_t GetEnergyOnVeto_total(Int_t i){return energyOnveto_total[i];}


       Double_t GetBetaKineticEnergy(){return BetaKineticEnergy;}
       Double_t GetProtonKineticEnergy(){return ProtonKineticEnergy;}
       Double_t GetRecoilKineticEnergy(){return RecoilKineticEnergy;}



        void SetEnergyOnGas_beam(Double_t val){energyOnGas_beam=val;}
	void SetEnergyOnGas_protons(Double_t val){energyOnGas_protons=val;}
	void SetEnergyOnGas_betas(Double_t val){energyOnGas_betas=val;}
	void SetEnergyOnGas_alphas(Double_t val){energyOnGas_alphas=val;}
	void SetEnergyOnGas_total(Double_t val){energyOnGas_total=val;}
	void SetEnergyOnDegrader(Double_t val){energyOnDegrader=val;}
	void SetEnergyOnVeto(Double_t val){ energyOnVeto=val;}
	void SetEnergyOnVeto_protons(Double_t val){ energyOnVeto_protons=val;}
	void SetEnergyOnVeto_betas(Double_t val){ energyOnVeto_betas=val;}
	void SetStepSumLengthOnGas_beam(Double_t val){  stepSumLengthOnGas_beam=val;}
	void SetStepSumLengthOnGas_protons(Double_t val){  stepSumLengthOnGas_protons=val;}
	void SetStepSumLengthOnGas_betas(Double_t val){  stepSumLengthOnGas_betas=val;}
	void SetStepSumLengthOnGas_alphas(Double_t val){  stepSumLengthOnGas_alphas=val;}
	void SetBeamLastPositionX(Double_t val){  beam_last_positionX=val;}
	void SetBeamLastPositionY(Double_t val){  beam_last_positionY=val;}
	void SetBeamLastPositionZ(Double_t val){  beam_last_positionZ=val;}
	void SetProtonLastPositionX(Double_t val){  proton_last_positionX=val;}
	void SetProtonLastPositionY(Double_t val){  proton_last_positionY=val;}
	void SetProtonLastPositionZ(Double_t val){  proton_last_positionZ=val;}
	void SetAlphaLastPositionX(Double_t val){  alpha_last_positionX=val;}
	void SetAlphaLastPositionY(Double_t val){  alpha_last_positionY=val;}
	void SetAlphaLastPositionZ(Double_t val){  alpha_last_positionZ=val;}
        void SetProtonMultiplicity(Int_t val){proton_multiplicity=val;}
        void SetBetaMultiplicity(Int_t val){beta_multiplicity=val;}
        void SetTotalMultiplicity(Int_t val){total_multiplicity=val;}
        void SetProtonVetoMultiplicity(Int_t val){proton_Vetomultiplicity=val;}
        void SetBetaVetoMultiplicity(Int_t val){beta_Vetomultiplicity=val;}
        void SetTotalVetoMultiplicity(Int_t val){total_Vetomultiplicity=val;}
	void SetEventID(Int_t val){  eventID =val;}
	void SetRunID(Int_t val){  runID =val;}

        void SetEnergyOnPad_protons(Int_t i, Double_t val){energyOnPad_protons[i]=val;}
        void SetEnergyOnPad_betas(Int_t i,  Double_t val){energyOnPad_betas[i]=val;}
       void SetEnergyOnPad_total(Int_t i,  Double_t val){energyOnPad_total[i]=val;}
       void SetEnergyOnVeto_protons(Int_t i, Double_t val){energyOnveto_protons[i]=val;}
        void SetEnergyOnVeto_betas(Int_t i,  Double_t val){energyOnveto_betas[i]=val;}
       void SetEnergyOnVeto_total(Int_t i,  Double_t val){energyOnveto_total[i]=val;}



       void SetBetaKineticEnergy(Double_t val){ BetaKineticEnergy=val;}
       void SetProtonKineticEnergy(Double_t val){ ProtonKineticEnergy=val;}
       void SetRecoilKineticEnergy(Double_t val){ RecoilKineticEnergy=val;}
 

ClassDef(ProtonDetectorData,1);
};

#endif /* PROTONDETECTORDATA_HH_ */
