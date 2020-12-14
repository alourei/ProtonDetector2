/*
 * ProtonDetectorPrimaryGeneratorAction.hh
 *
 *  Created on: Dec 4, 2013
 *      Author: perezlou
 */

#ifndef PROTONDETECTORPRIMARYGENERATORACTION_HH_
#define PROTONDETECTORPRIMARYGENERATORACTION_HH_

#include <globals.hh>
#include <G4LorentzVector.hh>
#include "MargotDataRecordTree.hh"
#include <G4VUserPrimaryGeneratorAction.hh>

#include <TFile.h>
#include <TH1F.h>



class G4ParticleGun;
class G4Event;
class G4GeneralParticleSource;
class G4SingleParticleSource;
class G4SPSEneDistribution;
class G4SPSPosDistribution;
class G4SPSAngDistribution;

class TFile;
class TH1F;



class ProtonDetectorPrimaryGeneratorAction: public G4VUserPrimaryGeneratorAction {
public:
	ProtonDetectorPrimaryGeneratorAction();
	ProtonDetectorPrimaryGeneratorAction(G4int numSources,G4String pName, G4String sourceType, G4double beamEnergy );
	~ProtonDetectorPrimaryGeneratorAction();
	void GeneratePrimaries(G4Event *anEvent);
	void SetBeamData(G4String theFileName);

private:

	G4ParticleGun *particleGun;
	G4GeneralParticleSource *theArchitect;
	G4SingleParticleSource *theSource[2];
	G4int dataEvent;
	G4double particleEnergy;

	//For fragment beams
	G4int *Z_beam;
	G4int *A_beam;
	G4double *Energy_beam;
	G4double *FWHM_Energy_beam;
	G4double *x_mean_beam;
	G4double *FWHM_x_beam;
	G4double *y_mean_beam;
	G4double *FWHM_y_beam;

	G4int numberOfFragments;
	G4double normalization;

	G4double *intensity_beam;
	G4double *probability_beam;
	G4double *probability_limits;

	MargotDataRecordTree* MargotDataOutPG;

	G4int numberOfSources;
	G4String sourceType;
	G4String beamName;

	G4SPSEneDistribution *GPS_ParticleEnergy;
	G4SPSPosDistribution *GPS_ParticlePosition;
	G4SPSAngDistribution *GPS_ParticleMomentum;

	TFile *rootFile;
	TH1F *theHistogram;


};

#endif /* PROTONDETECTORPRIMARYGENERATORACTION_HH_ */
