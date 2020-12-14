/*
 * DPLBetaPlusDecay.hh
 *
 *  Created on: Feb 18, 2014
 *      Author: perezlou
 */

#ifndef DPLBETAPLUSDECAY_HH_
#define DPLBETAPLUSDECAY_HH_

#include <G4VRestDiscreteProcess.hh>

#include "globals.hh"
#include "G4VParticleChange.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4VParticleChange.hh"

#include "G4Material.hh"
#include "G4UnitsTable.hh"


class DPLBetaPlus_Decay: public G4VRestDiscreteProcess {
public:
	DPLBetaPlus_Decay(const G4String& processName);
	~DPLBetaPlus_Decay();

	  // These functions are the virtual functions overriden from G4VRestDiscreteProcess

	 G4bool IsApplicable(const G4ParticleDefinition& aParticle);
	  // Decides if process applicable or not.

	 G4double GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize,
	                           G4ForceCondition* condition);
	  // Overrides function in base class.
	  // Returns MeanFreePath for a particle in a given material!
	  // Invoked by Process Manager.
	  // Setting to Infinite (DBL_MAX) for first version

	  G4VParticleChange* PostStepDoIt(const G4Track &aTrack, const G4Step &aStep);
	  // This is the important function where you define the process
	  // Returns steps, track and secondary particles
	  // Invoked by Process Manager.
	  // Does nothing in first version. Later perhaps add decay in-flight?

	  G4double GetMeanLifeTime(const G4Track& aTrack, G4ForceCondition* condition);
	  // Overrides function in base class.
	  // Define Parent Nucleus MeanLifeTime from known or measured half-life

	  G4VParticleChange* AtRestDoIt(const G4Track &aTrack, const G4Step &aStep);
	  // This is the important function where you define the process
	  // Returns steps, track and secondary particles
	  // In this version, uses Phase3 to generate Daughter Nucleus, positron and neutrino!

	  void SetBetaDecayData(G4String theFileName, G4String theType, G4double theQValue, G4double theHalfLife);  // Default = 95ms

	  void SetParentNucleus(G4int theMass, G4int theCharge); // Default = 23Al

	  void SetQvalue_Daughter(G4double value){QValue_daughter= value;} // Default = 23Al

          void SetMeanLifeTime(G4double tao){TaoMeanLifeTime=tao;}
          void SetMeanLifeTimeDaughter(G4double tao){TaoMeanLifeTime_Daughter=tao;}

private:

	  G4double Pi;

	  G4double QValue;      // Calculates internally if daughter nucleus is excited
	  G4double QValue_daughter;      // Calculates internally if daughter nucleus is excited

	  G4String DecayType;   // Can be = "Branch" (for multiple decays with branching ratios)
	                        // or "Single"

	  G4double TaoMeanLifeTime; // tao = 1/lambda
	  G4double TaoMeanLifeTime_Daughter; // tao = 1/lambda

	  G4int A_Parent;     // Stores Mass of Parent Nucleus as set by user in Physics List
	  G4int Z_Parent;     // Stores Charge of Parent Nucleus as set by user in Physics List

	  G4int A_Part;        // Stores Mass of Parent Nucleus
	  G4int Z_Part;        // Stores Charge of Parent Nucleus
	  G4int A_Rec;         // Stores Mass of Daughter Nucleus
	  G4int Z_Rec;         // Stores Charge of Daughter Nucleus

	  G4double Mass_Recoil;
	  G4double Mass_Positron;
	  G4double Mass_Neutrino;

	  G4LorentzVector P_Recoil;
	  G4LorentzVector P_Positron;
	  G4LorentzVector P_Neutrino;

	  G4int NumberOfLines;
	  G4double Normalization;
	  G4double *DaughterExEng;
	  G4double *BetaIntensityRaw;
	  G4double *ProbRaw;
	  G4double *ProbLimit;

	  DPLBetaPlus_Decay& operator=(const DPLBetaPlus_Decay &right);
	  // Hide normal constructor
	  DPLBetaPlus_Decay();
	  // Copy constructor
	  DPLBetaPlus_Decay(const DPLBetaPlus_Decay&);

	  G4double CalculateMeanLifeTime(G4double HalfLife); // Calculates lambda MeanLifeTime

	  G4double GetDaughterExEng(); // Randomly chooses Daughter ExEng

	  G4double Absolute(G4double Num);      // Gives absolute value of number
	  G4double VectModulus(G4ThreeVector Mom);
	  G4double ScalarProduct(G4LorentzVector V1, G4ThreeVector V2);

	  G4LorentzVector LorentzBoost(G4LorentzVector Mom, G4ThreeVector Beta);

	  G4double dNdp3(G4double Mass1, G4double Mass2, G4double Mass3, G4double Eng, G4double P3);

	  void ThreeBodyPhaseSpace(G4double Qvalue);


};

#endif /* DPLBETAPLUSDECAY_HH_ */
