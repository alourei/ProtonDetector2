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
//
// $Id: BTRBetaPlus_CustomDecayAl23.hh,  2010/03/15 21:07:56 roeder Exp $
// GEANT4 tag $Name: geant4-09-02.p01 $
//
//-----------------------------------------------------------------------
// Name : BTRBetaPlus_CustomDecayAl23.hh
//
// Description : A GEANT4 VRestDiscreteProcess that Generates Beta+
//               Decays in a given nucleus. Does not include EC.
//               -> This version -> customized for 23Al Beta Decay
//
// Method      : Uses "Phase3" 3 body phase space from FM Marques to
//               decay a parent nucleus into the daughter + e+ + v.
//               -> This version -> includes beta+ decay branching 
//               ratios to 3 20Na states.
//
// Assumption  : v1: Assume lifetime of parent nucleus is >> flight time.
//               Thus we assume MeanFreePath infinite (DBL_MAX) and
//               no PostStepDoIt process.
//               Impliment Beta Decay in MeanLifeTime Process and 
//               AtRestDoIt processes.
//
// First Version : 15 March 2010 - BTR
//-----------------------------------------------------------------------
//
// 16 Jun 2009 - geant4.9.2.p01 -> Modified AtRestDoIt so that only the 
// user-defined parent nucleus Beta Decays. Others are left in status 
// "fStopButAlive" so that other processes can be applied to them.
// This fixes the problem in geant4.9.2 where the ionIonisation E&M process
// can only be used with geant4 defined "generic ions" and not with user
// defined particles.
//
// Usage (in physics list) :
// G4String BetaDecayProcessName = "BTRBetaPlusDecay";
// BTRBetaPlus_CustomDecayAl23* theBetaPlusDecay = new BTRBetaPlus_CustomDecayAl23(BetaDecayProcessName);
// theBetaPlusDecay->SetParentNucleus(23,13);            // e.g. 23Al (A=23,Z=13)
// theBetaPlusDecay->SetBetaDecayData("Branch",0.,470.*ms); // (Type,Single Beta Energy, Mean Life Time)
// pManager->AddRestProcess(theBetaPlusDecay);
//
//-----------------------------------------------------------------------
// version - 3/15/2010 - BTR
//
//


#ifndef BTRBetaPlus_CustomDecayAl23_hh 
#define BTRBetaPlus_CustomDecayAl23_hh 

#include "globals.hh"
#include "G4VRestDiscreteProcess.hh"
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


class BTRBetaPlus_CustomDecayAl23 : public G4VRestDiscreteProcess
{
public:
  // constructor
  BTRBetaPlus_CustomDecayAl23(const G4String& processName);

  // destructor
  ~BTRBetaPlus_CustomDecayAl23();

public:
  // These functions are the virtual functions overriden from G4VRestDiscreteProcess

 G4bool IsApplicable(const G4ParticleDefinition& aParticle);
  // Decides if process applicable or not. Eventually make only valid for Neutrons
  // at some energy. 
 
  G4double GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize,
                           G4ForceCondition* condition);
  // Overrides function in base class. 
  // Returns MeanFreePath for a particle in a given material!
  // Invoked by Process Manager.
  // Setting to Infinite (DBL_MAX) for first version

  G4VParticleChange* PostStepDoIt(const G4Track &aTrack, const G4Step &aStep);
  // This is the important function where you define the process
  // Returns steps, track and secondary particles
  // In this test, changes direction of neutron.
  // Invoked by Process Manager.
  // Does nothing in first version. Later perhaps add decay in-flight?

  G4double GetMeanLifeTime(const G4Track& aTrack, G4ForceCondition* condition);
  // Overrides function in base class. 
  // Define Parent Nucleus MeanLifeTime from known or measured half-life

  G4VParticleChange* AtRestDoIt(const G4Track &aTrack, const G4Step &aStep);
  // This is the important function where you define the process
  // Returns steps, track and secondary particles
  // In this version, uses Phase3 to generate Daughter Nucleus, positron and neutrino!

  void SetBetaDecayDataAl23(G4String theFileName, G4String theType, G4double theQValue, G4double theHalfLife);  // Default = 95ms

  void SetParentNucleus(G4int theMass, G4int theCharge); // Default = 23Al

private:

  // Hide assignment operator as private 
  BTRBetaPlus_CustomDecayAl23& operator=(const BTRBetaPlus_CustomDecayAl23 &right);
  // Hide normal constructor 
  BTRBetaPlus_CustomDecayAl23();
  // Copy constructor
  BTRBetaPlus_CustomDecayAl23(const BTRBetaPlus_CustomDecayAl23&);

private:

  // These functions are included in the class to generate 
  // MeanLifeTime
  // 3Body Phase Space (secondaries)

  G4double CalculateMeanLifeTime(G4double HalfLife); // Calculates lambda MeanLifeTime

  G4double GetDaughterExEngAl23(); // Randomly chooses Daughter ExEng

  G4double Absolute(G4double Num);      // Gives absolute value of number
  G4double VectModulus(G4ThreeVector Mom);
  G4double ScalarProduct(G4LorentzVector V1, G4ThreeVector V2);

  G4LorentzVector LorentzBoost(G4LorentzVector Mom, G4ThreeVector Beta);

  G4double dNdp3(G4double Mass1, G4double Mass2, G4double Mass3, G4double Eng, G4double P3);

  void ThreeBodyPhaseSpace(G4double Qvalue);

  // Class Variables

  G4double Pi;

  G4double QValueAl23;      // Calculates internally if daughter nucleus is excited

  G4String DecayTypeAl23;   // Can be = "Branch" (for multiple decays with branching ratios) 
                        // or "Single"

  G4double TaoMeanLifeTimeAl23; // tao = 1/lambda

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

  G4int NumberOfLinesAl23;
  G4double DaughterExEngAl23[16];	
  G4double BetaIntensityRawAl23[16];
  G4double NormalizationAl23;
  G4double ProbRawAl23[16];
  G4double ProbLimitAl23[16];

};
#endif
