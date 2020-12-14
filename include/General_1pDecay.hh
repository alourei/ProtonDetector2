//--------------------------------------------------------------
//--------------------------------------------------------------
//
// General_1pDecay.hh 
//
// Description - a Discrete Process in G4 used to model Breakup
//               into Recoil(A,Z) + 1p
//
// updated 30 October 2008
//
// a class derived from G4VDiscreteProcess
//
// Written by: Brian Roeder, LPC Caen, email - roeder@lpccaen.in2p3.fr
// Modified :  Brian Roeder, TAMU, email - broeder@comp.tamu.edu 
//
// Program based on treatment given in bzb by JL Lecouey and 
// FM Marques. Also used certain functions written by them.
// Modified for use with C++ and GEANT4 by BTR.
//
//---------------------------------------------------------------
//---------------------------------------------------------------
// Call in your physics list with :
//
// #include "G4C9.hh"
// #include "General_1pDecay.hh"
//
//  pManager = G4C9::C9()->GetProcessManager();
//
//   G4double Rel_Eng = 0.5*MeV; // Breakup Eng of C9*  
//   G4double C9_Mass = 9.*931.494*MeV;  // Nuclear Mass of C9   
//
//   G4String theB8pProcessName = "General_1pDecay";
//   General_1pDecay* theB8pProcess = new General_1pDecay(theC8pProcessName);
//   theB8pProcess->SetRelativeEnergy("NotRandom",Rel_Eng,0.0);
//
//   theB8pProcess->SetGoldhaberParameters(true,9.,C9_Mass,90.*MeV);
//
//   pManager->AddDiscreteProcess(theB8pProcess);
//
//-------------------------------------------------------------
//-------------------------------------------------------------
//
// 8/23/07 - Breaks up 7He into 6He + n. Process is forced in
// GetMeanFreePath() (perhaps later can put in a short lifetime?)
// In PostStepDoIt(), we kill the 7He particle and produce the 
// two secondaries, 6He and n. 
// 
// 10/08/07 - Added Goldhaber and Japanese Kick functions, as in
// bzb. Use a kick of ~90 MeV/c. Set Beam/Kick in Physics List with
// SetGoldhaberParameters() function (see below).
// 
// 10/24/07 - Generalized so it can be used for other Frag->Rec+1n 
// systems. IsApplicable updated for Parent Nuclei 13Be and 16B in addition
// to 7He.
//
// 7 Oct. 2008 - Removed dependance on randCLHEP.h block. Random numbers 
// generated with G4UniformRand(). Added function randBW() to generate
// Breit-Wigner functions (instead of Delta functions).
//
// 30 Oct. 2008 - Modified into General_1pDecay, replacing neutrons
// with protons according to variable names, masses and secondaries
//

#ifndef General_1pDecay_hh 
#define General_1pDecay_hh 

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
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


class General_1pDecay : public G4VDiscreteProcess
{
public:

  // constructor
  General_1pDecay(const G4String& processName="General_1pDecay");

  // destructor
  ~General_1pDecay();

public:
  // These are the functions in the Process: Derived from G4VDiscreteProcess.hh

  G4bool IsApplicable(const G4ParticleDefinition& aParticle);
  // Decides if process applicable or not. Eventually make only valid for Neutrons
  // at some energy.

  G4double GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize,
                           G4ForceCondition* condition);
  // Overrides function in base class. Need to figure out how it plugs in!
  // Says returns MeanFreePath for a particle in a given material!
  // Invoked by Process Manager (Therefore necessary to get Process going!)

  G4VParticleChange* PostStepDoIt(const G4Track &aTrack, const G4Step &aStep);
  // This is the important function where you define the process
  // Returns steps, track and secondary particles
  // In this test, changes direction of neutron.
  // Invoked by Process Manager (Therefore necessary to get Process going!)

  void SetParentNucleus(G4int theMass, G4int theCharge);

  G4double Absolute(G4double num);
  G4double randBW(G4double center, G4double width);
  G4double VectModulus(G4ThreeVector Mom);
  G4double VectModulus(G4LorentzVector Mom);
  G4double ScalarProduct(G4ThreeVector V1, G4ThreeVector V2);
  G4double ScalarProduct(G4LorentzVector V1, G4ThreeVector V2);
  G4double ScalarProduct(G4LorentzVector V1, G4LorentzVector V2);

  G4LorentzVector LorentzBoost(G4LorentzVector Mom, G4ThreeVector Beta);

  // Produces BacktoBack recoil/neutron event in the center of mass.
  void BacktoBack(G4double DeltaM);

  void CMFrameToLabFrame();

    // Goldhaber Momentum Kick Functions - Added 2 Oct. 2007 - BTR

  G4double GaussianRandom(G4double FWHM);
  G4double GoldhaberFWHM(G4double A_Beam, G4double A_frag, G4double sigma0);
  void GoldhaberKick();

  // Use this function to enable Goldhaber Momentum Kick in Physics List.
  // Default setting has no kick! (cond = false) 
  // A_Beam is mass number of beam (e.g. 8. for 8He)
  // M_beam is mass of beam in units of MeV/c2
  void SetGoldhaberParameters(G4bool cond, G4double A_Beam, G4double M_beam, G4double FWHM);

  // Sets RelativeEnergy between particles (i.e. resonance energy)
  void SetRelativeEnergy(G4String Ran, G4double value, G4double Width);

private:

 // Hide assignment operator as private 
  General_1pDecay& operator=(const General_1pDecay &right);

  // Copy constructor
  General_1pDecay(const General_1pDecay&);

  // Other Pre-Defined Class Variables and Functions are below!

private:

  G4double Pi;

  G4bool GoldCond;     // enable Goldhaber if set to "true"
  G4double A_Beam;
  G4double Mass_beam;
  G4double Mom_FWHM;

  G4String RanCond;    // enable random E_rel (for acceptance) set to "random"
  G4double E_rel;
  G4double ResWidth;

  G4int A_Parent;
  G4int Z_Parent;

  // Storage Variables for Recoil
  G4int A_Rec;        // Mass Number of Recoil Nucleus
  G4int Z_Rec;        // Charge of Recoil Nucleus

  // Momentum 4-vectors for Fragment, recoil and neutron
  G4LorentzVector P_Frag;
  G4LorentzVector P_Recoil;
  G4LorentzVector P_Proton;

  // Mass Storage Variables
  G4int A_Frag;
  G4int Z_Frag;
  G4double Mass_Fragment;
  G4double Mass_Proton;
  G4double Mass_Recoil;

};
#endif
