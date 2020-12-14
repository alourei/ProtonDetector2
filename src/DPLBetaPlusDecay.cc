/*
 * DPLBetaPlusDecay.cc
 *
 *  Created on: Feb 18, 2014
 *      Author: perezlou
 */

#include "DPLBetaPlusDecay.hh"
#include <G4TouchableHandle.hh>
#include <cmath>


DPLBetaPlus_Decay::DPLBetaPlus_Decay(const G4String& processName):G4VRestDiscreteProcess(processName) {

	// Default 23Al


	Pi = CLHEP::pi;

	G4double HalfLife = 470.*ms;
	G4double HalfLifeDaughter = 2.572*s;
	QValue = 12.243*MeV;   // 23Al g.s. -> 23Mg g.s. Q-value

	TaoMeanLifeTime = CalculateMeanLifeTime(HalfLife);
	TaoMeanLifeTime_Daughter =CalculateMeanLifeTime(HalfLifeDaughter) ; // tao = 1/lambda

	DecayType = "Single";

	A_Parent = 23;
	Z_Parent = 13;

	G4cout << "Constructor for DPLBetaPlus_Decay process started! " << G4endl;
	G4cout << "An AtRest Model for 3-body Beta Plus Decay" << G4endl;

}

DPLBetaPlus_Decay::~DPLBetaPlus_Decay() {

	  delete [] DaughterExEng;
	  delete [] BetaIntensityRaw;
	  delete [] ProbRaw;
	  delete [] ProbLimit;

}

void DPLBetaPlus_Decay::SetBetaDecayData(G4String theFileName, G4String theType, G4double theQValue, G4double theHalfLife){

	 TaoMeanLifeTime = CalculateMeanLifeTime(theHalfLife);

	  DecayType = theType;

	  QValue = theQValue;

	  if(DecayType != "Branch")
	    {
	      return;
	    }

	  G4cout << "The Half-Life T1/2 is now = " << G4BestUnit(theHalfLife,"Time") << G4endl;
	  G4cout << "The MeanLifeTime Tao is now = " << G4BestUnit(TaoMeanLifeTime,"Time") << G4endl;
	  G4cout << "The Decay scheme type for Al23 is set to : " << DecayType << G4endl;
	  G4cout << "The QValue for Z="<<Z_Parent<<"A="<<A_Parent<<" is set to = " << G4BestUnit(QValue,"Energy") << G4endl;

	  system("echo $G4BETADECAYSIDETAL23");

	  if(!getenv("G4BETADECAYSIDETAL23"))
	     {
	       G4cerr << "Please setenv G4BETADECAYSIDETAL23 to point to the Z="<<Z_Parent<<"A="<<A_Parent<<"branching data file" << G4endl;
	 	  G4ExceptionDescription description;
	 	  description << "Environment variable G4BETADECAYSIDETAL23 not set! ";
	       G4Exception("DPLBetaPlus_Decay::SetBetaDecayData()","DecayData_001",FatalException,description);
	     }
	   G4String  DirName = getenv("G4BETADECAYSIDETAL23");

	   G4String FileName = DirName+"/"+theFileName;
	   std::fstream theFile;
	   theFile.open(FileName, std::fstream::in );
	   theFile >> NumberOfLines;

	   DaughterExEng = new G4double[NumberOfLines];
	   BetaIntensityRaw = new G4double[NumberOfLines];
	   ProbRaw = new G4double[NumberOfLines];
	   ProbLimit = new G4double[NumberOfLines];

	   if(theFile.good())
	      {
	        G4cout << "Loading Data For FileName = " << FileName << G4endl;
	        for(G4int i=0; i<NumberOfLines; i++)
	          {
	   	 G4double theEnergy = 0.;    // Nucleus (23Mg) excitation eng.
	   	 G4double theIntensity = 0.;  // absolute (not in percent)

	   	 theFile >> theEnergy >> theIntensity;

	   	 G4cout << theEnergy << "  " << theIntensity << G4endl;

	   	 DaughterExEng[i] = theEnergy;
	   	 BetaIntensityRaw[i] = theIntensity;
	   	 Normalization += theIntensity;
	          }
	        G4cout << "Successfully Loaded !  Normalization = " << Normalization << G4endl;
	        theFile.close();

	        for(G4int i=0; i<NumberOfLines; i++)
	          {
	   	 ProbRaw[i] = BetaIntensityRaw[i]/Normalization;
	          }
	      }

	   else
	     {
	       G4cerr << "File = " << FileName << " not found or in improper format." << G4endl;
	    	 G4ExceptionDescription description;
	  	 description<< "File not found!!";
	       G4Exception("DPLBetaPlus_Decay::SetBetaDecayData()","DecayData_002",FatalException,description);
	   }

	      // Setup Limits and Distances

	      for(G4int i=0; i<NumberOfLines; i++)
		{
		  if(i == 0)
		    {ProbLimit[i] = ProbRaw[i];}
		  else
		    {ProbLimit[i] = ProbLimit[i-1]+ProbRaw[i];}
		}
}

G4double DPLBetaPlus_Decay::GetDaughterExEng(){

    G4double RecRandnum = G4UniformRand();
     // Search for Prob Interval ----

    G4int ProbInt = 0;

    for(G4int k=0; k<NumberOfLines; k++)
      {
	 if(k == 0 && RecRandnum < ProbLimit[k] && ProbRaw[k] > 0.)
	   { ProbInt = k; }
	 else if( k > 0 && RecRandnum >= ProbLimit[k-1] && RecRandnum < ProbLimit[k] && ProbRaw[k] > 0.)
	   { ProbInt = k; }
      }

    G4double EnergyOut = DaughterExEng[ProbInt];

 return EnergyOut;

}

G4double DPLBetaPlus_Decay::CalculateMeanLifeTime(G4double theHalfLife){

	  G4double Lambda = log(2.)/theHalfLife;

	  G4double Tao = 1./Lambda;

	  return Tao;
}

void DPLBetaPlus_Decay::SetParentNucleus(G4int theMass,G4int theCharge){

	  A_Parent = theMass;
	  Z_Parent = theCharge;

	  G4cout << "The Parent Nucleus for the Beta Decay is A = " << A_Parent << " and Z = "
		 << Z_Parent << G4endl;

}

G4double DPLBetaPlus_Decay::Absolute(G4double Num){
	 if(Num < 0.)
	    {Num *= -1.;}

	 return Num;
}

G4double DPLBetaPlus_Decay::VectModulus(G4ThreeVector Mom){

	  // Gives Modulus (normal) of a 3-vector
	  G4double V_Modulus = sqrt(pow(Mom[0],2)+pow(Mom[1],2)+pow(Mom[2],2));
	  return V_Modulus;
}

G4double DPLBetaPlus_Decay::ScalarProduct(G4LorentzVector V1,G4ThreeVector V2){

	 // Gives Scalar Product of a Lorentz Vector and a ThreeVector
	  G4double theScalarProduct = V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2];
	  return theScalarProduct;
}

G4LorentzVector DPLBetaPlus_Decay::LorentzBoost(G4LorentzVector Mom,G4ThreeVector Beta){

	 // Boost routine originally written in Fortran by F.Miguel Marques
	  // "Adds" a Momentum vector in the rest frame to the Beta
	  // thereby giving a Lorentz Transformation of two 4-vectors
	  // Translated to C++/GEANT4 by : BTR

	  G4double Norm_Beta = VectModulus(Beta);

	  if(Norm_Beta == 1.)
	    {
	      G4cout << " Error! Vector |Beta| is at Speed of Light! - Norm_Beta = " << Norm_Beta << G4endl;
	      exit(0);
	    }
	  else if(Norm_Beta == 0.)
	    {
	      G4cout << "No Boost, |Beta| = 0" << G4endl;
	      return Mom;
	    }

	  G4double gamma = 1./sqrt(1.-pow(Norm_Beta,2));    // Calculates gamma factor
	  // G4cout << "gamma = " << gamma << G4endl;

	  // Now add the Lorentz Vectors by Scalar Dot Product
	  // and Matrix Transform

	  G4double theScalarAmp = ScalarProduct(Mom,Beta);
	  G4double P_Parallel = theScalarAmp/Norm_Beta;
	  G4double E_total = Mom[3];

	  G4ThreeVector P_Perp;

	  for(G4int i=0; i<3; i++)
	    {P_Perp[i] = Mom[i]-(P_Parallel*Beta[i]/Norm_Beta);}

	  P_Parallel = gamma*(P_Parallel+(Norm_Beta*E_total));
	  Mom[3] = gamma*(E_total+theScalarAmp);
	  // Update Output Momentum Vector

	  for(G4int i=0; i<3; i++)
	    {Mom[i] = P_Perp[i]+(P_Parallel*Beta[i]/Norm_Beta);}

	  return Mom;
}

G4double DPLBetaPlus_Decay::dNdp3(G4double Mass1,G4double Mass2,G4double Mass3,G4double Eng,G4double P3){

	  // Momentum Probability Distribution Function

	  G4double Eng3 = sqrt(pow(Mass3,2)+pow(P3,2));

	  G4double Fact1 = pow(P3,2)/Eng3;

	  G4double Fact2 = (pow(Eng,2)+pow(Mass3,2)-2.*Eng*Eng3-pow((Mass1+Mass2),2))*
			   (pow(Eng,2)+pow(Mass3,2)-2.*Eng*Eng3-pow((Mass1-Mass2),2));

	  G4double Denom = pow(Eng,2)+pow(Mass3,2)-2.*Eng*Eng3;

	  G4double dNdp3 = Fact1*sqrt(Absolute(Fact2))/Denom;

	  return dNdp3;
}

void DPLBetaPlus_Decay::ThreeBodyPhaseSpace(G4double Q){

	  const G4int NDiv = 200;

	  G4double Mass1 = Mass_Recoil;
	  G4double Mass2 = Mass_Positron;
	  G4double Mass3 = Mass_Neutrino;

	  // Four - Vectors for Phase-Space Generation
	  G4LorentzVector P1;
	  G4LorentzVector P2;
	  G4LorentzVector P3;
	  G4ThreeVector Beta;

	  G4double Energy = Mass1+Mass2+Mass3+Q;

	  G4double P3M = sqrt(Absolute((pow(Energy,2)-pow((Mass1+Mass2+Mass3),2))*
				       (pow(Energy,2)-pow((Mass1+Mass2-Mass3),2))))/(2.*Energy);

	  // Generate Distribution Function

	  G4double F_Old = 0.;
	  G4double P_Old = 0.;
	  G4double StepSize = P3M/(static_cast<G4double>(NDiv));

	  G4double P_Max = P3M - StepSize;
	  G4double F_Max = dNdp3(Mass1,Mass2,Mass3,Energy,P_Max);

	  while(F_Max > F_Old)
	    {
	      P_Old = P_Max;
	      F_Old = F_Max;
	      P_Max = P_Old-StepSize;
	      F_Max = dNdp3(Mass1,Mass2,Mass3,Energy,P_Max);
	    }
	  F_Max = dNdp3(Mass1,Mass2,Mass3,Energy,(P_Max+StepSize/2.));

	  G4double P3_Mod;
	  G4double Func3;

	 do{
		P3_Mod = P3M*G4UniformRand();
	        Func3 = F_Max*G4UniformRand();
	    }
	 while(Func3 > dNdp3(Mass1,Mass2,Mass3,Energy,P3_Mod));

	 // Generate an isotropic theta, phi
	  G4double theta = acos(1.-2.*G4UniformRand());
	  G4double phi = 2.*Pi*G4UniformRand();

	  // Create Momentum vectors for Particle3

	  P3[0] = P3_Mod*sin(theta)*cos(phi);
	  P3[1] = P3_Mod*sin(theta)*sin(phi);
	  P3[2] = P3_Mod*cos(theta);
	  P3[3] = sqrt((pow(Mass3,2))+(pow(P3_Mod,2)));

	  // Continue with a two-body phase space in (12) center of mass

	  G4double Mass12 = sqrt(pow((Energy-P3[3]),2) - pow(P3_Mod,2));

	  G4double P1_Mod = sqrt(Absolute((pow(Mass12,2)-pow((Mass1+Mass2),2))*
	                                  (pow(Mass12,2)-pow((Mass1-Mass2),2))))/(2.*Mass12);

	  theta = acos(1.-2.*G4UniformRand());
	  phi = 2.*Pi*G4UniformRand();

	  P1[0] = P1_Mod*sin(theta)*cos(phi);
	  P1[1] = P1_Mod*sin(theta)*sin(phi);
	  P1[2] = P1_Mod*cos(theta);
	  P1[3] = sqrt((pow(Mass1,2))+(pow(P1_Mod,2)));

	  for(G4int i=0;i<3;i++)
	    {P2[i] = - P1[i];}
	  P2[3] = sqrt((pow(Mass2,2))+(pow(P1_Mod,2)));

	  // Now put V1 and V2 into frame of V3

	  for(G4int i=0;i<3;i++)
	    { Beta[i] = -P3[i]/(Energy-P3[3]); }

	  G4LorentzVector P1_Final;
	  G4LorentzVector P2_Final;
	  G4LorentzVector P3_Final;

	  P1_Final = LorentzBoost(P1,Beta);
	  P2_Final = LorentzBoost(P2,Beta);
	  P3_Final = P3;
	  // Do nothing with P3_Final (P3)

	  // Update Particle Vectors
	  P_Recoil = P1_Final;
	  P_Positron = P2_Final;
	  P_Neutrino = P3_Final;
}

G4bool DPLBetaPlus_Decay::IsApplicable(const G4ParticleDefinition& particle){

	  // All particles, other than nuclei, are rejected by default.
	  G4String ParticleName = particle.GetParticleName();

	  if (particle.GetParticleType() == "nucleus")
	    {
	    G4cout << "IsApplicable() : Particle " << ParticleName
	           << " valid for BetaPlus3BodyDecay Process!" << G4endl;
	    return true;
	    }
	  else
	    {
	      G4cout << "IsApplicable() : DPLBetaPlus_Decay can not be applied to " << ParticleName << G4endl;
	      G4cout << "This Particle Type is not a Nucleus! " << G4endl;
	      return false;
	    }
	  return false;
}

G4double DPLBetaPlus_Decay::GetMeanFreePath(const G4Track &aTrack,G4double previousStepSize,G4ForceCondition *condition){

	  // This Process is not called unless this is a "Discrete Process"
	  // Set to Do Nothing in that case.
	 G4cout << "DPLBetaPlus_Decay GetMeanFreePath called!" << G4endl;

	 *condition=NotForced;
	  //*condition=Forced;  // Forces Process to Occur.
	  previousStepSize = 0.1*mm;

	 return DBL_MAX;
}

G4VParticleChange *DPLBetaPlus_Decay::PostStepDoIt(const G4Track &aTrack,const G4Step &aStep){

	  // Needed only if GetMeanFreePath != DBL_MAX !
	  // This Process is not called unless this is also called as a "Discrete Process"

	  G4cout << "PostStepDoIt called !" << G4endl;
	  G4cout << "Does Nothing!! " << G4endl;

	 return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);

}


G4double DPLBetaPlus_Decay::GetMeanLifeTime(const G4Track &aTrack, G4ForceCondition *condition){

	*condition=NotForced;
	 //*condition=Forced;  // Forces Process to Occur.

	 // G4cout << "DPLBetaPlus_Decay GetMeanLifeTime called!" << G4endl;
	 // G4cout << "Tao was : " << G4BestUnit(TaoMeanLifeTime,"Time") << G4endl;

	  const G4DynamicParticle* Particle = aTrack.GetDynamicParticle();
	  const G4ParticleDefinition* ParticleDef = Particle->GetDefinition();

	  Z_Part = static_cast<int>(ParticleDef->GetPDGCharge());   // Returns Charge of Particle

	  /*
	  if(Z_Part == Z_Parent)
	    {return TaoMeanLifeTimeAl23;}
	  else
	    {return DBL_MAX;}
	  */

	  // For a daughter decay included!
	  if(Z_Part == Z_Parent || Z_Part == (Z_Parent-1))
	    {return TaoMeanLifeTime;}
	  else
	    {return DBL_MAX;}


	 return DBL_MAX;
	}

G4VParticleChange *DPLBetaPlus_Decay::AtRestDoIt(const G4Track &aTrack, const G4Step &aStep){

	  const G4DynamicParticle* Particle = aTrack.GetDynamicParticle();
	  const G4ParticleDefinition* ParticleDef = Particle->GetDefinition();


	  A_Part = ParticleDef->GetBaryonNumber();                  // Returns int. mass number of Particle
	  Z_Part = static_cast<int>(ParticleDef->GetPDGCharge());   // Returns Charge of Particle

	  G4String theParticleName = ParticleDef->GetParticleName();

	  G4StepPoint* thePostStepPoint = aStep.GetPostStepPoint();
	  G4ThreeVector thePosition = thePostStepPoint->GetPosition();

	  //G4double Particle_KE = Particle->GetKineticEnergy();
	  G4double GlobalTime = aTrack.GetGlobalTime();
	  //G4double LocalTime = aTrack.GetLocalTime();
	  G4double ProperTime = aTrack.GetProperTime();

	  G4TouchableHandle theTouchable = aTrack.GetTouchableHandle();

	  // 6/16/09 - Modified for use with geant4.9.2.p01
	  // Use SetParentNucleus() method inside physics list to choose which" generic ion"
	  // beta decays. All other ions do nothing

	  // 1/19/2011
	  // Modify this "if" block also if decaying the daughter nuclei!
	  //if(Z_Part > Z_Parent || Z_Part < (Z_Parent))
	  if(Z_Part > Z_Parent || Z_Part < (Z_Parent-1))
	    {
	      // Stop the At Rest Ion - Do Nothing
	      aParticleChange.ProposeTrackStatus(fStopButAlive);
	      aParticleChange.ProposeEnergy(0.);
	      aParticleChange.ProposePosition(thePosition);
	      aParticleChange.ProposeGlobalTime(GlobalTime);
	      aParticleChange.ProposeProperTime(ProperTime);
	      return pParticleChange;
	    }

	  // 3/2/2011
	  // Beta Decay only 30% of the 23Mg, else do nothing (due to lifetime)!
	  if(Z_Part == 12 && A_Part == 23)
	    {
	      G4double Mg23RandChoice = G4UniformRand();
	      if(Mg23RandChoice > 0.25)
		{
		  // Stop the At Rest Ion - Do Nothing
		  aParticleChange.ProposeTrackStatus(fStopAndKill);
		  aParticleChange.ProposeEnergy(0.);
		  aParticleChange.ProposePosition(thePosition);
		  aParticleChange.ProposeGlobalTime(GlobalTime);
		  aParticleChange.ProposeProperTime(ProperTime);
		  return pParticleChange;
		}
	      else
		{
		  // Continue
 		}
	    }


	  // 3/15/2010
	  // Custom Beta Decay for Al23 --------------------------
	  // Calculate Transition Energy based on Branching Ratios, etc.

	  // Below produces Daughter nucleus in various excitation states
	  // according to intensities given in the loaded data file.
	  // Apply further physics processes for further decays after this one!

	  // 2/17/2010 - Had Error in Transition Energy - for Beta+ decay,
	  // need to account for 2Me mass to create positron out.
	  // Corrected here so can continue to use Qvalue(g.s.)->(g.s.) in PhysicsList!

	  G4double TransEng = QValue-(2.*0.511*MeV);  // same for single energy Beta Decay
	  G4double ExEngRec = 0.*MeV;

	  if(DecayType== "Branch" && Z_Part == Z_Parent)
	    {
	      ExEngRec = GetDaughterExEng();  // Returns excitation energy of Daugther Nucleus
	      TransEng = QValue - (2.*0.511*MeV) - ExEngRec;
	    }
	  else
	    {
	      // 01/19/2011 - Include also daughter decay for 23Mg
	      ExEngRec = 0.;
	      TransEng = QValue_daughter - (2.*0.511*MeV); // For 23Mg g.s. -> 23Na g.s.
	    }

	  // For Beta+ Decay, A_Rec = A_Part, Z_Rec = Z_Part-1.
	 A_Rec = A_Part;
	 Z_Rec = Z_Part-1;

	 // Generate Daughter Nucleus
	    G4DynamicParticle* theRec = new G4DynamicParticle;
	    // GetIon() Method Arguements are GetIon(Charge,Mass,ExcitationEng)
	    G4ParticleDefinition* theRecDefinition = G4ParticleTable::GetParticleTable()->GetIon(Z_Rec,A_Rec,ExEngRec);
	    theRec->SetDefinition(theRecDefinition);
	    Mass_Recoil = theRec->GetMass();
	    //G4cout << "The Recoil Mass is : " << Mass_Recoil << " MeV " << G4endl;

	// Generate a Positron
	    G4DynamicParticle* thePositron = new G4DynamicParticle;
	    G4ParticleDefinition* thePositronDefinition = G4Positron::Positron();
	    thePositron->SetDefinition(thePositronDefinition);
	    Mass_Positron = thePositron->GetMass();

	// Generate a Neutrino
	    G4DynamicParticle* theNeutrino = new G4DynamicParticle;
	    G4ParticleDefinition* theNeutrinoDefinition = G4NeutrinoE::NeutrinoE();
	    theNeutrino->SetDefinition(theNeutrinoDefinition);
	    Mass_Neutrino = theNeutrino->GetMass();

	 // Generate 4-Vectors for Products using ThreeBodyPhaseSpace Routine

	    ThreeBodyPhaseSpace(TransEng);

	// Kill the Parent nucleus track
	   aParticleChange.ProposeTrackStatus(fStopAndKill); //Kills Parent Nucleus
	   aParticleChange.ProposeEnergy(0.);
	   aParticleChange.ProposePosition(thePosition);

	// Set Momentrum 4-Vectors for Products

	   theRec->Set4Momentum(P_Recoil);
	   thePositron->Set4Momentum(P_Positron);
	   theNeutrino->Set4Momentum(P_Neutrino);

	// Set Final Time according to TaoMeanLifeTime

	   G4double RanTime = G4UniformRand();
	   if(RanTime < 1.e-5) // Prevents errors with log(0)
	     {RanTime = 1.e-5;}
	   G4double TimeAtDecay;

	   if(Z_Part == Z_Parent)
	     {TimeAtDecay = (-log(RanTime)*TaoMeanLifeTime);}

	   GlobalTime += TimeAtDecay;

	 // Set final Tracks in Motion!
	    G4Track* theRecTrack = new G4Track(theRec,GlobalTime,thePosition);
	    theRecTrack->SetTouchableHandle(theTouchable);

	    G4Track* thePositronTrack = new G4Track(thePositron,GlobalTime,thePosition);
	    thePositronTrack->SetTouchableHandle(theTouchable);

	    G4Track* theNeutrinoTrack = new G4Track(theNeutrino,GlobalTime,thePosition);
	    theNeutrinoTrack->SetTouchableHandle(theTouchable);

	    aParticleChange.SetNumberOfSecondaries(3);
	    aParticleChange.AddSecondary(theRecTrack);
	    aParticleChange.AddSecondary(thePositronTrack);
	    aParticleChange.AddSecondary(theNeutrinoTrack);

	    // 01/19/2011 --- If you don't follow the recoils, only interested in the betas
	    /*
	    aParticleChange.SetNumberOfSecondaries(2);
	    aParticleChange.AddSecondary(thePositronTrack);
	    aParticleChange.AddSecondary(theNeutrinoTrack);
	    */

	  return pParticleChange;

}
