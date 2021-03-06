/*
 * DPL1pdecay.cc
 *
 *  Created on: Feb 19, 2014
 *      Author: perezlou
 */

#include <DPL1pdecay.hh>
#include <cmath>
#include <iomanip>


DPL_1pDecay::DPL_1pDecay(const G4String& processName) : G4VDiscreteProcess(processName) {

	  Pi = CLHEP::pi;

	  A_Parent = 23; // Default is 23Mg
	  Z_Parent = 12;

	  GoldCond = false;
	  RanCond = "NotRandom";
	  E_rel = 0.500*MeV;    // Normal rel. eng., changed with SetRelativeEnergy Function
	  Excitation_Energy = 0.*MeV;    //Excitation energy of the daughter
	  ResWidth = 0.0;       // Default delta function, changed with SetRelativeEnergy Function

	  G4cout << "Constructor for DPL_1pDecay process was called! " << G4endl;
	  G4cout << "A simple relativistic model for 1 proton breakup!" << G4endl;
}

DPL_1pDecay::~DPL_1pDecay() {
	// TODO Auto-generated destructor stub
	   delete [] DaughterExEng;
	   delete [] IntensityRaw;
	   delete [] ProbRaw;
	   delete [] ProbLimit;

}

G4double DPL_1pDecay::Absolute(G4double num){
	// Replaces abs() command in c++ cmath library
	  if(num < 0.)
	    { num *= -1.;}
	  return num;
}

G4double DPL_1pDecay::randBW(G4double center=1.0,G4double width=0.2){

	// creates a random Breit-Wigner Distribution
   // values outside of 10*width thrown out
	  G4double highlimit = center+10.*width;
	  G4double lowlimit = center-10.*width;
	  G4double ranBW=center;
	  do
	    {ranBW = CLHEP::RandBreitWigner::shoot(center,width);}
	  while(ranBW < lowlimit || ranBW > highlimit);

	  return ranBW;
}

G4double DPL_1pDecay::VectModulus(G4ThreeVector Mom){

	// Gives Modulus (normal) of a 3-vector
	  G4double V_Modulus = sqrt(pow(Mom[0],2)+pow(Mom[1],2)+pow(Mom[2],2));
	  return V_Modulus;
}

G4double DPL_1pDecay::VectModulus(G4LorentzVector Mom){

	  // Gives Modulus (normal) of a 3-vector within a four vector
	  G4double V_Modulus = sqrt(pow(Mom[0],2)+pow(Mom[1],2)+pow(Mom[2],2));
	  return V_Modulus;

}

G4double DPL_1pDecay::ScalarProduct(G4ThreeVector V1,G4ThreeVector V2){

	// Gives Scalar Product of two three Vectors
	  G4double theScalarProduct = V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2];
	  return theScalarProduct;
}

G4double DPL_1pDecay::ScalarProduct(G4LorentzVector V1, G4ThreeVector V2){
  // Gives Scalar Product of a Lorentz Vector and a ThreeVector
  G4double theScalarProduct = V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2];
  return theScalarProduct;
}

G4double DPL_1pDecay::ScalarProduct(G4LorentzVector V1, G4LorentzVector V2){
  // Gives Scalar Product of two three Vectors within two Lorentz Vectors
  G4double theScalarProduct = V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2];
  return theScalarProduct;
}

G4LorentzVector DPL_1pDecay::LorentzBoost(G4LorentzVector Mom, G4ThreeVector Beta){
  // Boost routine originally written in Fortran by F.Miguel Marques
  // "Adds" a Momentum vector in the rest frame to the Beta
  // thereby giving a Lorentz Transformation of two 4-vectors
  // Translated to C++/GEANT4 by : BTR

  G4double Norm_Beta = VectModulus(Beta);

  if(Norm_Beta == 1.)
    {
      G4cout << " Error! Vector |Beta| is at Speed of Light! - Norm_Beta = " << Norm_Beta << G4endl;
	  G4ExceptionDescription description;
	  description << "Beta=1!";
      G4Exception("DPL_1pDecay::LorentzBoost()","Boost_001",EventMustBeAborted,description);
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

void DPL_1pDecay::BacktoBack(G4double delta_M)
{
  // C++ version of "Phase2()" written in Fortran by F.Miguel Marques
  // Updated slightly by BTR with C++/GEANT setup in mind *** 8/24/2007
  // uses CLHEP lib. Include randCLHEP.h header file.

  //G4cout << "Running BacktoBack() Function " << G4endl;

  G4double Mass1 = Mass_Proton;
  G4double Mass2 = Mass_Recoil;

  if(ResWidth > 0.)
    {delta_M = randBW(delta_M,ResWidth);}

  G4double Energy = Mass1+Mass2+delta_M;

  G4double P_Mag = sqrt(Absolute((pow(Energy,2)-pow((Mass1+Mass2),2))*(pow(Energy,2)-pow((Mass1-Mass2),2))))/(2.*Energy);

  // Generate an isotropic theta, phi
  G4double theta;
  //G4double OpenAngle = 2.6*deg;
  //theta = acos(1.+(cos(OpenAngle)-1.)*G4UniformRand());
  //theta = acos(1.-2.*G4UniformRand());
  theta =Pi/2;

  G4double phi = 2.*Pi*G4UniformRand();

  // Create Momentum vectors for Proton, Recoil in center of mass

  P_Proton[0] = P_Mag*sin(theta)*cos(phi);
  P_Proton[1] = P_Mag*sin(theta)*sin(phi);
  P_Proton[2] = P_Mag*cos(theta);
  P_Proton[3] = sqrt((pow(Mass1,2))+(pow(P_Mag,2)));

  for(G4int i=0; i<3; i++)
    {
     P_Recoil[i] = -P_Proton[i];
    }
  P_Recoil[3] = sqrt((pow(Mass2,2))+(pow(P_Mag,2)));

  // End of BacktoBack() function
}


void DPL_1pDecay::CMFrameToLabFrame()
{
  // Takes CM distribition in BacktoBack() and Boosts to Lab Frame
  // Follows routine written in Fortran for Bzb written by: JL Lecouey
  // Translated and updated for C++/GEANT by : BTR

  //G4cout << " The Frag Momentum Direction (P_Frag) is : " << P_Frag << G4endl;

  G4ThreeVector Beta_Frag;    // Turn Beam Momentum into normalized velocity vector
  for(G4int i=0;i<3;i++)
    {Beta_Frag[i] = P_Frag[i]/P_Frag[3];}

// Now add Beam Lorentz Vector to Momentum Vectors for the neutron and recoil

  //G4cout<<"Before BOOST!! "<<P_Proton <<" "<<P_Recoil<<" BETA "<<Beta_Frag<<G4endl;

  P_Proton = LorentzBoost(P_Proton, Beta_Frag);
  P_Recoil = LorentzBoost(P_Recoil, Beta_Frag);

  //G4cout<<"After BOOST!! "<<P_Proton <<" "<<P_Recoil<<" BETA "<<Beta_Frag<<G4endl;

}

void DPL_1pDecay::SetRelativeEnergy(G4String Ran, G4double value, G4double Width)
{
  RanCond = Ran;
  E_rel = value;
  ResWidth = Width;
}


// Functions for Goldhaber Momentum Kick - Added 2 Oct. 2007, BTR

G4double DPL_1pDecay::GaussianRandom(G4double FWHM)
{
  // Translated Random Gaussian Function Generator
  // center is at 0.
  // Originally written by FM Marques

  G4double GRand = 0.;

  if(FWHM > 0.)
    {
      G4double r1 = G4UniformRand();
      G4double r2 = G4UniformRand();
      if(r2 == 0.)
	{ r2 = 1e-9; }
      G4double cons = sqrt(-2.*log(r2));
      G4double x = cons*cos(2.*Pi*r1);
      GRand = FWHM*x/2.3548;
    }
  return GRand;
}

G4double DPL_1pDecay::GoldhaberFWHM(G4double A_part, G4double A_frag, G4double sigma0)
{
  // Masses used here should be integer mass numbers of beam and fragment!

  G4double fraction = (A_frag*(A_part-A_frag))/(A_part-1.);

  G4double GoldFWHM = 0.;

  if(fraction > 0.)
    {GoldFWHM = 2.3548*sigma0*sqrt(fraction);}
  else
    {GoldFWHM = 0.;}

  return GoldFWHM;
}

void DPL_1pDecay::GoldhaberKick()
{
  G4double GFWHM = GoldhaberFWHM(A_Beam,A_Frag,Mom_FWHM);

	G4LorentzVector P_Gold;

	for(G4int i=0;i<3;i++)
	  {
            P_Gold[i] = GaussianRandom(GFWHM);
	    //P_Gold[i] = randGauss(0.,GFWHM);
	  }

	P_Gold[3] = sqrt(pow(P_Gold[0],2)+pow(P_Gold[1],2)+pow(P_Gold[2],2)+pow((Mass_Fragment+E_rel),2));
        //G4cout << "P_Gold = " << P_Gold << G4endl;

	// Normalize a three vector
	G4ThreeVector B3_Gold;
	for(G4int i=0;i<3;i++)
	  {B3_Gold[i] = P_Gold[i]/P_Gold[3];}

	//G4cout << "B3_Gold = " << B3_Gold << G4endl;

        // Boost Particle Vector with Same Kick
	P_Proton = LorentzBoost(P_Proton,B3_Gold);
	P_Recoil = LorentzBoost(P_Recoil,B3_Gold);
}

void DPL_1pDecay::SetGoldhaberParameters(G4bool cond, G4double A_Mass_Num, G4double M_beam, G4double FWHM_set)
{
  GoldCond = cond;
  A_Beam = A_Mass_Num;
  Mass_beam = M_beam;
  Mom_FWHM = FWHM_set;

  if(GoldCond == true)
    {G4cout << "****************Goldhaber enabled!!!!" << G4endl;}
}

void DPL_1pDecay::SetParentNucleus(G4int theMass, G4int theCharge)
{
  A_Parent = theMass;
  Z_Parent = theCharge;

  G4cout << "The Parent Nucleus for the Proton Decay is A = " << A_Parent << " and Z = "
	 << Z_Parent << G4endl;
}


void DPL_1pDecay::SetProtonDecayData(G4String theFileName){

	G4cout<<"Getting Energies for protons"<<G4endl;

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
	   IntensityRaw = new G4double[NumberOfLines];
	   ProbRaw = new G4double[NumberOfLines];
	   ProbLimit = new G4double[NumberOfLines];
	   Normalization = 0.;

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
	   	 IntensityRaw[i] = theIntensity;
	   	 Normalization += theIntensity;
	          }
	        G4cout << "Successfully Loaded !  Normalization = " << Normalization << G4endl;
	        theFile.close();

	        for(G4int i=0; i<NumberOfLines; i++)
	          {
	   	 ProbRaw[i] = IntensityRaw[i]/Normalization;
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







void DPL_1pDecay::SetProtonDecayData2(G4String theFileName){

	G4cout<<"Getting Kinetic Energies for protons"<<G4endl;

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

	   //std::vector <G4double> ExEnergy;
	   G4double energy,proton,probability;
	   G4double previous_energy=0;
	   std::vector <G4double> temp1;
	   std::vector <G4double> temp2;
	   G4int j=0;

	   if(theFile.good())
	      {
		std::cout << "Loading Data For FileName = " << FileName << std::endl;
		while(theFile>>energy>>proton>>probability){
		  //cout<<j<<" "<<previous_energy<<" "<<energy<<" "<<proton<<" "<<probability<<endl;
		  if((energy == previous_energy )|| j==0){
		    if(ExEnergy.size()==0)
		      ExEnergy.push_back(energy);
		    //cout<<"HERE!"<<endl;
		    temp1.push_back(proton);
		    temp2.push_back(probability);
		    //cout<<temp1.at(j)<<endl;
		    j++;
		  }
		  else{
		    ExEnergy.push_back(energy);
		    Energy_table.push_back(temp1);      
		    Probability_table.push_back(temp2);
		    Limits_table.push_back(temp2);
		    temp1.clear();
		    temp2.clear();
		    temp1.push_back(proton);
		    temp2.push_back(probability);
		  }
		  //cout<<temp1.size()<<endl;
		  previous_energy=energy;
		}

	        theFile.close();

	      }

	   else
	     {
	       G4cerr << "File = " << FileName << " not found or in improper format." << G4endl;
	    	 G4ExceptionDescription description;
	  	 description<< "File not found!!";
	       G4Exception("DPL1p_Decay::SetProtonDecayData()","DecayData_002",FatalException,description);
	   }

	      // Setup Limits and Distances

	      ExEnergy.pop_back();

	      for(G4int j=0;j<ExEnergy.size();j++){
		for(G4int k=0;k<Energy_table[j].size();k++){ 
		  if(k == 0)
		    Limits_table[j][k] = Probability_table[j][k];
		  else
		    Limits_table[j][k] = Limits_table[j][k-1]+Probability_table[j][k];
		
		  std::cout<<ExEnergy[j]<<" "<< Energy_table[j][k]<<" "<<Probability_table[j][k]<<" "<<Limits_table[j][k]<<std::endl;
		}
		std::cout<<std::endl;
	      }
		std::cout << "Successfully Loaded !  Normalization = " << Normalization << std::endl;
	      
	      
}



G4double DPL_1pDecay::GetProtonEng(){

    G4double RecRandnum = G4UniformRand();
     // Search for Prob Interval ----

    //G4cout<<"RANDOM "<<RecRandnum<<G4endl;

    G4int ProbInt = 0;

    for(G4int k=0; k<NumberOfLines; k++)
      {
	 if(k == 0 && RecRandnum < ProbLimit[k] && ProbRaw[k] > 0.)
	   { ProbInt = k; }
	 else if( k > 0 && RecRandnum >= ProbLimit[k-1] && RecRandnum < ProbLimit[k] && ProbRaw[k] > 0.)
	   { ProbInt = k; }
	 //else
	// ProbInt = NumberOfLines;
      }

    G4double EnergyOut = DaughterExEng[ProbInt];

    //G4cout<<"Energy "<<ProbInt<<" "<<DaughterExEng[ProbInt]<<G4endl;

    	return EnergyOut;

}



G4double DPL_1pDecay::GetProtonEng2(G4double Energy){

  G4int index=0;
  
  //G4cout<<"HERE!!!!! "<<ExEnergy.size()<<G4endl;

  for(G4int h=0;h<ExEnergy.size();h++){

    //G4cout<<"HERE!!!!! "<<h<<" "<<ExEnergy.at(h)<<G4endl;

    if(Energy*1000 == ExEnergy.at(h))
      index=h;    
    else
      continue;
  }
  
  
  G4double RecRandnum = G4UniformRand();
  // Search for Prob Interval ----
  
  //G4cout<<"RANDOM "<<RecRandnum<<G4endl;
  
  G4int ProbInt = 0;
  
  for(G4int k=0;k<Energy_table[index].size();k++){ 
    
    if(k == 0 && RecRandnum <  Limits_table[index][k] && Probability_table[index][k] > 0.)
      { ProbInt = k; }
    else if( k > 0 && RecRandnum >= Limits_table[index][k-1] && RecRandnum < Limits_table[index][k] && Probability_table[index][k] > 0.)
      { ProbInt = k; }
  }

  G4double EnergyOut = Energy_table[index][ProbInt]/1000.;

  //G4cout<<"Energy "<<Energy<<" "<<ExEnergy[index]<<" "<<RecRandnum<<" "<<ProbInt<<" "<<Energy_table[index][ProbInt]<<" "<<EnergyOut<<G4endl;

    	return EnergyOut;

}





G4bool DPL_1pDecay::IsApplicable(const G4ParticleDefinition& particle)
{
// returns "true" if this model can be applied to
// the particle type.

  //G4cout << "IsApplicable loaded!" << G4endl;

  // All particles, other than nuclei, are rejected by default.
  G4String ParticleName = particle.GetParticleName();

  if (particle.GetParticleType() == "nucleus")
  {
    //G4cout << "IsApplicable() : Particle " << ParticleName
    //       << " valid for DPL_1pDecay Process!" << G4endl;
    return true;
  }
  else
    {
      G4cout << "IsApplicable() : DPL_1pDecay can not be applied to " << ParticleName << G4endl;
      G4cout << "This Particle Type is not a Nucleus! " << G4endl;
	  G4ExceptionDescription description;
	  description << "Wrong particle type!";
      G4Exception("DPL_1pDecay::IsApplicable","Applicable_001",EventMustBeAborted,description);
      return false;
    }
  return false;
}


G4double DPL_1pDecay::GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize,
                                      G4ForceCondition* condition)
{
  // G4cout << "DPL_1pDecay GetMeanFreePath called!" << G4endl;

  *condition=NotForced;
  //*condition=Forced;  // Forces Process to Occur.
  previousStepSize = 0.1*mm;

  // 6/16/2009 -
  // Fix for geant4.9.2.p01 so that you can use this decay model and use
  // ionIonisation to stop the recoil particles.

   const G4DynamicParticle* Fragment = aTrack.GetDynamicParticle();
   const G4ParticleDefinition* Particle = Fragment->GetDefinition();
   G4int theCurrentMass = Particle->GetBaryonNumber();
   G4int theCurrentCharge = static_cast<int>(Particle->GetPDGCharge());
   G4double ExEngParticle = ((const G4Ions*)(Particle))->GetExcitationEnergy();

   // 23Mg -> 22Na + p only if above proton separation threshold (7580.3 keV)
   // 3/15/2010 -> included possiblity to decay by gammas

   if(theCurrentMass != A_Parent || theCurrentCharge != Z_Parent || ExEngParticle <= 100.*keV)
     {return DBL_MAX;}

  return 1./DBL_MAX;   // Forces Process to occur immediately. Perhaps later change to something
                      // with a short-lifetime consideration?
}

G4VParticleChange *DPL_1pDecay::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep){

	  const G4DynamicParticle* Fragment = aTrack.GetDynamicParticle();
	  const G4ParticleDefinition* Particle = aTrack.GetDefinition();

	  if (!IsApplicable(*Particle))  // Check applicability
	   {
	     G4cerr<<"General_1pDecay::PostStepDoIt: Only to be used with invalid particle !"<<G4endl;
	     return 0;
	   }

	  P_Frag = Fragment->Get4Momentum();    // 4-Momentum for the Fragment
	  A_Frag = Particle->GetBaryonNumber(); // Returns int. mass number of Frag.
	  Z_Frag = static_cast<int>(Particle->GetPDGCharge());    // Returns Charge of Fragment

	  G4double ExEngParticle = ((const G4Ions*)(Particle))->GetExcitationEnergy();

	  G4bool GammaFlag = false;

	  //Proton Decay occurs only if Ex>Sp!
	  //If Ex>threshold read proton energies from file!

	  G4double Energy_threshold= Proton_Separation_Energy+Delta_Separation_Energy;

	  G4double Proton_emission= G4UniformRand();
	  G4double Gamma_pOverGamma_gamma=1;
	  //G4double Gamma_pOverGamma_gamma=0; 
	//for the 6390 level
	// if(ExEngParticle<6.280){
	// Gamma_pOverGamma_gamma=1e-12;
	//G4cout<<"ExEng "<<ExEngParticle<<" "<<Gamma_pOverGamma_gamma<<G4endl;
	//Gamma_pOverGamma_gamma=0;
	// }
	// if(ExEngParticle>6.280&&ExEngParticle<6.491){
	//     Gamma_pOverGamma_gamma=4e-4;
	// }
	// if(ExEngParticle>6.491&&ExEngParticle<7.050){
	//     Gamma_pOverGamma_gamma=1;
	// }
	// if(ExEngParticle>7.050){	
	//     Gamma_pOverGamma_gamma=1;
	// }
//	else{
//	    Gamma_pOverGamma_gamma=1;
//	}
	//Gamma_pOverGamma_gamma=1; //Emable one proton per decay
	//Proton_emission=1 ;//disable proton emission
	//G4cout<<Proton_emission<<" Exc. Energy "<<ExEngParticle<<" Gamma_pOverGamma_gamma-> "<<Gamma_pOverGamma_gamma<<G4endl;

	  if(ExEngParticle >= Energy_threshold &&Proton_emission <Gamma_pOverGamma_gamma){

		  //G4cout<<"ExEng "<<ExEngParticle<<" "<<Energy_threshold<<G4endl;

		  A_Rec= A_Frag -1;
		  Z_Rec= Z_Frag -1;
		  E_rel = GetProtonEng2(ExEngParticle);
		  //E_rel =ExEngParticle-Energy_threshold-Excitation_Energy;

		  //G4cout<<"E_rel "<<E_rel<<G4endl;

		  if(E_rel<0){
		       E_rel = ExEngParticle;
		       A_Rec = A_Frag;
		       Z_Rec = Z_Frag;
		       GammaFlag = true;
			   }

	  }
	  else {
	       // 23Al->23Mg+gamma (direct gamma only)
	       E_rel = ExEngParticle;
	       A_Rec = A_Frag;
	       Z_Rec = Z_Frag;
	       GammaFlag = true;
	   }


	  Mass_Fragment = Fragment->GetMass();  // Mass of Fragment in GeV (note: includes ex. eng. mass!)

	  G4double KinEng_Int = Fragment->GetKineticEnergy();
	  G4ThreeVector theMomentum = aTrack.GetMomentumDirection();
	  G4TouchableHandle theTouchable = aTrack.GetTouchableHandle();

	  G4StepPoint* thePreStepPoint = aStep.GetPreStepPoint();
	  G4ThreeVector thePosition = thePreStepPoint->GetPosition();

	  G4double ProperTime = aTrack.GetProperTime();
	  G4double GlobalTime = aTrack.GetGlobalTime();

	  if(KinEng_Int < 1e-9*MeV)
	    {
	      //G4cout << "Kin Energy = " << KinEng_Int << G4endl;
	      //G4cout << "Do Nothing, Kill Event! " << G4endl;
	     aParticleChange.ProposeEnergy(0.);
	     aParticleChange.ProposeTrackStatus(fStopButAlive);
	     aParticleChange.ProposeGlobalTime(GlobalTime);
	     aParticleChange.ProposeProperTime(ProperTime);
	     return pParticleChange;
	    }

	  // Generate a Secondary Charged Recoil Particle
	      G4DynamicParticle* theRec = new G4DynamicParticle;
	      // GetIon() Method Arguments are GetIon(Charge,Mass,ExcitationEng)
	      G4ParticleDefinition* theRecDefinition = G4ParticleTable::GetParticleTable()->GetIon(Z_Rec,A_Rec,Excitation_Energy);
	      theRec->SetDefinition(theRecDefinition);
	      Mass_Recoil = theRec->GetMass();
	      //G4cout << "The Recoil was a " << theRecDefinition->GetParticleName() << G4endl;
	      //G4cout << "The Recoil Mass is : " << Mass_Recoil << " MeV " << G4endl;

	      // Generate a Secondary Proton
	          G4DynamicParticle* theProton = new G4DynamicParticle;
	          G4ParticleDefinition* theProtonDefinition;
	          if(GammaFlag == true)
	            {theProtonDefinition = G4Gamma::Gamma();}
	          else
	            {theProtonDefinition = G4Proton::Proton();}
	          theProton->SetDefinition(theProtonDefinition);
	          Mass_Proton = theProton->GetMass();
	          //G4cout << "The Proton mass = " << Mass_Proton << G4endl;

	          // Event Generated as follows. First generate 2-Body phase space decay of Rec + proton
	          // in BacktoBack(), then Boost from Center of Mass Frame to lab Frame according to the
	          // beam energy.

	              // if(RanCond == "Random")            // For acceptance tests, else it is constant
	              //  { E_rel = (20.*G4UniformRand())*MeV; }
	              //G4double deltaM = Mass_Fragment + E_rel - Mass_Recoil - Mass_Proton;

	              //G4cout << "E_rel = " << E_rel << G4endl;
	              G4double deltaM = E_rel;

	              // Generate Rec + p decay in center of mass frame
	              // Now use class variables to carry momentum between methods - 4 Oct. 2007.
	              BacktoBack(deltaM);

	              // Add Momentum Kick from Goldhaber if GoldCond = true
	              if(GoldCond == true)
	                {GoldhaberKick();}

	              // Boost to Lab Frame
	              CMFrameToLabFrame();

	          // Set Secondary information!

	              theRec->Set4Momentum(P_Recoil);
	              theProton->Set4Momentum(P_Proton);

	           // Kill the Unbound nucleus track
	             aParticleChange.ProposeTrackStatus(fStopAndKill); //Kills Parent Unbound Nucleus
	             aParticleChange.ProposeEnergy(0.);
	             aParticleChange.ProposePosition(thePosition);

	           // Set final Tracks in Motion!
	              G4Track* theRecTrack = new G4Track(theRec,GlobalTime,thePosition);
	              theRecTrack->SetTouchableHandle(theTouchable);

	              G4Track* theProtonTrack = new G4Track(theProton,GlobalTime,thePosition);
	              theProtonTrack->SetTouchableHandle(theTouchable);

	              aParticleChange.SetNumberOfSecondaries(2);
	              aParticleChange.AddSecondary(theRecTrack);
	              aParticleChange.AddSecondary(theProtonTrack);

	              //  G4cout << "Made it to the end ! " << G4endl;

	              return pParticleChange;



}
