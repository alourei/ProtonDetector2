/*
 * ProtonDetectorPrimaryGeneratorAction.cc
 *
 *  Created on: Dec 4, 2013
 *      Author: perezlou
 */

#include <ProtonDetectorPrimaryGeneratorAction.hh>

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

#include "G4GeneralParticleSource.hh"
#include "G4SingleParticleSource.hh"
#include "G4SPSEneDistribution.hh"
#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"

#include "Randomize.hh"
#include "G4UnitsTable.hh"

#include "ProtonDetectorROOTAnalysis.hh"

ProtonDetectorPrimaryGeneratorAction::ProtonDetectorPrimaryGeneratorAction(G4int sources, G4String particleName, G4String type, G4double energy = 5.0*MeV):
numberOfSources(sources), sourceType(type)
{
	G4cout<< "##################################################################" <<G4endl
			<<"######  ProtonDetectorPrimaryGenerationAction::ProtonDetectorPrimaryGeneratorAction()  ########" <<G4endl;
	G4cout<< "##################################################################"<<G4endl;


	//MargotDataOutPG = MargotDataRecordTree::MargotPointer;
	G4int nParticle = 1;
	particleGun = new G4ParticleGun(nParticle);
	theArchitect = new G4GeneralParticleSource();
	theArchitect->SetNumberOfParticles(nParticle);
	theArchitect->SetMultipleVertex(true);
	for(G4int i = 0; i<numberOfSources -1; i++){
		theArchitect->AddaSource(1.0);
	}
	theArchitect->ListSource();
	beamName = particleName;
	particleEnergy = energy;
	dataEvent=0;
	GPS_ParticleEnergy = NULL;
	GPS_ParticlePosition = NULL;
	GPS_ParticleMomentum = NULL;

	rootFile=new TFile("Zpos_al23.root");
	theHistogram=(TH1F*)rootFile->Get("h_pos");

}

ProtonDetectorPrimaryGeneratorAction::~ProtonDetectorPrimaryGeneratorAction() {
	delete theArchitect;
	delete particleGun;

 delete [] Z_beam;
 delete [] A_beam;
 delete [] Energy_beam;
 delete [] FWHM_Energy_beam;
 delete [] x_mean_beam;
 delete [] FWHM_x_beam;
 delete [] y_mean_beam;
 delete [] FWHM_y_beam;

 delete [] intensity_beam;
 delete [] probability_beam;
 delete [] probability_limits;


}

void ProtonDetectorPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent){

	//G4cout<< "##################################################################" <<G4endl
	//		<<"######  ProtonDetectorPrimaryGenerationAction::GeneratePrimaries()  ########" <<G4endl;
	//G4cout<< "##################################################################"<<G4endl;


	dataEvent++;
	G4int Z,A;

	//Define the Beam Particle depending on particleName
	//TODO Change this for a messenger in future
	if(beamName == "Al23"){
		Z=13;
		A=23;
	}
	else if(beamName == "Mg20"){
		Z=12;
		A=20;
	}
	else if(beamName == "Cl31"){
		Z=17;
		A=31;
	}
	else if(beamName == "Si25"){
		Z=14;
		A=25;
	}

	else if(beamName == "Ar32"){
		Z=18;
		A=32;
	}
	else if(beamName == "proton"){
		Z=1;
		A=1;
	}

	//Creates a Beam Ion with no Excitation Energy

	if(beamName != "geantino" && ( sourceType != "GPS"|| sourceType != "General")){
		theArchitect->SetParticleDefinition(G4ParticleTable::GetParticleTable()->GetIon(Z,A,0.));
		particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->GetIon(Z,A,0.));
	}
	else{
	  ;
	  //G4ParticleTable *theParticleTable = G4ParticleTable::GetParticleTable();
		//theArchitect->SetParticleDefinition(theParticleTable->FindParticle("geantino"));
		//particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->GetIon(Z,A,0.));
		//theArchitect->SetParticleDefinition(G4ParticleTable::GetParticleTable()->GetIon(Z,A,0.));
	}


	  G4double momentum_X = 0;   // Default momentum is z-direction for
	  G4double momentum_Y = 0;   // Pencil beams
	  G4double momentum_Z = 1;

	  G4double theta = 0.;
	  G4double phi= 0.;

	  const G4double Pi = CLHEP::pi;
	  const G4double twopi = 2*Pi;

	  G4ThreeVector theStartPosition(0.,0.,-30.5*cm);


	if(sourceType != "GPS" && sourceType !="General"){

	  G4cout<<"Here I am"<<G4endl;

	theSource[0] = theArchitect->GetCurrentSource();


	GPS_ParticleEnergy = theSource[0]->GetEneDist();
	GPS_ParticleEnergy->SetEnergyDisType("Mono");
	GPS_ParticleEnergy->SetMonoEnergy(particleEnergy);

	GPS_ParticleMomentum = theSource[0]->GetAngDist();

	GPS_ParticlePosition = theSource[0]->GetPosDist();
	GPS_ParticlePosition ->SetPosDisType("Point");

	//****************************************************
	// Select type of beam
	//****************************************************
	// SourceType = "test" -> Beam on z-axis at beam energy given in the main()
	// SourceType = "pencil" -> Pencil beam with +- 5% energy variance
	// SourceType = "RealBeam" -> 5mm Gaussian in x-y beam spot with 1 deg. divergence


	  GPS_ParticlePosition->SetCentreCoords(theStartPosition);

	}
	  if(sourceType == "test")
	  {
	     GPS_ParticlePosition->SetCentreCoords(theStartPosition);
	     momentum_X = 0.;
	     momentum_Y = 0.;
	     momentum_Z = 1.;
	     //MargotDataOutPG->senddataPG(particleEnergy,theStartPosition);
	     GPS_ParticleMomentum = theSource[0]->GetAngDist();
	     G4ThreeVector v(momentum_X,momentum_Y,momentum_Z);
	     GPS_ParticleMomentum->SetParticleMomentumDirection(v);
	     theArchitect->GeneratePrimaryVertex(anEvent);
	    }

	  else if(sourceType == "pencil"){

		  //Straight Pencil Beam
		  momentum_X = 0.;
		  momentum_Y = 0.;
		  momentum_Z = 1.;
		  GPS_ParticlePosition->SetCentreCoords(theStartPosition);

		  // Now create Beam Energy Spread of +- 0.5 % !
		  //G4double energySpread = 0.01*particleEnergy;
		  G4double energySpread = 0.0*particleEnergy;
		  G4double energyLowLimit = particleEnergy-0.5*energySpread;
		  G4double deltaEnergy = energySpread*G4UniformRand();
		  G4double energyOut= energyLowLimit+deltaEnergy;
		  GPS_ParticleEnergy->SetMonoEnergy(energyOut);
		  //GPS_ParticleEnergy->SetMonoEnergy(energyOut);
		  G4ThreeVector v(momentum_X,momentum_Y,momentum_Z);
		  GPS_ParticleMomentum->SetParticleMomentumDirection(v);
		  theArchitect->GeneratePrimaryVertex(anEvent);
		  //MargotDataOutPG->senddataPG(energyOut,theStartPosition);
	  }
	  else if(sourceType == "RealBeam"){

	      // +- .5% energy Spread
	      // +- 2.9% energy Spread LISE++
		  G4double energySpread = 0.01*particleEnergy;
		  G4double energyLowLimit = particleEnergy-2.9*energySpread;
		  G4double deltaEnergy = energySpread*G4UniformRand();
		  G4double energyOut= energyLowLimit+deltaEnergy;
		  //GPS_ParticleEnergy->SetMonoEnergy(energyOut);
		  particleGun->SetParticleEnergy(energyOut);


	      // 5 mm gaussian random starting position in x-y
		  // Positions from LISE++ calculation
		  G4double FWHM_X = 22.6*mm;
		  G4double FWHM_Y = 4.804*mm;
		  G4double mean_X = -0.1005*mm;
		  //G4double mean_Y = 6.355*mm;
		  G4double mean_Y = 0*mm;
		  G4double conversionFactor = 2.35482;   // 2*sqrt(2*log(2)) - converts FWHM to Sigma
		  G4double sigmaX =  FWHM_X/conversionFactor;
		  G4double sigmaY =  FWHM_Y/conversionFactor;
		  theStartPosition[0]=CLHEP::RandGauss::shoot(mean_X,sigmaX)*mm;
		  theStartPosition[1]=CLHEP::RandGauss::shoot(mean_Y,sigmaY)*mm;
		  theStartPosition[2]=-27*cm;
		  //GPS_ParticlePosition->SetCentreCoords(theStartPosition);
		  particleGun->SetParticlePosition(theStartPosition);

	      // 1 deg. divergence (from starting point)
		  G4double beamDivergence = 1.*deg;
		  G4double ran = G4UniformRand();
		  // isotropic in cosine for OpenAngle (in degrees)
		  theta = acos(1+(cos(beamDivergence)-1)*ran);
		  ran = G4UniformRand();
		  phi = twopi*ran;

		  momentum_X = sin(theta)*cos(phi);
		  momentum_Y = sin(theta)*sin(phi);
		  momentum_Z = cos(theta);

		  G4ThreeVector v(momentum_X ,momentum_Y ,momentum_Z);
		  //GPS_ParticleMomentum->SetParticleMomentumDirection(v);
		  particleGun->SetParticleMomentumDirection(v);
		  //MargotDataOutPG->senddataPG(energyOut,theStartPosition);
		  particleGun->GeneratePrimaryVertex(anEvent);
	  }
	  else if(sourceType == "General"){

		  for(G4int i=0; i<numberOfSources;i++){
			  theArchitect->SetCurrentSourceto(i);
			  theSource[i]=theArchitect->GetCurrentSource();
			  theSource[i]->SetParticleDefinition(G4ParticleTable::GetParticleTable()->GetIon(Z,A,0.));
			  //particleEnergy = 0.4*eV;
			  //GPS_ParticleEnergy->SetEnergyDisType("Mono");
			  //GPS_ParticleEnergy->SetMonoEnergy(particleEnergy);
			  // 9 mm gaussian random starting position in x-y
			  G4double FWHM_X = 26.33*mm;
			  G4double FWHM_Y = 5.178*mm;

			  G4double mean_X = -2.513*mm;
			  G4double mean_Y =  0.02015*mm;

			  G4double conversionFactor = 2.35482;   // 2*sqrt(2*log(2)) - converts FWHM to Sigma
			  G4double sigmaX =  FWHM_X/conversionFactor;
			  G4double sigmaY =  FWHM_Y/conversionFactor;
			  G4double centroid = 0.*mm;
			  theStartPosition[0]=CLHEP::RandGauss::shoot(mean_X,sigmaX)*mm;
			  theStartPosition[1]=CLHEP::RandGauss::shoot(mean_Y,sigmaY)*mm;
			  // 24 mm gaussian random starting position in z 800 torr
			  G4double FWHMz = 0*mm;
			  G4double DeltaZ =  FWHMz*(1-2*G4UniformRand());
			  G4double sigmaz = FWHMz/conversionFactor;
			  G4double centroidz = -1.01*mm;
			  //theStartPosition[2]=CLHEP::RandGauss::shoot(centroidz,sigmaz)*mm;
			  theStartPosition[2]=theHistogram->GetRandom()*mm;

			  GPS_ParticlePosition = theSource[i]->GetPosDist();
			  GPS_ParticlePosition->SetPosDisType("Point");
			  GPS_ParticlePosition->SetCentreCoords(theStartPosition);


			  if(beamName == "proton"){
				  // Use isotropic momentum dist. for proton eff. checks
				  G4double OpenAngle = 180.*deg;
				  G4double ran = G4UniformRand();
				  // isotropic in cosine for OpenAngle (in degrees)
				  theta = acos(1+(cos(OpenAngle)-1)*ran);
				  ran = G4UniformRand();
				  phi = twopi*ran;
			  }
			  momentum_X = sin(theta)*cos(phi);
			  momentum_Y = sin(theta)*sin(phi);
			  momentum_Z = cos(theta);

			  G4ThreeVector v(momentum_X,momentum_Y,momentum_Z);

			  //MargotDataOutPG->senddataPG(particleEnergy,theStartPosition);

			  GPS_ParticleMomentum = theSource[i]->GetAngDist();
			  GPS_ParticleMomentum->SetParticleMomentumDirection(v);
		  }
		  theArchitect->GeneratePrimaryVertex(anEvent);

		  }
	  else if (sourceType == "FragmentBeam"){


		  G4double RecRandnum = G4UniformRand();

		  G4int ProbInt = 0;

		    for(G4int k=0; k<numberOfFragments; k++)
		      {
			 if(k == 0 && RecRandnum < probability_limits[k] && probability_beam[k] > 0.)
			   { ProbInt = k; }
			 else if( k > 0 && RecRandnum >= probability_limits[k-1] && RecRandnum < probability_limits[k] && probability_beam[k] > 0.)
			   { ProbInt = k; }
		      }

			 Z=Z_beam[ProbInt];
			 A=A_beam[ProbInt];



		  particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->GetIon(Z,A,0.));


		  G4double totalEnergy=Energy_beam[ProbInt]*A_beam[ProbInt];

		  // +- .5% energy Spread
	      // +- 2.9% energy Spread LISE++
		  G4double energySpread = FWHM_Energy_beam[ProbInt]/2.*A_beam[ProbInt];;
		  //G4double energyLowLimit = Energy_beam[ProbInt]-2.9*energySpread;
		  G4double deltaEnergy = energySpread*(1-2*G4UniformRand());
		  G4double energyOut= totalEnergy+deltaEnergy;
		  //GPS_ParticleEnergy->SetMonoEnergy(energyOut);
		  particleGun->SetParticleEnergy(energyOut);


	      // 5 mm gaussian random starting position in x-y
		  // Positions from LISE++ calculation
		  G4double FWHM_X = FWHM_x_beam[ProbInt]*mm;
		  G4double FWHM_Y = FWHM_y_beam[ProbInt]*mm;
		  G4double mean_X = x_mean_beam[ProbInt]*mm;
		  G4double mean_Y = y_mean_beam[ProbInt]*mm;
		  //G4double mean_Y = 0*mm;
		  G4double conversionFactor = 2.35482;   // 2*sqrt(2*log(2)) - converts FWHM to Sigma
		  G4double sigmaX =  FWHM_X/conversionFactor;
		  G4double sigmaY =  FWHM_Y/conversionFactor;
		  theStartPosition[0]=CLHEP::RandGauss::shoot(mean_X,sigmaX)*mm;
		  theStartPosition[1]=CLHEP::RandGauss::shoot(mean_Y,sigmaY)*mm;
		  theStartPosition[2]=-27*cm;
		  //GPS_ParticlePosition->SetCentreCoords(theStartPosition);
		  particleGun->SetParticlePosition(theStartPosition);

	      // 1 deg. divergence (from starting point)
		  G4double beamDivergence = 1.*deg;
		  G4double ran = G4UniformRand();
		  // isotropic in cosine for OpenAngle (in degrees)
		  theta = acos(1+(cos(beamDivergence)-1)*ran);
		  ran = G4UniformRand();
		  phi = twopi*ran;

		  momentum_X = sin(theta)*cos(phi);
		  momentum_Y = sin(theta)*sin(phi);
		  momentum_Z = cos(theta);

		  G4ThreeVector v(momentum_X ,momentum_Y ,momentum_Z);
		  //GPS_ParticleMomentum->SetParticleMomentumDirection(v);
		  particleGun->SetParticleMomentumDirection(v);
		  //MargotDataOutPG->senddataPG(energyOut,theStartPosition);
		  particleGun->GeneratePrimaryVertex(anEvent);

	  }
	  else if(sourceType == "GPS"){

	    // G4cout<<"GPS source"<<G4endl;

	    theArchitect->GeneratePrimaryVertex(anEvent);

	  }

	  if(gProtonDetectorROOTAnalysis){
		//G4cout<<"HERE!!"<<G4endl;
		//G4PrimaryVertex *theVertex=anEvent->GetPrimaryVertex();
		//G4cout<<"theVertex "<<theVertex<<G4endl;
		gProtonDetectorROOTAnalysis->GenerateBeam(anEvent);
	  }
}

void ProtonDetectorPrimaryGeneratorAction::SetBeamData(G4String theFileName){

	G4cout<<"Getting Energies for fragments"<<G4endl;

	system("echo $G4BETADECAYSIDETAL23");

	if(!getenv("G4BETADECAYSIDETAL23"))
	     {
	       //G4cerr << "Please setenv G4BETADECAYSIDETAL23 to point to the Z="<<Z_Parent<<"A="<<A_Parent<<"branching data file" << G4endl;
	 	  G4ExceptionDescription description;
	 	  description << "Environment variable G4BETADECAYSIDETAL23 not set! ";
	       G4Exception("ProtonDetectorPrimaryGeneratorAction::SetBeamata()","BeamData_001",FatalException,description);
	     }
	   G4String  DirName = getenv("G4BETADECAYSIDETAL23");

	   G4String FileName = DirName+"/"+theFileName;
	   std::fstream theFile;
	   theFile.open(FileName, std::fstream::in );
	   theFile >> numberOfFragments;

	   Z_beam = new G4int[numberOfFragments];
	   A_beam = new G4int[numberOfFragments];
	   Energy_beam = new G4double[numberOfFragments];
	   FWHM_Energy_beam = new G4double[numberOfFragments];
	   x_mean_beam = new G4double[numberOfFragments];
	   FWHM_x_beam = new G4double[numberOfFragments];
	   y_mean_beam = new G4double[numberOfFragments];
	   FWHM_y_beam = new G4double[numberOfFragments];

	   intensity_beam = new G4double[numberOfFragments];
	   probability_beam = new G4double[numberOfFragments];
	   probability_limits = new G4double[numberOfFragments];

	   normalization = 0.;

	   if(theFile.good()){

		   G4cout << "Loading Data For FileName = " << FileName << G4endl;
	        for(G4int i=0; i<numberOfFragments; i++)
	          {
	        	G4String theName;
	        	G4int theCharge=0;
	        	G4int theMass=0;
	        	G4double theEnergy = 0.;
	        	G4double theEnergy_Res = 0.;
	        	G4double thePositionX =0.;
	        	G4double thePositionX_Res=0;
	        	G4double thePositionY =0.;
	        	G4double thePositionY_Res=0;
	        	G4double theIntensity = 0.;

	   	 theFile >> theName >>theCharge >> theMass >> theEnergy >> theEnergy_Res >> thePositionX >> thePositionX_Res
	   	>> thePositionY >> thePositionY_Res >> theIntensity;

	   	 G4cout << theName << "  " << theIntensity << G4endl;


		   Z_beam[i]= theCharge;
		   A_beam[i]= theMass;
		   Energy_beam[i] = theEnergy;
		   FWHM_Energy_beam[i] = theEnergy_Res;
		   x_mean_beam[i] = thePositionX;
		   FWHM_x_beam[i] = thePositionX_Res;
		   y_mean_beam[i] = thePositionY;
		   FWHM_y_beam[i] = thePositionX_Res;
		   intensity_beam[i] = theIntensity;

		   normalization += theIntensity;
	          }
	        G4cout << "Successfully Loaded !  Normalization = " << normalization << G4endl;
	        theFile.close();

	        for(G4int i=0; i<numberOfFragments; i++)
	          {
	        	probability_beam[i] = intensity_beam[i]/normalization;
	          }

	   }
	   else
	     {
	       G4cerr << "File = " << FileName << " not found or in improper format." << G4endl;
	    	 G4ExceptionDescription description;
	  	 description<< "File not found!!";
	       G4Exception("ProtonDetectorPrimaryGeneratorAction::SetBeamData()","BeamData_002",FatalException,description);
	   }

	      // Setup Limits and Distances

	      for(G4int i=0; i<numberOfFragments; i++)
		{
		  if(i == 0)
		    {probability_limits[i] = probability_beam[i];}
		  else
		    {probability_limits[i] = probability_limits[i-1]+probability_beam[i];}
		}


}


