/*
 * ProtonDetectorROOTAnalysis.cc
 *
 *  Created on: Dec 10, 2013
 *      Author: perezlou
 */

#include <ProtonDetectorROOTAnalysis.hh>

#include<time.h>
//Geant4 includes
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4Track.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "G4TrackStatus.hh"
#include "G4Step.hh"
#include "G4Types.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Trajectory.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"

//another includes
#include <G4PrimaryVertex.hh>
#include <G4TrajectoryContainer.hh>
//ROOT includes
#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TProfile.h"

//ProtonDetector includes
#include <ProtonDetectorBeamInfo.hh>
#include <ProtonDetectorSimpleTrack.hh>
#include <ProtonDetectorTrack.hh>
#include <ProtonDetectorData.hh>

//global pointer to the ROOT analysis manager
ProtonDetectorROOTAnalysis *gProtonDetectorROOTAnalysis = (ProtonDetectorROOTAnalysis *)0;


ProtonDetectorROOTAnalysis::ProtonDetectorROOTAnalysis() {
	 LastDoItTime = (time_t)0;
	  if(gSystem) gSystem->ProcessEvents();

	  if(gProtonDetectorROOTAnalysis)
	    delete gProtonDetectorROOTAnalysis;
	  gProtonDetectorROOTAnalysis = this;

	  newDirName = new char[255];

	  innerPadFlag = "on"; //On for inner Pad
	  //innerPadFlag = "off";

	  storeTracksFlag = "off";
	  //storeSimpleTracksFlag = "on";
	  storeSimpleTracksFlag = "off";
	  simFile = 0;
	  eventTree=0;
	  tracksTree=0;


	InitAnalysis();

	OnceAWhileDoIt();


}

ProtonDetectorROOTAnalysis::~ProtonDetectorROOTAnalysis() {

	 //delete analMessenger;

	  simFile->Write();
	  simFile->Close();

	  if (gProtonDetectorROOTAnalysis == this)
	    gProtonDetectorROOTAnalysis = (ProtonDetectorROOTAnalysis *)0;


	  if (gSystem) gSystem->ProcessEvents();

}

void ProtonDetectorROOTAnalysis::InitAnalysis(){

	  G4cout << "##################################################################" << G4endl
		 << "###########    ProtonDetectorROOTAnalysis::InitAnalysis()    ####################" << G4endl;
	  G4cout << "################################################################## " << G4endl;

	  minStrideLength=0.01*mm;
	  //minStrideLength=0.0*mm;

	  //TFile for storing the info
	  if(!simFile){
	    simFile = new TFile("root_files/simFile.root","RECREATE");
	    simFile->cd();
	    eventTree = new TTree("The_Event_Tree","Event_Tree");
	    tracksTree = new TTree("The_Tracks_Tree","Tracks_Tree");

	  }

	  //The clones array of Crystal Hits

		  simpleTrackCA = new TClonesArray("ProtonDetectorSimpleTrack",100);
		  eventTree->Branch("simpleTrackData",&simpleTrackCA);




	  theTracks = new ProtonDetectorTrack();
	  tracksTree->Branch("trackData","ProtonDetectorTrack",&theTracks,128000,99);

	   //Now, simple track as a TClonesArray


	   theBeamInfo=new ProtonDetectorBeamInfo();

	   eventTree->Branch("beamInfo",&theBeamInfo,128000,99);

	   theData=new ProtonDetectorData();


	   eventTree->Branch("detectorData",&theData,128000,99);

	   //minStrideLength = 0.5 * mm; //default value for the minimum stride length

	   //simFile->cd();

	   OnceAWhileDoIt(true); // do it now


}

void ProtonDetectorROOTAnalysis::Construct(const G4VPhysicalVolume *theWorldVolume) {
  //
  // Things to do while contructing...
  //

  if (theWorldVolume) {;} /* keep the compiler "quiet" */
  if (gSystem) gSystem->ProcessEvents();

  OnceAWhileDoIt();
}

void ProtonDetectorROOTAnalysis::ConstructParticle(){
  //
  // Actions to perform in the analysis during the particle construction
  //

  if (gSystem) gSystem->ProcessEvents();

  OnceAWhileDoIt();
}


void ProtonDetectorROOTAnalysis::ConstructProcess(){
  //
  // Actions to perform in the analysis during the processes construction
  //

  if (gSystem) gSystem->ProcessEvents();

  OnceAWhileDoIt();
}


void ProtonDetectorROOTAnalysis::SetCuts() {
  //
  // Actions to perform in the analysis during the cut setting
  //

  if (gSystem) gSystem->ProcessEvents();

  OnceAWhileDoIt();
}


void ProtonDetectorROOTAnalysis::GenerateBeam(const G4Event *anEvent){
//
// Defining any beam related histogram or information in the output file
//

	//G4cout<<"ProtonDetectorROOTAnalysis::GenerateBeam()"<<G4endl;

	if (gSystem) gSystem->ProcessEvents();

	G4PrimaryVertex *theVertex= anEvent->GetPrimaryVertex();

	theBeamInfo->SetXInitial(theVertex->GetX0());
	theBeamInfo->SetYInitial(theVertex->GetY0());
	theBeamInfo->SetZInitial(theVertex->GetZ0());

	//G4cout<<"Initial Position of the Beam "<<theBeamInfo->GetZInitial()/mm<<G4endl;
	//G4cout<<"Initial Position of the Beam "<<theVertex->GetZ0()/mm<<G4endl;


	G4PrimaryParticle *theParticle = theVertex->GetPrimary();

	//theBeamInfo->SetCharge(theParticle->GetCharge());
	//theBeamInfo->SetMass(theParticle->GetMass());

	G4double theEnergy= theParticle->GetKineticEnergy()/MeV;

	theBeamInfo->SetEnergyEntrance(theEnergy);

	G4ThreeVector theMomentum=theParticle->GetMomentum();

	G4double thetaEntrance = theMomentum.theta()/deg;
	G4double phiEntrance = theMomentum.phi()/deg;

	theBeamInfo->SetThetaEntrance(thetaEntrance);
	theBeamInfo->SetPhiEntrance(phiEntrance);


	OnceAWhileDoIt();

}



void ProtonDetectorROOTAnalysis::BeginOfRunAction(const G4Run *aRun){
  //
  // Actions to perform in the analysis at the beginning of the run
  //
  if (gSystem) gSystem->ProcessEvents();
  //Storing the runID
  SetTheRunID(aRun->GetRunID());

  //going to the file!!!
  G4cout << "##################################################################" << G4endl
	 << "########  ProtonDetectorROOTAnalysis::BeginOfRunAction()    ############" << G4endl;
  G4cout << "########  New Run With Number " << aRun->GetRunID() << " Detected!!  ######" << G4endl;
  G4cout << "####  A new directory will be opened in the output ROOT file  ####" << G4endl;
  G4cout << "################################################################## " << G4endl;

  sprintf(newDirName,"Run_%02d",aRun->GetRunID());
  simFile->mkdir(newDirName,newDirName);
  simFile->cd(newDirName);

  //simFile->cd();

   OnceAWhileDoIt(true); // do it now

}

void ProtonDetectorROOTAnalysis::EndOfRunAction(const G4Run *aRun){
	 if (aRun) {;} /* keep the compiler "quiet" */
	  if (gSystem) gSystem->ProcessEvents();

	  G4cout << "##################################################################" << G4endl
		 << "########  ProtonDetectorROOTAnalysis::EndOfRunAction()    ############" << G4endl;
	  G4cout << "########  Run With Number " << aRun->GetRunID() << " finished!  ######" << G4endl;
	  G4cout << "################################################################## " << G4endl;

	  OnceAWhileDoIt(true); // do it now


}

void ProtonDetectorROOTAnalysis::BeginOfEventAction(const G4Event *anEvent){

	SetTheEventID(anEvent->GetEventID());
}

void ProtonDetectorROOTAnalysis::EndOfEventAction(const G4Event *anEvent){

	G4int hitsCollectionID =
	G4SDManager::GetSDMpointer()->GetCollectionID("gasCollection");

	G4int degCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("degCollection");

	G4int vetoCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("VetoCollection");

	G4HCofThisEvent* HCofEvent = anEvent->GetHCofThisEvent();

	G4TrajectoryContainer *theTrajContainer=anEvent->GetTrajectoryContainer();
	G4int NbTrajectories;	      
	if(theTrajContainer)
			NbTrajectories = theTrajContainer->entries();
	else		
	  NbTrajectories =0;
	//theVertex->Print();

	//G4cout<< "Nb of trajectories "<<NbTrajectories<<G4endl;

	ProtonDetectorGeantHitsCollection *hitsCollection =
			(ProtonDetectorGeantHitsCollection*)HCofEvent->GetHC(hitsCollectionID);

	ProtonDetectorGeantHitsCollection *degHitsCollection =
			(ProtonDetectorGeantHitsCollection*)HCofEvent->GetHC(degCollectionID);

	ProtonDetectorGeantHitsCollection *vetoHitsCollection =
			(ProtonDetectorGeantHitsCollection*)HCofEvent->GetHC(vetoCollectionID);

	//Number of ActarSimGasGeantHit (or steps) in the hitsCollection

	G4int NbHits = hitsCollection->entries();

	G4int NbHitsDeg = degHitsCollection->entries();
	G4int NbHitsVeto = vetoHitsCollection->entries();


	  //Let us create as many simpleTracks as the primaries...
	   simpleTrack = new ProtonDetectorSimpleTrack*[NbTrajectories];
	   for (G4int i=0;i<NbTrajectories;i++)
	     simpleTrack[i] = new ProtonDetectorSimpleTrack();



	G4int NbStrides = 0;

	simpleTrackCA->Clear();

	//G4cout << "Information on the collection..." << G4endl
	//<< "Number of GasGeantHits in the collection: " <<  NbHits
	//<< G4endl;

    G4int *strideOrdinal= new G4int[NbTrajectories];

	Double_t energyOnGas_beam=0;
	Double_t energyOnGas_protons=0;
	Double_t energyOnGas_betas=0;
	Double_t energyOnGas_alphas=0;
	Double_t energyOnGas_total=0;
	Double_t stepSumLengthOnGas_beam=0;
	Double_t stepSumLengthOnGas_protons=0;
	Double_t stepSumLengthOnGas_betas=0;
	Double_t stepSumLengthOnGas_alphas=0;
	Double_t beam_last_positionX=-999;
	Double_t beam_last_positionY=-999;
	Double_t beam_last_positionZ=-999;
	Double_t proton_last_positionX=-999;
	Double_t proton_last_positionY=-999;
	Double_t proton_last_positionZ=-999;
	Double_t alpha_last_positionX=-999;
	Double_t alpha_last_positionY=-999;
	Double_t alpha_last_positionZ=-999;

	Double_t energyOnDegrader=0;
	Double_t energyOnVeto=0;
	Double_t energyOnVeto_protons=0;
	Double_t energyOnVeto_betas=0;

	Int_t proton_multiplicity=0;
	Int_t beta_multiplicity=0;
	Int_t total_multiplicity=0;

	Int_t proton_Vetomultiplicity=0;
	Int_t beta_Vetomultiplicity=0;
	Int_t total_Vetomultiplicity=0;

	Int_t proton_padNumber=-1;
	Int_t beta_padNumber=-1;

	Double_t Energy_pad_protons[NPADS]={0};
	Double_t Energy_pad_betas[NPADS]={0};
	Double_t Energy_pad_total[NPADS]={0};

	Double_t Energy_veto_protons[NVETO]={0};
	Double_t Energy_veto_betas[NVETO]={0};
	Double_t Energy_veto_total[NVETO]={0};

	G4bool FirstBetaDetectorHit=true;
	G4bool FirstProtonDetectorHit=true;
	G4bool FirsRecoilDetectorHit=true;

	G4bool FirstBetaVetoHit=true;
	G4bool FirstProtonVetoHit=true;
	G4bool FirsRecoilVetoHit=true;
	
	G4double BetaKineticEnergyDetector=0;
	G4double ProtonKineticEnergyDetector=0;
	G4double RecoilKineticEnergyDetector=0;
	
	G4double BetaKineticEnergyVeto=0;
	G4double ProtonKineticEnergyVeto=0;
	G4double RecoilKineticEnergyVeto=0;

	G4double BetaKineticEnergy=0;
	G4double ProtonKineticEnergy=0;
	G4double RecoilKineticEnergy=0;

	//Loop on degrader Hits
	   for(G4int i=0; i<NbHitsDeg;i++){
		   energyOnDegrader+=(*degHitsCollection)[i]->GetEdep();

	   }

	//Loop on Veto Hits
	   for(G4int i=0; i<NbHitsVeto;i++){

	   Int_t PadNumber=(*vetoHitsCollection)[i]->GetDetID();
	   
	     if((*vetoHitsCollection)[i]->GetParticleID()==2212){
	       energyOnVeto_protons+=(*vetoHitsCollection)[i]->GetEdep();
	       if(FirstProtonVetoHit){
		 ProtonKineticEnergyVeto=(*vetoHitsCollection)[i]->GetParticleKineticEnergy();
		 FirstProtonVetoHit=false;		  
	       }
	       Energy_veto_protons[PadNumber]+=(*vetoHitsCollection)[i]->GetEdep();
	       Energy_veto_total[PadNumber]+=(*vetoHitsCollection)[i]->GetEdep();
	       energyOnVeto+=(*vetoHitsCollection)[i]->GetEdep();
	     }
	     if((*vetoHitsCollection)[i]->GetParticleID()==-11){
	       Energy_veto_betas[PadNumber]+=(*vetoHitsCollection)[i]->GetEdep();
	       if(FirstBetaVetoHit){
		 BetaKineticEnergyVeto=(*vetoHitsCollection)[i]->GetParticleKineticEnergy();
		 FirstBetaVetoHit=false;		  
	       }

	       Energy_veto_total[PadNumber]+=(*vetoHitsCollection)[i]->GetEdep();
	       energyOnVeto_betas+=(*vetoHitsCollection)[i]->GetEdep();
	       energyOnVeto+=(*vetoHitsCollection)[i]->GetEdep();
	     }
	     
	     if((*vetoHitsCollection)[i]->GetParticleCharge()>2){
	       Energy_veto_total[PadNumber]+=(*vetoHitsCollection)[i]->GetEdep();
	       energyOnVeto+=(*vetoHitsCollection)[i]->GetEdep();
	       if(FirsRecoilVetoHit){
		 RecoilKineticEnergyVeto=(*vetoHitsCollection)[i]->GetParticleKineticEnergy();
		 FirsRecoilVetoHit=false;		  
	       }
	       
	     }
	   }

    for(G4int i=0; i<NbHits;i++){
      for(G4int j=0;j<NbTrajectories;j++){
	
	if(storeSimpleTracksFlag=="on"){
	  
	  if((*hitsCollection)[i]->GetTrackID()==(j+1) && (*hitsCollection)[i]->GetParticleCharge()!=0){
	  //if((*hitsCollection)[i]->GetTrackID()==(j+1)){

	  if(simpleTrack[j]->GetNumberSteps() == 0) {
	    //the first step in the stride!
	    simpleTrack[j]->SetXPre((*hitsCollection)[i]->GetPrePos().x());
	    simpleTrack[j]->SetYPre((*hitsCollection)[i]->GetPrePos().y());
	    simpleTrack[j]->SetZPre((*hitsCollection)[i]->GetPrePos().z());
	    simpleTrack[j]->SetXPost((*hitsCollection)[i]->GetPostPos().x());
	    simpleTrack[j]->SetYPost((*hitsCollection)[i]->GetPostPos().y());
	    simpleTrack[j]->SetZPost((*hitsCollection)[i]->GetPostPos().z());
	    simpleTrack[j]->SetEnergyStride  ((*hitsCollection)[i]->GetEdep());
	    simpleTrack[j]->SetParticleCharge((*hitsCollection)[i]->GetParticleCharge());
	    simpleTrack[j]->SetParticleMass((*hitsCollection)[i]->GetParticleMass());
	    simpleTrack[j]->SetParticleID((*hitsCollection)[i]->GetParticleID());
	    simpleTrack[j]->SetStrideLength((*hitsCollection)[i]->GetStepLength());
	    simpleTrack[j]->SetTimePre((*hitsCollection)[i]->GetPreToF());
	    simpleTrack[j]->SetTimePost((*hitsCollection)[i]->GetPostToF());
	    simpleTrack[j]->SetNumberSteps(1);
	    simpleTrack[j]->SetTrackID((*hitsCollection)[i]->GetTrackID());
	    simpleTrack[j]->SetParentTrackID((*hitsCollection)[i]->GetParentID());
	    simpleTrack[j]->SetEventID(GetTheEventID());
	    simpleTrack[j]->SetRunID(GetTheRunID());
	    simpleTrack[j]->SetStrideOrdinal(strideOrdinal[j]);
	    simpleTrack[j]->SetVolumeName((*hitsCollection)[i]->GetDetName());
			 
	    simpleTrack[j]->SetProcessName((*hitsCollection)[i]->GetProcessName());
	  }
	  else{
	    simpleTrack[j]->SetXPost((*hitsCollection)[i]->GetPostPos().x());
	    simpleTrack[j]->SetYPost((*hitsCollection)[i]->GetPostPos().y());
	    simpleTrack[j]->SetZPost((*hitsCollection)[i]->GetPostPos().z());
	    simpleTrack[j]->SetTimePost((*hitsCollection)[i]->GetPostToF());
	    simpleTrack[j]->SetEnergyStride(simpleTrack[j]->GetEnergyStride() +
					    (*hitsCollection)[i]->GetEdep());
	    simpleTrack[j]->SetStrideLength(simpleTrack[j]->GetStrideLength() +
					    (*hitsCollection)[i]->GetStepLength());
	    simpleTrack[j]->SetNumberSteps(simpleTrack[j]->GetNumberSteps()+1);
	  }
	  
	  if(simpleTrack[j]->GetStrideLength() > minStrideLength ){
	    
	    //Fill the Clones Array
	    new((*simpleTrackCA)[NbStrides])ProtonDetectorSimpleTrack(*simpleTrack[j]);
	    
	    NbStrides++;
	    strideOrdinal[j]++;
	    simpleTrack[j]->Reset();
	    
	  }
	  }
	
	}
      }


      
      
      if((*hitsCollection)[i]->GetParticleCharge()>2){
	energyOnGas_beam+=(*hitsCollection)[i]->GetEdep();
	//energyOnGas_total+=(*hitsCollection)[i]->GetEdep();
	if(FirsRecoilDetectorHit){
	  RecoilKineticEnergyDetector=(*hitsCollection)[i]->GetParticleKineticEnergy();
	  FirsRecoilDetectorHit=false;		  
	}
	
	stepSumLengthOnGas_beam+= (*hitsCollection)[i]->GetStepLength();
	beam_last_positionX=(*hitsCollection)[i]->GetPostPos().x();
	beam_last_positionY=(*hitsCollection)[i]->GetPostPos().y();
	beam_last_positionZ=(*hitsCollection)[i]->GetPostPos().z();
      }
      
      else if((*hitsCollection)[i]->GetParticleID()==2212){
	energyOnGas_protons+=(*hitsCollection)[i]->GetEdep();
	energyOnGas_total+=(*hitsCollection)[i]->GetEdep();
	if(FirstProtonDetectorHit){
	  ProtonKineticEnergyDetector=(*hitsCollection)[i]->GetParticleKineticEnergy();
	  FirstProtonDetectorHit=false;		  
	}

	proton_last_positionX=(*hitsCollection)[i]->GetPostPos().x();
	proton_last_positionY=(*hitsCollection)[i]->GetPostPos().y();
	proton_last_positionZ=(*hitsCollection)[i]->GetPostPos().z();
	stepSumLengthOnGas_protons+= (*hitsCollection)[i]->GetStepLength();
	G4String detName=(*hitsCollection)[i]->GetDetName();
	Int_t PadNumber=(*hitsCollection)[i]->GetDetID();

	if(innerPadFlag == "on"){
	  if(detName == "InnerPad" && PadNumber==0){
	    Energy_pad_protons[0]+=(*hitsCollection)[i]->GetEdep(); 
	    Energy_pad_total[0]+=(*hitsCollection)[i]->GetEdep();
	  }
	  else if (detName == "OuterPad"){
	    Energy_pad_protons[PadNumber+1]+=(*hitsCollection)[i]->GetEdep(); 
	    Energy_pad_total[PadNumber+1]+=(*hitsCollection)[i]->GetEdep();
	  }
	}
	else{
	  if (detName == "OuterPad"){
	    Energy_pad_protons[PadNumber]+=(*hitsCollection)[i]->GetEdep(); 
	    Energy_pad_total[PadNumber]+=(*hitsCollection)[i]->GetEdep();
	  }
	}
	}
    
      else if((*hitsCollection)[i]->GetParticleID()==-11){
	energyOnGas_betas+=(*hitsCollection)[i]->GetEdep();
	energyOnGas_total+=(*hitsCollection)[i]->GetEdep();
	stepSumLengthOnGas_betas+= (*hitsCollection)[i]->GetStepLength();
	if(FirstBetaDetectorHit){
	  BetaKineticEnergy=(*hitsCollection)[i]->GetParticleKineticEnergy();
	  FirstBetaDetectorHit=false;		  
	}

	G4String detName=(*hitsCollection)[i]->GetDetName();
	Int_t PadNumber=(*hitsCollection)[i]->GetDetID();

	if(innerPadFlag == "on"){

	  if(detName == "InnerPad" && PadNumber==0){
	    Energy_pad_betas[0]+=(*hitsCollection)[i]->GetEdep();
	    Energy_pad_total[0]+=(*hitsCollection)[i]->GetEdep();
	  }
	else if (detName == "OuterPad"){
	  Energy_pad_betas[PadNumber+1]+=(*hitsCollection)[i]->GetEdep();
	  Energy_pad_total[PadNumber+1]+=(*hitsCollection)[i]->GetEdep();
	} 
	}

	else{
	  if (detName == "OuterPad"){
	    Energy_pad_betas[PadNumber]+=(*hitsCollection)[i]->GetEdep();
	    Energy_pad_total[PadNumber]+=(*hitsCollection)[i]->GetEdep();
	  }
	}
      }
      else if((*hitsCollection)[i]->GetParticleCharge()==2){
	energyOnGas_alphas+=(*hitsCollection)[i]->GetEdep();
	energyOnGas_total+=(*hitsCollection)[i]->GetEdep();
	alpha_last_positionX=(*hitsCollection)[i]->GetPostPos().x();
	alpha_last_positionY=(*hitsCollection)[i]->GetPostPos().y();
	alpha_last_positionZ=(*hitsCollection)[i]->GetPostPos().z();
	stepSumLengthOnGas_alphas+= (*hitsCollection)[i]->GetStepLength();
      }
      
      
      //G4cout<<"Hits Track ID "<<(*hitsCollection)[j]->GetTrackID()<<" "<<(*hitsCollection)[j]->GetParentID()
      //		<<" "<<(*hitsCollection)[j]->GetParticleID()<<G4endl;
      
    }

    

    //G4cout<<"FillingTree..."<<G4endl;

    if(BetaKineticEnergyDetector>BetaKineticEnergyVeto)
      BetaKineticEnergy=BetaKineticEnergyDetector;
    else
      BetaKineticEnergy=BetaKineticEnergyVeto;
    
    if(ProtonKineticEnergyDetector>ProtonKineticEnergyVeto)
      ProtonKineticEnergy=ProtonKineticEnergyDetector;
    else
      ProtonKineticEnergy=ProtonKineticEnergyVeto;
    
    if(RecoilKineticEnergyDetector>RecoilKineticEnergyVeto)
      RecoilKineticEnergy=RecoilKineticEnergyDetector;
    else
      RecoilKineticEnergy=RecoilKineticEnergyVeto;
    
    theData->SetEnergyOnGas_beam(energyOnGas_beam);
    theData->SetEnergyOnGas_protons(energyOnGas_protons);
    theData->SetBetaKineticEnergy(BetaKineticEnergy);
    theData->SetProtonKineticEnergy(ProtonKineticEnergy);
    theData->SetRecoilKineticEnergy(RecoilKineticEnergy);
    for(Int_t j=0;j<NPADS;j++){
      theData->SetEnergyOnPad_protons(j,Energy_pad_protons[j]);
      theData->SetEnergyOnPad_betas(j,Energy_pad_betas[j]);
      theData->SetEnergyOnPad_total(j,Energy_pad_total[j]);
      if(Energy_pad_protons[j]>0)
	proton_multiplicity++;
      if(Energy_pad_betas[j]>0)
	beta_multiplicity++;
      if(Energy_pad_total[j]>0)
	total_multiplicity++;
    }

        for(Int_t j=0;j<NVETO;j++){
      theData->SetEnergyOnVeto_protons(j,Energy_veto_protons[j]);
      theData->SetEnergyOnVeto_betas(j,Energy_veto_betas[j]);
      theData->SetEnergyOnVeto_total(j,Energy_veto_total[j]);
     if(Energy_veto_protons[j]>0)
	proton_Vetomultiplicity++;
      if(Energy_veto_betas[j]>0)
	beta_Vetomultiplicity++;
      if(Energy_veto_total[j]>0)
	total_Vetomultiplicity++;
    }


    theData->SetEnergyOnGas_betas(energyOnGas_betas);
    theData->SetEnergyOnGas_alphas(energyOnGas_alphas);
    theData->SetEnergyOnGas_total(energyOnGas_total);
    theData->SetEnergyOnDegrader(energyOnDegrader);

    theData->SetEnergyOnVeto(energyOnVeto);
    theData->SetEnergyOnVeto_protons(energyOnVeto_protons);
    theData->SetEnergyOnVeto_betas(energyOnVeto_betas);

    theData->SetStepSumLengthOnGas_beam(stepSumLengthOnGas_beam);
    theData->SetStepSumLengthOnGas_protons(stepSumLengthOnGas_protons);
    theData->SetStepSumLengthOnGas_betas(stepSumLengthOnGas_betas);
    theData->SetStepSumLengthOnGas_alphas(stepSumLengthOnGas_alphas);
    theData->SetBeamLastPositionX(beam_last_positionX);
    theData->SetBeamLastPositionY(beam_last_positionY);
    theData->SetBeamLastPositionZ(beam_last_positionZ);
    theData->SetProtonLastPositionX(proton_last_positionX);
    theData->SetProtonLastPositionY(proton_last_positionY);
    theData->SetProtonLastPositionZ(proton_last_positionZ);
    theData->SetAlphaLastPositionX(alpha_last_positionX);
    theData->SetAlphaLastPositionY(alpha_last_positionY);
    theData->SetAlphaLastPositionZ(alpha_last_positionZ);
    theData->SetProtonMultiplicity(proton_multiplicity);
    theData->SetBetaMultiplicity(beta_multiplicity);
    theData->SetTotalMultiplicity(total_multiplicity);
    theData->SetProtonVetoMultiplicity(proton_Vetomultiplicity);
    theData->SetBetaVetoMultiplicity(beta_Vetomultiplicity);
    theData->SetTotalVetoMultiplicity(total_Vetomultiplicity);

    theData->SetRunID(GetTheRunID());
    theData->SetEventID(GetTheEventID());

    eventTree->Fill();

    delete [] strideOrdinal;
    //for(G4int k=0;k<NbTrajectories;k++)
    //delete [] simpleTrack[k];

}



void ProtonDetectorROOTAnalysis::UserSteppingAction(const G4Step* aStep){

	  G4Track* myTrack = aStep->GetTrack();
	  G4ThreeVector prePoint = aStep->GetPreStepPoint()->GetPosition();
	  G4ThreeVector postPoint = aStep->GetPostStepPoint()->GetPosition();

	  if(storeTracksFlag == "on") {
	      theTracks->SetXCoord(postPoint.x());
	      theTracks->SetYCoord(postPoint.y());
	      theTracks->SetZCoord(postPoint.z());
	      theTracks->SetXPreCoord(prePoint.x());
	      theTracks->SetYPreCoord(prePoint.y());
	      theTracks->SetZPreCoord(prePoint.z());
	      theTracks->SetEnergyStep(aStep->GetTotalEnergyDeposit());
	      theTracks->SetTrackID(myTrack->GetTrackID());
	      theTracks->SetParentTrackID(myTrack->GetParentID());
	      theTracks->SetEventID(GetTheEventID());
	      theTracks->SetRunID(GetTheRunID());
	      tracksTree->Fill();
	    }

}








void ProtonDetectorROOTAnalysis::OnceAWhileDoIt(const G4bool DoItNow) {
  //
  //
  //
  time_t Now = time(0); // get the current time (measured in seconds)
  if ( (!DoItNow) && (LastDoItTime > (Now - 10)) ) return; // every 10 seconds
  LastDoItTime = Now;

  if (gSystem) gSystem->ProcessEvents();

}


