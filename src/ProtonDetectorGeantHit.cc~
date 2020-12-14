/*
 * ProtonDetectorGeantHit.cc
 *
 *  Created on: Dec 6, 2013
 *      Author: perezlou
 */

#include <ProtonDetectorGeantHit.hh>

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<ProtonDetectorGeantHit> ProtonDetectorGeantHitAllocator;


ProtonDetectorGeantHit::ProtonDetectorGeantHit() {
	// TODO Auto-generated constructor stub

}

ProtonDetectorGeantHit::~ProtonDetectorGeantHit() {
	// TODO Auto-generated destructor stub
}


ProtonDetectorGeantHit::ProtonDetectorGeantHit(const ProtonDetectorGeantHit& right) : G4VHit() {
  //
  // Copy constructor
  //
  trackID = right.trackID;
  parentID = right.parentID;
  edep = right.edep;
  particleCharge = right.particleCharge;
  particleMass = right.particleMass;
  particleID = right.particleID;
  prePos = right.prePos;
  postPos = right.postPos;
  detName = right.detName;
  detID = right.detID;
  preToF = right.preToF;
  postToF = right.postToF;
  stepLength = right.stepLength;
}


const ProtonDetectorGeantHit& ProtonDetectorGeantHit::operator=(const ProtonDetectorGeantHit& right){
  //
  // Operator =
  //
  trackID = right.trackID;
  parentID = right.parentID;
  edep = right.edep;
  particleCharge = right.particleCharge;
  particleMass = right.particleMass;
  particleID = right.particleID;
  prePos = right.prePos;
  postPos = right.postPos;
  detName = right.detName;
  detID = right.detID;
  preToF = right.preToF;
  postToF = right.postToF;
  stepLength = right.stepLength;

  return *this;
}


G4int ProtonDetectorGeantHit::operator==(const ProtonDetectorGeantHit& right) const{
  //
  // Operator ==
  //
  return (this==&right) ? 1 : 0;
}


void ProtonDetectorGeantHit::Draw(){
  //
  // Draws the Hit. A clear red point on the Hit position
  //
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(prePos);
    circle.SetScreenSize(4);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}


void ProtonDetectorGeantHit::Print(){
  //
  // Prints full information about the calGeantHit
  //
  G4cout << "##################################################################"
	 << G4endl
	 << "############     ProtonDetectorGeantHit::Print()     ################" << G4endl
	 << "trackID: " << trackID
	 << "parentID: " << parentID
	 << ", detID: " << detID
	 << ", detName: " << detName << G4endl;
  G4cout << "edep: " << edep  / MeV << " MeV"
	 << ", prePos: " << prePos
	 << ", postPos: " << postPos
	 << ", stepLength: " << stepLength  / mm << " mm"
	 << ", preToF: " << preToF  / ns << " ns"
	 << ", posToF: " << postToF  / ns << " ns" 	 << G4endl;
  G4cout << "##################################################################"
	 << G4endl;

}
