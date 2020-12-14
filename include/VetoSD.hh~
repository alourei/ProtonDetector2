/*
 * ProtonDetectorSD.hh
 *
 *  Created on: Dec 6, 2013
 *      Author: perezlou
 */

#ifndef PROTONDETECTORSD_HH_
#define PROTONDETECTORSD_HH_

#include <G4VSensitiveDetector.hh>
#include <ProtonDetectorGeantHit.hh>

class G4Step;
class G4HCofThisEvent;

class ProtonDetectorSD: public G4VSensitiveDetector {
public:
	ProtonDetectorSD(G4String name);
	~ProtonDetectorSD();

	void Initialize(G4HCofThisEvent*);
	G4bool ProcessHits(G4Step*,G4TouchableHistory*);
	void EndOfEvent(G4HCofThisEvent*);


private:

	ProtonDetectorGeantHitsCollection* hitsCollection; //Geant step-like hits collect.

};

#endif /* PROTONDETECTORSD_HH_ */
