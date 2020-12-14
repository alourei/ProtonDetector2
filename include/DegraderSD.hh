/*
 * DegraderSD.hh
 *
 *  Created on: Dec 20, 2013
 *      Author: perezlou
 */

#ifndef DEGRADERSD_HH_
#define DEGRADERSD_HH_

#include <ProtonDetectorGeantHit.hh>

class G4Step;
class G4HCofThisEvent;
#include <G4VSensitiveDetector.hh>

class DegraderSD: public G4VSensitiveDetector {
public:
	DegraderSD(G4String name);
	~DegraderSD();

	void Initialize(G4HCofThisEvent*);
	G4bool ProcessHits(G4Step*,G4TouchableHistory*);
	void EndOfEvent(G4HCofThisEvent*);


private:

	ProtonDetectorGeantHitsCollection* hitsCollection; //Geant step-like hits collect.

};

#endif /* DEGRADERSD_HH_ */
