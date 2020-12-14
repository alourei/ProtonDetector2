/*
 * ProtonDetectorSteppingAction.hh
 *
 *  Created on: Dec 13, 2013
 *      Author: perezlou
 */

#ifndef PROTONDETECTORSTEPPINGACTION_HH_
#define PROTONDETECTORSTEPPINGACTION_HH_

#include <G4UserSteppingAction.hh>
#include "globals.hh"

class ProtonDetectorConstruction;
class ProtonDetectorEventAction;

class ProtonDetectorSteppingAction: public G4UserSteppingAction {
public:
	ProtonDetectorSteppingAction();
	~ProtonDetectorSteppingAction();
	void UserSteppingAction(const G4Step*);
private:
	ProtonDetectorConstruction *detector;
	ProtonDetectorEventAction *eventAction;
};

#endif /* PROTONDETECTORSTEPPINGACTION_HH_ */
