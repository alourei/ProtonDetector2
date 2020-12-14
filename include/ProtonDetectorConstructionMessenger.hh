/*
 * ProtonDetectorConstructionMessenger.hh
 *
 *  Created on: Dec 17, 2013
 *      Author: perezlou
 */

#ifndef PROTONDETECTORCONSTRUCTIONMESSENGER_HH_
#define PROTONDETECTORCONSTRUCTIONMESSENGER_HH_

#include <globals.hh>
#include <G4UImessenger.hh>

class ProtonDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWithoutParameter;


class ProtonDetectorConstructionMessenger: public G4UImessenger {

private:
	ProtonDetectorConstruction *ProtonDetector;
	G4UIdirectory*             ProtonDetectorDir;
	G4UIdirectory*             detDir;

	G4UIcmdWithAString*        detectorGeometryCmd;
	G4UIcmdWithAString* 	   degraderIncludedFlagCmd;
        G4UIcmdWithAString*        gasMaterCmd;
	G4UIcmdWithADoubleAndUnit* xGasBoxCmd;
	G4UIcmdWithADoubleAndUnit* yGasBoxCmd;
	G4UIcmdWithADoubleAndUnit* zGasBoxCmd;
	G4UIcmdWithADoubleAndUnit* radiusGasTubCmd;
	G4UIcmdWithADoubleAndUnit* lengthGasTubCmd;

	G4UIcmdWithADoubleAndUnit* radiusInnerPadCmd;
	G4UIcmdWithADoubleAndUnit* radiusOuterPadCmd;

	//Degrader settings
	G4UIcmdWith3VectorAndUnit* DegraderPosCmd;
	G4UIcmdWithADoubleAndUnit* DegraderThicknessCmd;
	G4UIcmdWithADoubleAndUnit* DegraderAngleCmd;


	G4UIcmdWithoutParameter*   updateCmd;
	G4UIcmdWithoutParameter*   printCmd;



public:
	ProtonDetectorConstructionMessenger(ProtonDetectorConstruction*);
	~ProtonDetectorConstructionMessenger();
	void SetNewValue(G4UIcommand*, G4String);
	//G4String GetCurrentValue(G4UIcommand*);
};

#endif /* PROTONDETECTORCONSTRUCTIONMESSENGER_HH_ */
