#ifndef KINEMATICSEXCEPTIONS_H
#define KINEMATICSEXCEPTIONS_H

#include <exception>
#include <stdexcept>
/*
	ELossException
	This is an exception that is thrown by the Energy loss calculator. It is not specific to one particular
	location in EnergyLoss.cpp/.h, however there are really only a couple of ways that these can fail.
	Either your nucleus is not well defined (i.e. doe not exist in the Eloss_Tables.h) or you've tried to run
	it without defining the proper information. See the SetTargetComponentsand GetElectronicStoppingPower functions
	for locations where this could be thrown
*/
struct ELossException : public std::exception {
	const char* what() const noexcept {
		return "Failure to calculate particle energy loss. See KinematicsExceptions.h for documentation.";
	};
};

/*
	MassException
	This is an exception thrown when the MassLookup is queried for an isotope for which it does not have a defined
	mass. The masses are defined in ./etc/mass.txt, which is a condensed version of the AMDC Mass Evaluation from 2017.
*/
struct MassException : public std::exception {
	const char* what() const noexcept {
		return "Unable to find a given isotopic mass. See Kinematics.h for documentation.";
	};
};

/*
	MassFileException
	This is an exception thrown when the MassLookup cannot find the ./etc/mass.txt file.
*/
struct MassFileException : public std::exception {
	const char* what() const noexcept {
		return "Unable to find ./etc/mass.txt. Check that it is present.";
	};
};

/*
	ReactionLayerException
	This is an exception thrown when the ReactionSystem cannot find a good layer in its LayeredTarget
	to set as the layer where the reaction takes place. A good layer is a layer which contains the isotope which
	corresponds to the isotope defined as the "target" in the reaction equation.
*/
struct ReactionLayerException : public std::exception {
	const char* what() const noexcept {
		return "Unable to find a valid layer for reaction in the target. See KinematicsExceptions.h for documentation.";
	};
};

/*
	QValueException
	This is an exception thrown when the Reaction caluclates a decay Q-value that is negative, indicating
	that there is not enough energy to create the decay products.
*/
struct QValueException : public std::exception {
	const char* what() const noexcept {
		return "Q-value is negative for decay calculation. See KinematicsExceptions.h for documentation.";
	};
};



#endif