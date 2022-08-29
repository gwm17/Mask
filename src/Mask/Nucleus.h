/*
	Nucleus.h
	Nucleus is a derived class of Vec4. A nucleus is the kinematics is essentially a 4 vector with the
	additional properties of the number of total nucleons (A), the number of protons (Z), a ground state mass,
	an exctitation energy, and an isotopic symbol.

	--GWM Jan 2021
*/
#ifndef NUCLEUS_H
#define NUCLEUS_H

#include <string>
#include <vector>
#include "Math/Vector4D.h"
#include "MassLookup.h"

namespace Mask {

	struct Nucleus
	{
		void SetVec4Spherical(double theta, double phi, double p, double E)
		{
			vec4.SetPxPyPzE(std::sin(theta)*std::cos(phi)*p,
							std::sin(theta)*std::sin(phi)*p,
							std::cos(theta)*p,
							E
						   );
		}

		double GetKE() const
		{
			return vec4.E() - vec4.M();
		}

		double GetExcitationEnergy() const
		{
			return vec4.M() - groundStateMass;
		}

		uint32_t Z = 0; 
		uint32_t A = 0;
		double groundStateMass = 0.0;
		std::string isotopicSymbol = "";
		double thetaCM = 0.0;
		ROOT::Math::PxPyPzEVector vec4;

		bool isDetected = false;
		double detectedKE = 0.0;
		double detectedTheta = 0.0;
		double detectedPhi = 0.0;
	};

	Nucleus CreateNucleus(int z, int a);

	bool EnforceDictionaryLinked();

};

#endif
