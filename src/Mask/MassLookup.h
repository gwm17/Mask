/*

MassLookup.h
Generates a map for isotopic masses using AMDC data; subtracts away
electron mass from the atomic mass by default. Creates a static global instance
of this map (MASS) for use throughout code it is included into.

Written by G.W. McCann Aug. 2020

Converted to true singleton to simplify usage -- Aug. 2021 GWM
*/
#ifndef MASS_LOOKUP_H
#define MASS_LOOKUP_H

#include <fstream>
#include <string>
#include <unordered_map>

namespace Mask {

	class MassLookup
	{
	public:

		struct KeyPair
		{
			uint32_t Z;
			uint32_t A;

			//Use szudzik pairing method to make unqiue key out of two unsigned ints. Use size_t as extra safety.
			std::size_t GetID()
			{
				return Z >= A ? Z*Z + Z + A : A*A + Z;
			}
		};

		~MassLookup();
		double FindMass(uint32_t Z, uint32_t A);
		double FindMassU(uint32_t Z, uint32_t A) { return FindMass(Z, A)/u_to_mev; }
		std::string FindSymbol(uint32_t Z, uint32_t A);
	
		static MassLookup& GetInstance() { return *s_instance; }
	
	private:
		MassLookup();

		static MassLookup* s_instance;
		std::unordered_map<std::size_t, double> massTable;
		std::unordered_map<std::size_t, std::string> elementTable;
	
		//constants
		static constexpr double u_to_mev = 931.4940954;
		static constexpr double electron_mass = 0.000548579909;
	    
	};

}

#endif
