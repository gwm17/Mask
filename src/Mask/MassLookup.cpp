/*

MassLookup.h
Generates a map for isotopic masses using AMDC data; subtracts away
electron mass from the atomic mass by default. Creates a static global instance
of this map (MASS) for use throughout code it is included into.

Written by G.W. McCann Aug. 2020

*/
#include "MassLookup.h"
#include "KinematicsExceptions.h"
#include <sstream>

namespace Mask {

	MassLookup* MassLookup::s_instance = new MassLookup();

	MassLookup::MassLookup()
	{
		std::ifstream massfile("etc/mass.txt");
		KeyPair key;
		if(massfile.is_open())
		{
			std::string junk, element;
			double atomicMassBig, atomicMassSmall, isotopicMass;
			getline(massfile,junk);
			getline(massfile,junk);
			while(massfile>>junk)
			{
				massfile>>key.Z>>key.A>>element>>atomicMassBig>>atomicMassSmall;
				isotopicMass = (atomicMassBig + atomicMassSmall*1e-6 - key.Z*electron_mass)*u_to_mev;
				massTable[key.GetID()] = isotopicMass;
				elementTable[key.GetID()] = std::to_string(key.A) + element;
			}
		}
		else
			throw MassFileException();
	}
	
	MassLookup::~MassLookup() {}
	
	//Returns nuclear mass in MeV
	double MassLookup::FindMass(uint32_t Z, uint32_t A)
	{
		KeyPair key({Z, A});
		auto data = massTable.find(key.GetID());
		if(data == massTable.end())
			throw MassException();
	
		return data->second;
	}
	
	//returns element symbol
	std::string MassLookup::FindSymbol(uint32_t Z, uint32_t A)
	{
		KeyPair key({Z, A});
		auto data = elementTable.find(key.GetID());
		if(data == elementTable.end())
			throw MassException();
	
		return data->second;
	}

}