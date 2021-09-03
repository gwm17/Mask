/*

MassLookup.h
Generates a map for isotopic masses using AMDC data; subtracts away
electron mass from the atomic mass by default. Creates a static global instance
of this map (MASS) for use throughout code it is included into.

Written by G.W. McCann Aug. 2020

*/
#include "MassLookup.h"
#include "KinematicsExceptions.h"


/*
  Read in AMDC mass file, preformated to remove excess info. Here assumes that by default
  the file is in a local directory etc/
*/

MassLookup* MassLookup::s_instance = nullptr;

MassLookup::MassLookup() {
	std::ifstream massfile("./etc/mass.txt");
	if(massfile.is_open()) {
		std::string junk, A, element;
		int Z;
		double atomicMassBig, atomicMassSmall, isotopicMass;
		getline(massfile,junk);
		getline(massfile,junk);
		while(massfile>>junk) {
			massfile>>Z>>A>>element>>atomicMassBig>>atomicMassSmall;
			isotopicMass = (atomicMassBig + atomicMassSmall*1e-6 - Z*electron_mass)*u_to_mev;
			std::string key = "("+std::to_string(Z)+","+A+")";
			massTable[key] = isotopicMass;
			elementTable[Z] = element;
		}
	} else {
	  throw MassFileException();
	}
}

MassLookup::~MassLookup() {}

//Returns nuclear mass in MeV
double MassLookup::FindMass(int Z, int A) {
  	std::string key = "("+std::to_string(Z)+","+std::to_string(A)+")";
  	auto data = massTable.find(key);
  	if(data == massTable.end()) {
  	  	throw MassException();
  	}
  	return data->second;
}

//returns element symbol
std::string MassLookup::FindSymbol(int Z, int A) {
  	auto data = elementTable.find(Z);
  	if(data == elementTable.end()) {
  	  	throw MassException();
  	}
  	std::string fullsymbol = std::to_string(A) + data->second;
  	return fullsymbol;
}
