/*
	Nucleus.h
	Nucleus is a derived class of Vec4. A nucleus is the kinematics is essentially a 4 vector with the
	additional properties of the number of total nucleons (A), the number of protons (Z), a ground state mass,
	an exctitation energy, and an isotopic symbol.

	--GWM Jan 2021
*/
#ifndef NUCLEUS_H
#define NUCLEUS_H

#include "Vec4.h"
#include <string>

namespace Mask {

class Nucleus : public Vec4 {
public:
	Nucleus();
	Nucleus(int Z, int A);
	Nucleus(int Z, int A, double px, double py, double pz, double E);
	~Nucleus();
	bool SetIsotope(int Z, int A);
	inline int GetZ() const { return m_z; };
	inline int GetA() const { return m_a; };
	inline double GetExcitationEnergy() const { return GetInvMass() - m_gs_mass; };
	inline double GetGroundStateMass() const { return m_gs_mass; };
	inline const char* GetIsotopicSymbol() const { return m_symbol.c_str(); };

	inline Nucleus& operator=(const Nucleus& rhs) {
		SetIsotope(rhs.GetZ(), rhs.GetA());
		SetVectorCartesian(rhs.GetPx(), rhs.GetPy(), rhs.GetPz(), rhs.GetE());
		return *this;
	};

	//Conservation of nucleons and momentum
	inline Nucleus operator+(const Nucleus& daughter) {
		return Nucleus(GetZ()+daughter.GetZ(), GetA()+daughter.GetA(), GetPx()+daughter.GetPx(), GetPy()+daughter.GetPy(), GetPz()+daughter.GetPz(), GetE()+daughter.GetE()); 
	};
	inline Nucleus operator-(const Nucleus& daughter) {
		return (GetZ() - daughter.GetZ()) <= 0 || (GetA() - daughter.GetA()) <= 0 ? Nucleus() :
				Nucleus(GetZ()-daughter.GetZ(), GetA() - daughter.GetA(), GetPx()-daughter.GetPx(), GetPy()-daughter.GetPy(), GetPz()-daughter.GetPz(), GetE()-daughter.GetE());

	};

private:
	int m_z, m_a;
	double m_gs_mass;
	std::string m_symbol;

};

};

#endif