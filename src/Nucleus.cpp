#include "Nucleus.h"
#include "MassLookup.h"

namespace Mask {

Nucleus::Nucleus () :
	Vec4(), m_z(0), m_a(0), m_gs_mass(0), m_symbol("")
{
}

Nucleus::Nucleus(int Z, int A) :
	Vec4(), m_z(Z), m_a(A)
{
	m_gs_mass = MASS.FindMass(Z, A);
	m_symbol = MASS.FindSymbol(Z, A);
	SetVectorCartesian(0,0,0,m_gs_mass);
}

Nucleus::Nucleus(int Z, int A, double px, double py, double pz, double E) :
	Vec4(px, py, pz, E), m_z(Z), m_a(A)
{
	m_gs_mass = MASS.FindMass(Z, A);
	m_symbol = MASS.FindSymbol(Z, A);
}

Nucleus::~Nucleus() {}

bool Nucleus::SetIsotope(int Z, int A) {
	if(Z>A) return false;
	
	m_z = Z;
	m_a = A;
	m_gs_mass = MASS.FindMass(Z, A);
	m_symbol = MASS.FindSymbol(Z, A);
	SetVectorCartesian(0,0,0,m_gs_mass);
	return true;
}

};