#include "Nucleus.h"
#include "MassLookup.h"

Nucleus::Nucleus () :
	G4Vec(), m_z(0), m_a(0), m_gs_mass(0), m_symbol("")
{
}

Nucleus::Nucleus(int Z, int A) :
	G4Vec(), m_z(Z), m_a(A)
{
	m_gs_mass = MASS.FindMass(Z, A);
	m_symbol = MASS.FindSymbol(Z, A);
	SetVectorCartesian(0,0,0,m_gs_mass);
}

Nucleus::Nucleus(int Z, int A, double px, double py, double pz, double E) :
	G4Vec(px, py, pz, E), m_z(Z), m_a(A)
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