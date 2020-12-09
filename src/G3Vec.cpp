#include "G3Vec.h"

G3Vec::G3Vec() {
	m_data[0] = 0.;
	m_data[1] = 0.;
	m_data[2] = 0.;
}

G3Vec::G3Vec(double x, double y, double z) {
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
}

G3Vec::~G3Vec() {}

void G3Vec::SetVectorCartesian(double x, double y, double z) {
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
}

void G3Vec::SetVectorSpherical(double r, double theta, double phi) {
	m_data[0] = r*std::cos(phi)*std::sin(theta);
	m_data[1] = r*std::sin(phi)*std::sin(theta);
	m_data[2] = r*std::cos(theta);
}

double G3Vec::Dot(const G3Vec& rhs) const {
	return GetX()*rhs.GetX() + GetY()*rhs.GetY() + GetZ()*rhs.GetZ();
}

//Unimplemented
G3Vec G3Vec::Cross(const G3Vec& rhs) const {
	return G3Vec(0.,0.,0.);
}