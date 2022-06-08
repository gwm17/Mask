/*
	Class which represents a 4-momentum vector. Can perform vector addition, subtraction, dot product
	and generate a boost vector to its rest frame as well as apply a boost to itself.

	--GWM Dec 2020.
	NOTE: uses (-,-,-,+) metric (same as ROOT convention)
*/
#include "Vec4.h"


namespace Mask {

	Vec4::Vec4() {
		for(auto& val: m_data)
			val = 0.0;
		for(auto& val: m_boost)
			val = 0.0;
	}
	
	Vec4::Vec4(double px, double py, double pz, double E) {
		m_data[0] = px;
		m_data[1] = py;
		m_data[2] = pz;
		m_data[3] = E;
		CalcBoostToCM();
	}
	
	Vec4::~Vec4() {}
	
	void Vec4::SetVectorCartesian(double px, double py, double pz, double E) {
		m_data[0] = px;
		m_data[1] = py;
		m_data[2] = pz;
		m_data[3] = E;
	
		CalcBoostToCM();
	}
	
	void Vec4::SetVectorSpherical(double theta, double phi, double p, double E) {
		m_data[0] = p*cos(phi)*sin(theta);
		m_data[1] = p*sin(phi)*sin(theta);
		m_data[2] = p*cos(theta);
		m_data[3] = E;
		CalcBoostToCM();
	}
	
	void Vec4::CalcBoostToCM() {
		m_boost[0] = m_data[0]/m_data[3];
		m_boost[1] = m_data[1]/m_data[3];
		m_boost[2] = m_data[2]/m_data[3];
	}
	
	void Vec4::ApplyBoost(const double* beta) {
		double beta2 = beta[0]*beta[0] + beta[1]*beta[1] + beta[2]*beta[2];
		double gamma  = 1.0/std::sqrt(1.0 - beta2);
		double bdotp = beta[0]*m_data[0] + beta[1]*m_data[1] + beta[2]*m_data[2];
		double gfactor = beta2>0.0 ? (gamma - 1.0)/beta2 : 0.0;
	
		SetVectorCartesian(GetPx()+gfactor*bdotp*beta[0]+gamma*beta[0]*GetE(),
			      		   GetPy()+gfactor*bdotp*beta[1]+gamma*beta[1]*GetE(),
			    		   GetPz()+gfactor*bdotp*beta[2]+gamma*beta[2]*GetE(),
			      		   gamma*(GetE() + bdotp));
	}
	
	double Vec4::Dot(const Vec4& rhs) const {
		return GetE()*rhs.GetE() - GetPx()*rhs.GetPx() - GetPy()*rhs.GetPy() - GetPz()*rhs.GetPz();
	}
	
	Vec4 Vec4::Cross(const Vec4& rhs) const {
		return Vec4();
	}

}