/*
	Class which represents a 4-momentum vector. Can perform vector addition, subtraction, dot product
	and generate a boost vector to its rest frame as well as apply a boost to itself.

	--GWM Dec 2020.
*/
#ifndef VEC4_H
#define VEC4_H

#include <cmath>

namespace Mask {

class Vec4 {
public:
	Vec4();
	Vec4(double px, double py, double pz, double E);
	virtual ~Vec4();
	void SetVectorCartesian(double px, double py, double pz, double E);
	void SetVectorSpherical(double theta, double phi, double p, double E);

	inline double GetE() const {return m_data[3];};
	inline double GetPx() const {return m_data[0];};
	inline double GetPy() const {return m_data[1];};
	inline double GetPz() const {return m_data[2];};
	inline double GetP() const {return std::sqrt(m_data[0]*m_data[0] + m_data[1]*m_data[1] + m_data[2]*m_data[2]);};
	inline double GetPxy() const {return std::sqrt(m_data[0]*m_data[0] + m_data[1]*m_data[1]); };
	inline double GetTheta() const {return GetPxy() == 0.0 && GetPz() == 0.0 ? 0.0 : Atan2(GetPxy(), GetPz());};
	inline double GetPhi() const {
		double phi = Atan2(GetPy(), GetPx());
		if(phi<0) phi += 2.0*M_PI;
		return GetPx() == 0.0 && GetPy() == 0.0 ? 0.0 : phi;
	};

	inline double GetInvMass() const {return std::sqrt(GetE()*GetE() - GetP()*GetP());};
	inline double GetKE() const {return GetE() - GetInvMass();};
	inline const double* GetBoost() const {return &m_boost[0];};

	void ApplyBoost(const double* boost);

	//Only intended for use in looping access!
	inline const double operator[] (int index) const {return index>3 || index < 0 ? 0.0 : m_data[index];};

	inline Vec4& operator=(const Vec4& rhs) {SetVectorCartesian(rhs.GetPx(), rhs.GetPy(), rhs.GetPz(), rhs.GetE()); return *this;};
	inline Vec4 operator+(const Vec4& rhs) const {return Vec4(rhs.GetPx()+GetPx(), rhs.GetPy()+GetPy(), rhs.GetPz()+GetPz(), rhs.GetE()+GetE());};
	inline Vec4 operator-(const Vec4& rhs) const {return Vec4(rhs.GetPx()-GetPx(), rhs.GetPy()-GetPy(), rhs.GetPz()-GetPz(), rhs.GetE()-GetE());};
	
	double Dot(const Vec4& rhs) const;
	Vec4 Cross(const Vec4& rhs) const;

private:
	void CalcBoostToCM();

	//use instead of std::atan2. Better controll over x=0
	inline double Atan2(double y, double x) const { 
		if(x != 0) return std::atan2(y, x);
		else if( y > 0 ) return M_PI/2.0;
		else if( y < 0 ) return -M_PI/2.0;
		else return 0.0;
	};

	double m_data[4];
	double m_boost[3];

};

};

#endif