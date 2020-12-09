#ifndef G4VEC_H
#define G4VEC_H

#include <cmath>

class G4Vec {
public:
	G4Vec();
	G4Vec(double px, double py, double pz, double E);
	virtual ~G4Vec();
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

	inline G4Vec& operator=(const G4Vec& rhs) {SetVectorCartesian(rhs.GetPx(), rhs.GetPy(), rhs.GetPz(), rhs.GetE()); return *this;};
	inline G4Vec operator+(const G4Vec& rhs) const {return G4Vec(rhs.GetPx()+GetPx(), rhs.GetPy()+GetPy(), rhs.GetPz()+GetPz(), rhs.GetE()+GetE());};
	inline G4Vec operator-(const G4Vec& rhs) const {return G4Vec(rhs.GetPx()-GetPx(), rhs.GetPy()-GetPy(), rhs.GetPz()-GetPz(), rhs.GetE()-GetE());};
	double Dot(const G4Vec& rhs) const;
	G4Vec Cross(const G4Vec& rhs) const;

private:
	void CalcBoostToCM();
	inline double Atan2(double y, double x) const { 
		if(x != 0) return std::atan2(y, x);
		else if( y > 0 ) return M_PI/2.0;
		else if( y < 0 ) return -M_PI/2.0;
		else return 0.0;
	};
	double m_data[4];
	double m_boost[3];

};

#endif