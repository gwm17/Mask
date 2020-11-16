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
	inline double GetTheta() const {return GetP() == 0.0 ? 0.0 : acos(GetPz()/GetP());};
	inline double GetPhi() const {
		if(GetPx() == 0) return M_PI/2.0;
		double phi = std::atan(GetPy()/GetPx());
		if(GetPx()<0) phi += M_PI;
		else if(GetPy()<0) phi += 2.0*M_PI;
		return phi;
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
	double m_data[4];
	double m_boost[3];

};

#endif