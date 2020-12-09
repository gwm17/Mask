#ifndef G3VEC_H
#define G3VEC_H

#include <cmath>

class G3Vec {
public:
	G3Vec();
	G3Vec(double x, double y, double z);
	~G3Vec();

	void SetVectorCartesian(double x, double y, double z);
	void SetVectorSpherical(double r, double theta, double phi);
	inline double GetX() const { return m_data[0]; };
	inline double GetY() const { return m_data[1]; };
	inline double GetZ() const { return m_data[2]; };
	inline double GetRho() const { return std::sqrt(std::pow(m_data[0], 2.0) + std::pow(m_data[1], 2.0)); };
	inline double GetR() const { return std::sqrt(std::pow(m_data[0], 2.0) + std::pow(m_data[1], 2.0) + std::pow(m_data[2], 2.0)); }
	inline double GetTheta() const { return Atan2(GetRho(), GetZ()); };
	inline double GetPhi() const {
		double phi = Atan2(GetY(), GetX());
		if(phi < 0) phi += M_PI*2.0;
		return phi;
	};

	inline const double operator[](int index) const { return index>2 || index<0 ? 0.0 : m_data[index]; };
	inline G3Vec& operator=(const G3Vec& rhs) { SetVectorCartesian(rhs.GetX(), rhs.GetY(), rhs.GetZ()); return *this; };
	inline G3Vec operator+(const G3Vec& rhs) const { return G3Vec(this->GetX()+rhs.GetX(), this->GetY()+rhs.GetY(), this->GetZ()+rhs.GetZ()); };
	inline G3Vec operator-(const G3Vec& rhs) const { return G3Vec(this->GetX()-rhs.GetX(), this->GetY()-rhs.GetY(), this->GetZ()-rhs.GetZ()); };


	double Dot(const G3Vec& rhs) const;
	G3Vec Cross(const G3Vec& rhs) const;



private:

	inline double Atan2(double y, double x) const {
		if(x != 0.0) return std::atan2(y, x);
		else if(y > 0.0) return M_PI/2.0;
		else if(y < 0.0) return -M_PI/2.0;
		else return 0.0;
	}

	double m_data[3];

};

#endif