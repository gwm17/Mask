#ifndef GROTATION_H
#define GROTATION_H

#include "G3Vec.h"

class GXRotation {
public:
	GXRotation();
	GXRotation(double ang);
	~GXRotation();
	G3Vec Rotate(const G3Vec& vector);
	inline void SetAngle(double ang) { m_angle = ang; GenerateMatrix(); };
	inline GXRotation GetInverse() { return GXRotation(-m_angle); };
	inline G3Vec operator*(const G3Vec& vector) {
		double x = m_matrix[0][0]*vector[0] + m_matrix[0][1]*vector[1] + m_matrix[0][2]*vector[2];
		double y = m_matrix[1][0]*vector[0] + m_matrix[1][1]*vector[1] + m_matrix[1][2]*vector[2];
		double z = m_matrix[2][0]*vector[0] + m_matrix[2][1]*vector[1] + m_matrix[2][2]*vector[2];
		return G3Vec(x, y, z);
	};

private:
	void GenerateMatrix();
	double m_angle;
	double m_matrix[3][3];
};

class GYRotation {
public:
	GYRotation();
	GYRotation(double ang);
	~GYRotation();
	G3Vec Rotate(const G3Vec& vector);
	inline void SetAngle(double ang) { m_angle = ang; GenerateMatrix(); };
	inline GYRotation GetInverse() { return GYRotation(-m_angle); };
	inline G3Vec operator*(const G3Vec& vector) {
		double x = m_matrix[0][0]*vector[0] + m_matrix[0][1]*vector[1] + m_matrix[0][2]*vector[2];
		double y = m_matrix[1][0]*vector[0] + m_matrix[1][1]*vector[1] + m_matrix[1][2]*vector[2];
		double z = m_matrix[2][0]*vector[0] + m_matrix[2][1]*vector[1] + m_matrix[2][2]*vector[2];
		return G3Vec(x, y, z);
	};

private:
	void GenerateMatrix();
	double m_angle;
	double m_matrix[3][3];
};

class GZRotation {
public:
	GZRotation();
	GZRotation(double ang);
	~GZRotation();
	G3Vec Rotate(const G3Vec& vector);
	inline void SetAngle(double ang) { m_angle = ang; GenerateMatrix();};
	inline GZRotation GetInverse() { return GZRotation(-m_angle); };
	inline G3Vec operator*(const G3Vec& vector) {
		double x = m_matrix[0][0]*vector[0] + m_matrix[0][1]*vector[1] + m_matrix[0][2]*vector[2];
		double y = m_matrix[1][0]*vector[0] + m_matrix[1][1]*vector[1] + m_matrix[1][2]*vector[2];
		double z = m_matrix[2][0]*vector[0] + m_matrix[2][1]*vector[1] + m_matrix[2][2]*vector[2];
		return G3Vec(x, y, z);
	};

private:
	void GenerateMatrix();
	double m_angle;
	double m_matrix[3][3];
};

#endif