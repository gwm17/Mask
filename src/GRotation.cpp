#include "GRotation.h"

GXRotation::GXRotation() :
m_angle(0)
{
	GenerateMatrix();
}

GXRotation::GXRotation(double angle) :
m_angle(angle)
{
	GenerateMatrix();
}

GXRotation::~GXRotation() {}

void GXRotation::GenerateMatrix() {
	m_matrix[0][0] = 1.0; m_matrix[0][1] = 0.0; m_matrix[0][2] = 0.0;
	m_matrix[1][0] = 0.0; m_matrix[1][1] = std::cos(m_angle); m_matrix[1][2] = -std::sin(m_angle);
	m_matrix[2][0] = 0.0; m_matrix[2][1] = std::sin(m_angle); m_matrix[2][2] = std::cos(m_angle);
}

GYRotation::GYRotation() :
m_angle(0)
{
	GenerateMatrix();
}

GYRotation::GYRotation(double angle) :
m_angle(angle)
{
	GenerateMatrix();
}

GYRotation::~GYRotation() {}

void GYRotation::GenerateMatrix() {
	m_matrix[0][0] = std::cos(m_angle); m_matrix[0][1] = 0.0; m_matrix[0][2] = -std::sin(m_angle);
	m_matrix[1][0] = 0.0; m_matrix[1][1] = 1.0; m_matrix[1][2] = 0.0;
	m_matrix[2][0] = std::sin(m_angle); m_matrix[2][1] = 0.0; m_matrix[2][2] = std::cos(m_angle);
}


GZRotation::GZRotation() :
m_angle(0)
{
	GenerateMatrix();
}

GZRotation::GZRotation(double angle) :
m_angle(angle)
{
	GenerateMatrix();
}

GZRotation::~GZRotation() {}

void GZRotation::GenerateMatrix() {
	m_matrix[0][0] = std::cos(m_angle); m_matrix[0][1] = -std::sin(m_angle); m_matrix[0][2] = 0.0;
	m_matrix[1][0] = std::sin(m_angle); m_matrix[1][1] = std::cos(m_angle); m_matrix[1][2] = 0.0;
	m_matrix[2][0] = 0.0; m_matrix[2][1] = 0.0; m_matrix[2][2] = 1.0;
}