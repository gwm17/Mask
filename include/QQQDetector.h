/*
	QQQDetector.h
	Class implementing geometry for QQQDetector where the detector is perpendicular to the beam axis.
	Detector is first generated centered on the x-axis (phi=0)

	Coordinate convention : +z is downstream, -z is upstream. +y is vertically up in the lab.
*/
#ifndef QQQDETECTOR_H
#define QQQDETECTOR_H

#include <TRandom3.h>
#include <cmath>

#include "Vec3.h"
#include "Rotation.h"

class QQQDetector {
public:
	QQQDetector(double R_in, double R_out, double deltaPhi, double phiCentral, double z, double x=0, double y=0);
	~QQQDetector();
	inline Mask::Vec3 GetRingCoordinates(int ringch, int corner) { return m_ringCoords[ringch][corner]; };
	inline Mask::Vec3 GetWedgeCoordinates(int wedgech, int corner) { return m_wedgeCoords[wedgech][corner]; };
	Mask::Vec3 GetTrajectoryCoordinates(double theta, double phi);
	std::pair<int, int> GetTrajectoryRingWedge(double theta, double phi);
	Mask::Vec3 GetHitCoordinates(int ringch, int wedgech);

	inline int GetNumberOfRings() { return nrings; };
	inline int GetNumberOfWedges() { return nwedges; };

private:

	inline bool CheckChannel(int ch) { return (ch >=0 && ch < nrings); };
	inline bool CheckCorner(int corner) { return (corner >=0 && corner < 4); };

	void CalculateCorners();
	Mask::Vec3 TransformCoordinates(Mask::Vec3& vector) { return m_ZRot*vector + m_translation; };

	double m_Rinner, m_Router, m_deltaR, m_deltaPhi, m_deltaPhi_per_wedge, m_phiCentral;
	std::vector<std::vector<Mask::Vec3>> m_ringCoords, m_wedgeCoords;
	Mask::Vec3 m_translation;
	Mask::ZRotation m_ZRot;

	static constexpr int nrings = 16;
	static constexpr int nwedges = 16;
	static constexpr double deg2rad  = M_PI/180.0;
};

#endif