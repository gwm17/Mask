#ifndef SABREDETECTOR_H
#define SABREDETECTOR_H

#include <vector>
#include <cmath>

#include "G3Vec.h"
#include "GRotation.h"

class SabreDetector {
public:

	SabreDetector();
	SabreDetector(double Rin, double Rout, double deltaPhi_flat, double phiCentral, double tiltFromVert, double zdist, double xdist=0, double ydist=0);
	~SabreDetector();
	inline G3Vec GetRingFlatCoords(int ch, int corner) { return CheckRingLocation(ch, corner) ? m_ringCoords_flat[ch][corner] : G3Vec(); };
	inline G3Vec GetWedgeFlatCoords(int ch, int corner) { return CheckWedgeLocation(ch, corner) ? m_wedgeCoords_flat[ch][corner] : G3Vec(); };
	inline G3Vec GetRingTiltCoords(int ch, int corner) { return CheckRingLocation(ch, corner) ? m_ringCoords_tilt[ch][corner] : G3Vec(); };
	inline G3Vec GetWedgeTiltCoords(int ch, int corner) { return CheckWedgeLocation(ch, corner) ? m_wedgeCoords_tilt[ch][corner] : G3Vec(); };
	G3Vec GetTrajectoryCoordinates(double theta, double phi);
	G3Vec GetHitCoordinates(int ringch, int wedgech);

	inline int GetNumberOfWedges() { return m_nWedges; };
	inline int GetNumberOfRings() { return m_nRings; };
	inline double GetInnerRadius() { return m_Rinner; };
	inline double GetOuterRadius() { return m_Router; };
	inline double GetPhiCentral() { return m_phiCentral; };
	inline double GetTiltAngle() { return m_tilt; };
	inline G3Vec GetTranslation() { return m_translation; };


private:

	static constexpr int m_nRings = 16;
	static constexpr int m_nWedges = 8;
	static constexpr double deg2rad = M_PI/180.0;
	static constexpr double POSITION_TOL = 0.0001;
	static constexpr double ANGULAR_TOL = M_PI/180.0;

	void CalculateCorners();
	inline G3Vec TransformToTiltedFrame(G3Vec& vector) { return (m_ZRot*(m_YRot*vector)) + m_translation; };

	inline bool CheckRingChannel(int ch) { return (ch<m_nRings && ch>=0) ? true : false; };
	inline bool CheckWedgeChannel(int ch) { return (ch<m_nWedges && ch >=0) ? true : false; };
	inline bool CheckCorner(int corner) { return (corner < 4 && corner >=0) ? true : false; };
	inline bool CheckRingLocation(int ch, int corner) { return CheckRingChannel(ch) && CheckCorner(corner); };
	inline bool CheckWedgeLocation(int ch, int corner) { return CheckWedgeChannel(ch) && CheckCorner(corner); };

	inline bool CheckPositionEqual(double val1,double val2) { return fabs(val1-val2) > POSITION_TOL ? false : true; };
	inline bool CheckAngleEqual(double val1,double val2) { return fabs(val1-val2) > ANGULAR_TOL ? false : true; };

	inline bool IsInside(double r, double phi) { 
		double phi_1 = m_deltaPhi_flat/2.0;
		double phi_2 = M_PI*2.0 - m_deltaPhi_flat/2.0;
		return (((r > m_Rinner && r < m_Router) || CheckPositionEqual(r, m_Rinner) || CheckPositionEqual(r, m_Router)) && (phi > phi_2 || phi < phi_1 || CheckAngleEqual(phi, phi_1) || CheckAngleEqual(phi, phi_2)));
	};

	double m_Router, m_Rinner, m_deltaPhi_flat, m_phiCentral, m_tilt;
	G3Vec m_translation;
	GYRotation m_YRot;
	GZRotation m_ZRot;
	double m_deltaR_flat, m_deltaR_flat_ring;

	std::vector<std::vector<G3Vec>> m_ringCoords_flat, m_wedgeCoords_flat;
	std::vector<std::vector<G3Vec>> m_ringCoords_tilt, m_wedgeCoords_tilt;

};

#endif