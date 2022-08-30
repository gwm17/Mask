/*

  Class which represents a single MMM detector in the SABRE array at FSU. Origial code by KGH, re-written by
  GWM.

  Distances in meters, angles in radians.

  The channel arrays have four points, one for each corner. The corners are
  as follows, as if looking BACK along beam (i.e. from the target's pov):

  0---------------------1
  |                     |
  |                     |      x
  |                     |      <-----
  |                     |      		|
  |                     |      		|
  3---------------------2      		y
                               (z is hence positive along beam direction) 

  The channel numbers, also as looking back from target pov, are:

  >> rings are 0 -- 15 from inner to outer:

    15 -------------------
    14 -------------------
    13 -------------------
       .
       .
       .
     2 -------------------
     1 -------------------
     0 -------------------

  >> wedges are 0 -- 7 moving counterclockwise:

      7 6 ... 1 0
     | | |   | | |
     | | |   | | |
     | | |   | | |
     | | |   | | |
     | | |   | | |
     | | |   | | |


  >> Note that the detector starts centered on the x-axis (central phi = 0) untilted, and then is rotated to wherever the frick
  	 it is supposed to go; phi = 90 is centered on y axis, pointing down towards the bottom of the scattering chamber

  -- GWM, Dec 2020; based on the og code from kgh

*/

#ifndef SABREDETECTOR_H
#define SABREDETECTOR_H

#include <vector>
#include <cmath>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/RotationZ.h"
#include "Math/RotationY.h"
#include "Math/Translation3D.h"

class SabreDetector {
public:

	SabreDetector();
	SabreDetector(int detID, double phiCentral, double tiltFromVert, double zdist, double xdist=0, double ydist=0);
	~SabreDetector();

	/*Return coordinates of the corners of each ring/wedge in SABRE*/
	const ROOT::Math::XYZPoint& GetRingFlatCoords(int ch, int corner) const { return m_flatRingCoords[ch][corner]; }
	const ROOT::Math::XYZPoint& GetWedgeFlatCoords(int ch, int corner) const { return m_flatWedgeCoords[ch][corner]; }
	const ROOT::Math::XYZPoint& GetRingTiltCoords(int ch, int corner) const { return m_tiltRingCoords[ch][corner]; }
	const ROOT::Math::XYZPoint& GetWedgeTiltCoords(int ch, int corner) const { return m_tiltWedgeCoords[ch][corner]; }

	ROOT::Math::XYZPoint GetTrajectoryCoordinates(double theta, double phi);
	std::pair<int, int> GetTrajectoryRingWedge(double theta, double phi);
	ROOT::Math::XYZPoint GetHitCoordinates(int ringch, int wedgech);

	int GetNumberOfWedges() { return s_nWedges; }
	int GetNumberOfRings() { return s_nRings; }
	ROOT::Math::XYZVector GetNormTilted() { return m_zRotation*(m_yRotation*m_norm); }
	int GetDetectorID() { return m_detectorID; }

private:
	void CalculateCorners();

	/*Performs the transformation to the tilted,rotated,translated frame of the SABRE detector*/
	ROOT::Math::XYZPoint Transform(const ROOT::Math::XYZPoint& vector) { return m_translation*(m_zRotation*(m_yRotation*vector)); };

	/*Determine if a given channel/corner combo is valid*/
	bool CheckRingChannel(int ch) { return (ch<s_nRings && ch>=0) ? true : false; };
	bool CheckWedgeChannel(int ch) { return (ch<s_nWedges && ch >=0) ? true : false; };
	bool CheckCorner(int corner) { return (corner < 4 && corner >=0) ? true : false; };
	bool CheckRingLocation(int ch, int corner) { return CheckRingChannel(ch) && CheckCorner(corner); };
	bool CheckWedgeLocation(int ch, int corner) { return CheckWedgeChannel(ch) && CheckCorner(corner); };

	/*
		For all of the calculations, need a limit precision to determine if values are actually equal or not
		Here the approx. size of the strip spacing is used as the precision.
	*/
	bool CheckPositionEqual(double val1,double val2) { return fabs(val1-val2) > s_positionTol ? false : true; };
	bool CheckAngleEqual(double val1,double val2) { return fabs(val1-val2) > s_angularTol ? false : true; };

	/*Determine if a hit is within the bulk detector*/
	bool IsInside(double r, double phi)
	{ 
		double phi_1 = s_deltaPhiTotal/2.0;
		double phi_2 = M_PI*2.0 - s_deltaPhiTotal/2.0;
		return (((r > s_innerR && r < s_outerR) || CheckPositionEqual(r, s_innerR) || CheckPositionEqual(r, s_outerR)) 
				&& (phi > phi_2 || phi < phi_1 || CheckAngleEqual(phi, phi_1) || CheckAngleEqual(phi, phi_2)));
	};

	/*
		For a given radius/phi are you inside of a given ring/wedge channel,
		or are you on the spacing between these channels
	*/
	bool IsRing(double r, int ringch)
	{
		double ringtop = s_innerR + s_deltaR*(ringch + 1);
		double ringbottom = s_innerR + s_deltaR*(ringch);
		return (r>ringbottom && r<ringtop); 
	};

	inline bool IsRingTopEdge(double r, int ringch)
	{
		double ringtop = s_innerR + s_deltaR*(ringch + 1);
		return CheckPositionEqual(r, ringtop); 
	};

	inline bool IsRingBottomEdge(double r, int ringch)
	{
		double ringbottom = s_innerR + s_deltaR*(ringch);
		return CheckPositionEqual(r, ringbottom); 
	};

	inline bool IsWedge(double phi, int wedgech)
	{
		double wedgetop = -s_deltaPhiTotal/2.0 + s_deltaPhi*(wedgech+1);
		double wedgebottom = -s_deltaPhiTotal/2.0 + s_deltaPhi*(wedgech);
		return ((phi>wedgebottom && phi<wedgetop));
	};

	inline bool IsWedgeTopEdge(double phi, int wedgech)
	{
		double wedgetop = -s_deltaPhiTotal/2.0 + s_deltaPhi*(wedgech+1);
		return CheckAngleEqual(phi, wedgetop);
	}

	inline bool IsWedgeBottomEdge(double phi, int wedgech)
	{
		double wedgebottom = -s_deltaPhiTotal/2.0 + s_deltaPhi*(wedgech);
		return CheckAngleEqual(phi, wedgebottom);
	}

	/*Class data*/
	double m_centerPhi;
	double m_tilt;
	ROOT::Math::Translation3D m_translation;
	ROOT::Math::RotationY m_yRotation;
	ROOT::Math::RotationZ m_zRotation;
	ROOT::Math::XYZVector m_norm;
	int m_detectorID;

	std::vector<std::vector<ROOT::Math::XYZPoint>> m_flatRingCoords, m_flatWedgeCoords;
	std::vector<std::vector<ROOT::Math::XYZPoint>> m_tiltRingCoords, m_tiltWedgeCoords;

	/*Class constants*/
	static constexpr double s_deg2rad = M_PI/180.0;
	static constexpr int s_nRings = 16;
	static constexpr int s_nWedges = 8;
	static constexpr double s_outerR = 0.1351;
	static constexpr double s_innerR = 0.0326;
	static constexpr double s_deltaR = (s_outerR - s_innerR) / s_nRings;
	static constexpr double s_deltaPhiTotal = 54.4 * s_deg2rad;
	static constexpr double s_deltaPhi = (s_deltaPhiTotal / s_nWedges);
	/*These are implicitly the width of the spacing between detector active strips*/
	static constexpr double s_positionTol = 0.0001; //0.1 mm position tolerance
	static constexpr double s_angularTol = 0.1*M_PI/180.0; // 0.1 degree angular tolerance
};


#endif
