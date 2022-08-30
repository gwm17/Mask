#ifndef SX3_DETECTOR_H
#define SX3_DETECTOR_H

// +z is along beam axis
// +y is vertically "downward" in the lab frame

//angles must be in radians, but distances can be whatever
//PROVIDED all input distances are the same

//Front strips from largest y to smallest y

//Back strips from lowest z to highest z

#include <cmath>
#include <vector>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/RotationZ.h"
#include "Mask/RandomGenerator.h"

struct StripHit
{
	int front_strip_index=-1;
	int back_strip_index=-1;
	double front_ratio=0.0;
};

class StripDetector
{
public:
  
	StripDetector(double centerPhi, double centerZ, double centerRho);
	~StripDetector();
	const ROOT::Math::XYZPoint& GetFrontStripCoordinates(int stripch, int corner) const { return m_frontStripCoords[stripch][corner]; }
	const ROOT::Math::XYZPoint& GetBackStripCoordinates(int stripch, int corner) const { return m_backStripCoords[stripch][corner]; }
	const ROOT::Math::XYZPoint& GetRotatedFrontStripCoordinates(int stripch, int corner) const 
	{ 
		return m_rotFrontStripCoords[stripch][corner];
	}
	const ROOT::Math::XYZPoint& GetRotatedBackStripCoordinates(int stripch, int corner) const 
	{ 
		return m_rotBackStripCoords[stripch][corner];
	}
	ROOT::Math::XYZVector GetNormRotated() const { return m_zRotation * m_norm; }

	void SetPixelSmearing(bool isSmearing) { m_isSmearing = isSmearing; }

	ROOT::Math::XYZPoint GetHitCoordinates(int front_stripch, double front_strip_ratio);
	StripHit GetChannelRatio(double theta, double phi);

private:
	bool ValidChannel(int f) { return ((f >= 0 && f < s_nStrips) ? true : false); };
	bool ValidRatio(double r) { return ((r >= -1 && r <= 1) ? true : false); };
	void CalculateCorners();

	double m_centerPhi; //assuming det centered above x-axis (corresponds to zero phi)
	double m_centerZ;
	double m_centerRho; //perpendicular radius from axis

	std::vector<std::vector<ROOT::Math::XYZPoint>> m_frontStripCoords, m_backStripCoords;
	std::vector<std::vector<ROOT::Math::XYZPoint>> m_rotFrontStripCoords, m_rotBackStripCoords;

	ROOT::Math::XYZVector m_norm;

	ROOT::Math::RotationZ m_zRotation;

	std::uniform_real_distribution<double> m_uniformFraction;

	bool m_isSmearing;

	//Units in meters
	static constexpr double s_nStrips = 4; //Same for front and back
	static constexpr double s_nCorners = 4;
	static constexpr double s_totalLength = 0.075; //length of front strips
	static constexpr double s_backStripLength = s_totalLength / s_nStrips;
	static constexpr double s_totalWidth = 0.04; //width of back strips
	static constexpr double s_frontStripWidth = s_totalWidth / s_nStrips;
	static constexpr double s_deg2rad = M_PI/180.0;

};

#endif
