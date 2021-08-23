#ifndef __STRIP_DETECTOR_H
#define __STRIP_DETECTOR_H

// +z is along beam axis
// +y is vertically "upward" in the lab frame

//angles must be in radians, but distances can be whatever
//PROVIDED all input distances are the same

//Front strips from lowest y to highest y

//Back strips from lowest z to highest z

#include <TRandom3.h>
#include <cmath>
#include <vector>

#include "Vec3.h"
#include "Rotation.h"

class StripDetector {
public:
  
	StripDetector(int ns, double len, double wid, double cphi, double cz, double crho);
	~StripDetector();
	inline Mask::Vec3 GetFrontStripCoordinates(int stripch, int corner) { return front_strip_coords[stripch][corner]; };
	inline Mask::Vec3 GetBackStripCoordinates(int stripch, int corner) { return back_strip_coords[stripch][corner]; };
	inline Mask::Vec3 GetRotatedFrontStripCoordinates(int stripch, int corner) { return rotated_front_strip_coords[stripch][corner]; };
	inline Mask::Vec3 GetRotatedBackStripCoordinates(int stripch, int corner) { return rotated_back_strip_coords[stripch][corner]; };
	inline void SetRandomNumberGenerator(TRandom3* random) { m_random = random; };
	Mask::Vec3 GetHitCoordinates(int front_stripch, double front_strip_ratio);
	std::pair<int,double> GetChannelRatio(double theta, double phi);

private:
	void CalculateCorners();

	int num_strips;
	static constexpr int num_corners = 4;

	double length; //common to all strips, hence total
	double total_width;
	double front_strip_width; //assuming equal widths
	double back_strip_length; //assuming equal widths

	double center_phi; //assuming det centered above x-axis (corresponds to zero phi)
	double center_z;
	double center_rho; //perpendicular radius from axis

	std::vector<std::vector<Mask::Vec3>> front_strip_coords, back_strip_coords;
	std::vector<std::vector<Mask::Vec3>> rotated_front_strip_coords, rotated_back_strip_coords;

	Mask::ZRotation zRot;

	TRandom3* m_random; //Not owned by StripDetector!

	inline bool ValidChannel(int f) { return ((f >= 0 && f < num_strips) ? true : false); };
	inline bool ValidRatio(double r) { return ((r >= -1 && r <= 1) ? true : false); };
};

#endif
