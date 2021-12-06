#include "StripDetector.h"
#include <iostream>

/*
  Corner layout for each strip in the un-rotated frame
  0--------------------------1
  |                          |              
  |                          |              
  |                          |              
  |                          |              
  |                          |              x
  2--------------------------3    z<--------X
																						|
																						|
																						|
																						y
*/

StripDetector::StripDetector(int ns, double len, double wid, double cphi, double cz, double crho) :
	m_norm_unrot(1.0,0.0,0.0), m_uniform_fraction(0.0, 1.0), rndmFlag(false)
{

	num_strips = ns;
	
	length = std::fabs(len);
	total_width = std::fabs(wid);
	front_strip_width = total_width/num_strips;
	back_strip_length = length/num_strips;

	while (cphi < 0) cphi += 2*M_PI;
	center_phi = cphi;
	center_z   = cz;
	center_rho = std::fabs(crho);

	zRot.SetAngle(center_phi);

	front_strip_coords.resize(num_strips);
	back_strip_coords.resize(num_strips);
	rotated_front_strip_coords.resize(num_strips);
	rotated_back_strip_coords.resize(num_strips);
	for(int i=0; i<num_strips; i++) {
		front_strip_coords[i].resize(num_corners);
		back_strip_coords[i].resize(num_corners);
		rotated_front_strip_coords[i].resize(num_corners);
		rotated_back_strip_coords[i].resize(num_corners);
	}
	CalculateCorners();
	
}

StripDetector::~StripDetector() {}

void StripDetector::CalculateCorners() {
	double y_min, y_max, z_min, z_max; 
	for (int s=0; s<num_strips; s++) {
		y_max = total_width/2.0 - front_strip_width*s;
		y_min = total_width/2.0 - front_strip_width*(s+1);
		z_max = center_z + length/2.0;
		z_min = center_z - length/2.0;
		front_strip_coords[s][2] = Mask::Vec3(center_rho, y_max, z_max);
		front_strip_coords[s][3] = Mask::Vec3(center_rho, y_max, z_min);
		front_strip_coords[s][0] = Mask::Vec3(center_rho, y_min, z_max);
		front_strip_coords[s][1] = Mask::Vec3(center_rho, y_min, z_min);

		z_max = (center_z - length/2.0) + (s+1)*back_strip_length;
		z_min = (center_z - length/2.0) + (s)*back_strip_length;
		y_max = total_width/2.0;
		y_min = -total_width/2.0;
		back_strip_coords[s][2] = Mask::Vec3(center_rho, y_max, z_max);
		back_strip_coords[s][3] = Mask::Vec3(center_rho, y_max, z_min);
		back_strip_coords[s][0] = Mask::Vec3(center_rho, y_min, z_max);
		back_strip_coords[s][1] = Mask::Vec3(center_rho, y_min, z_min);
	}

	for(int s=0; s<num_strips; s++) {
		rotated_front_strip_coords[s][0] = zRot*front_strip_coords[s][0];
		rotated_front_strip_coords[s][1] = zRot*front_strip_coords[s][1];
		rotated_front_strip_coords[s][2] = zRot*front_strip_coords[s][2];
		rotated_front_strip_coords[s][3] = zRot*front_strip_coords[s][3];

		rotated_back_strip_coords[s][0] = zRot*back_strip_coords[s][0];
		rotated_back_strip_coords[s][1] = zRot*back_strip_coords[s][1];
		rotated_back_strip_coords[s][2] = zRot*back_strip_coords[s][2];
		rotated_back_strip_coords[s][3] = zRot*back_strip_coords[s][3];
	}
}

Mask::Vec3 StripDetector::GetHitCoordinates(int front_stripch, double front_strip_ratio) {

	if (!ValidChannel(front_stripch) || !ValidRatio(front_strip_ratio)) return Mask::Vec3(0,0,0);

	double y;
	if(rndmFlag) {
		y = -total_width/2.0 + (front_stripch + m_uniform_fraction(Mask::RandomGenerator::GetInstance().GetGenerator()))*front_strip_width;
	} else {
		y = -total_width/2.0 + (front_stripch+0.5)*front_strip_width;
	}

	//recall we're still assuming phi=0 det:
	Mask::Vec3 coords(center_rho, y, front_strip_ratio*(length/2) + center_z);
  
	//NOW rotate by appropriate phi
	return zRot*coords;

}

StripHit StripDetector::GetChannelRatio(double theta, double phi) {

	StripHit hit;
	while (phi < 0) phi += 2*M_PI;

	//to make the math easier (and the same for each det), rotate the input phi
	//BACKWARD by the phi of the det, s.t. we are looking along x-axis
	phi -= center_phi;

	if (phi > M_PI) phi -= 2*M_PI;

	//then we can check easily whether it even hit the detector in phi
	double det_max_phi = atan2(total_width/2,center_rho);
	double det_min_phi = -det_max_phi;
  
	if (phi < det_min_phi || phi > det_max_phi) return hit;

	//for theta it's not so simple, so we have to go through the typical plane-intersect method
	//first thing's first: we have a fixed x for the entire detector plane:
	double xhit = center_rho;
	//thus we find the corresponding y and z for that fixed x, given the input theta and phi:
	double yhit = xhit*tan(phi);
	double zhit = sqrt(xhit*xhit+yhit*yhit)/tan(theta);

	for (int s=0; s<num_strips; s++) {
		if (xhit >= front_strip_coords[s][0].GetX() && xhit <= front_strip_coords[s][0].GetX() && //Check min and max x (constant in flat)
			yhit >= front_strip_coords[s][1].GetY() && yhit <= front_strip_coords[s][2].GetY() && //Check min and max y
			zhit >= front_strip_coords[s][1].GetZ() && zhit <= front_strip_coords[s][0].GetZ()) //Check min and max z
		{
			hit.front_strip_index = s;
			hit.front_ratio = (zhit-center_z)/(length/2);
			break;
		}
	}

	for (int s=0; s<num_strips; s++) {
		if (xhit >= back_strip_coords[s][0].GetX() && xhit <= back_strip_coords[s][0].GetX() && //Check min and max x (constant in flat)
			yhit >= back_strip_coords[s][1].GetY() && yhit <= back_strip_coords[s][2].GetY() && //Check min and max y
			zhit >= back_strip_coords[s][1].GetZ() && zhit <= back_strip_coords[s][0].GetZ()) //Check min and max z
		{
			hit.back_strip_index = s;
			break;
		}
	}

	return hit;
}