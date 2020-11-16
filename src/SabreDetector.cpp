#include "SabreDetector.h"

#include <cmath>
#include <iostream>

using namespace std;

/*
  Distances in meters, angles in radians.

  The channel arrays have four points, one for each corner. The corners are
  as follows, as if looking BACK along beam (i.e. from the target's pov):

  0---------------------1
  |                     |
  |                     |      y
  |                     |      ^
  |                     |      |
  |                     |      |
  3---------------------2      -----> x
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

  ***note that these are for ARRAY storage, and may not necessarily
     correspond to the PHYSICAL channels; this will need to be taken
     into account if used in actual data analysis

  -- kgh, March 2020

*/

SabreDetGeometry::SabreDetGeometry(double iRinner_flat,
			     double iRouter_flat,
			     double ideltaPhi_flat,
			     double ibeamPhi_central,
			     double itiltFromVertical,
			     double idistFromTarget,
           double xoffset,
           double yoffset) :
  NUMRINGS(16), NUMWEDGES(8) {

  rbc = 0; //ring bottom channel
  rtc = NUMRINGS-1; //ring top channel
  wrc = 0; //wedge right channel
  wlc = NUMWEDGES-1; //wedge left channel

  Rinner_flat = iRinner_flat;
  Router_flat = iRouter_flat;
  deltaR_flat = iRouter_flat - iRinner_flat;
  deltaR_flat_ring = deltaR_flat/NUMRINGS;

  deltaPhi_flat = ideltaPhi_flat;
  deltaPhi_flat_wedge = deltaPhi_flat/NUMWEDGES;

  beamPhi_central = ibeamPhi_central;
  tiltFromVertical = itiltFromVertical;

  //this distance from target is from the CENTRAL point of
  //the detector, i.e. the center of the "circle" which the
  //detector forms. NOTE: these are all assumed to be negative!
  ZdistFromTarget = idistFromTarget;
  XdistFromTarget = xoffset;
  YdistFromTarget = yoffset;

  random = new TRandom3();
  random->SetSeed();

  ringch_flat_cart = new CartCoords*[NUMRINGS];
  wedgech_flat_cart = new CartCoords*[NUMWEDGES];
  ringch_tilted_cart = new CartCoords*[NUMRINGS];
  wedgech_tilted_cart = new CartCoords*[NUMWEDGES];
  for (int i=0; i<NUMRINGS; i++) {
    ringch_flat_cart[i] = new CartCoords[4];
    ringch_tilted_cart[i] = new CartCoords[4];
  }

  for (int i=0; i<NUMWEDGES; i++) {
    wedgech_flat_cart[i] = new CartCoords[4];
    wedgech_tilted_cart[i] = new CartCoords[4];
  }


  //Define rotation matricies
  XRotMatrix = new double*[3];
  ZRotMatrix = new double*[3];
  XRotMatrixInv = new double*[3];
  ZRotMatrixInv = new double*[3];
  for (int i=0; i<3; i++) {
    XRotMatrix[i] = new double[3];
    ZRotMatrix[i] = new double[3]; 
    XRotMatrixInv[i] = new double[3];
    ZRotMatrixInv[i] = new double[3];
  }

  XRotMatrix[0][0] = 1;
  XRotMatrix[0][1] = 0; XRotMatrix[0][2] = 0; XRotMatrix[1][0] = 0; XRotMatrix[2][0] = 0;
  XRotMatrix[1][1] = cos(tiltFromVertical); XRotMatrix[1][2] = -sin(tiltFromVertical);
  XRotMatrix[2][1] = sin(tiltFromVertical); XRotMatrix[2][2] = cos(tiltFromVertical);

  ZRotMatrix[0][0] = cos(beamPhi_central); ZRotMatrix[0][1] = -sin(beamPhi_central);
  ZRotMatrix[1][0] = sin(beamPhi_central); ZRotMatrix[1][1] = cos(beamPhi_central);
  ZRotMatrix[2][2] = 1;
  ZRotMatrix[2][0] = 0; ZRotMatrix[2][1] = 0; ZRotMatrix[0][2] = 0; ZRotMatrix[1][2] = 0;

  //Inverse of a rotation is its transpose

  XRotMatrixInv[0][0] = XRotMatrix[0][0]; XRotMatrixInv[0][1] = XRotMatrix[1][0]; XRotMatrixInv[0][2] = XRotMatrix[2][0];
  XRotMatrixInv[1][0] = XRotMatrix[0][1]; XRotMatrixInv[1][1] = XRotMatrix[1][1]; XRotMatrixInv[1][2] = XRotMatrix[2][1];
  XRotMatrixInv[2][0] = XRotMatrix[0][2]; XRotMatrixInv[2][1] = XRotMatrix[1][2]; XRotMatrixInv[2][2] = XRotMatrix[2][2];

  ZRotMatrixInv[0][0] = ZRotMatrix[0][0]; ZRotMatrixInv[0][1] = ZRotMatrix[1][0]; ZRotMatrixInv[0][2] = ZRotMatrix[2][0];
  ZRotMatrixInv[1][0] = ZRotMatrix[0][1]; ZRotMatrixInv[1][1] = ZRotMatrix[1][1]; ZRotMatrixInv[1][2] = ZRotMatrix[2][1];
  ZRotMatrixInv[2][0] = ZRotMatrix[0][2]; ZRotMatrixInv[2][1] = ZRotMatrix[1][2]; ZRotMatrixInv[2][2] = ZRotMatrix[2][2];
  //
 
  /* SETTING FLAT GEOMETRY */

  //RINGS:

  for (int rch=0; rch<NUMRINGS; rch++) {

    //top left:
    ringch_flat_cart[rch][0].x = (Rinner_flat + (rch+1)*deltaR_flat_ring)*sin(deltaPhi_flat/2); //x
    ringch_flat_cart[rch][0].y = (Rinner_flat + (rch+1)*deltaR_flat_ring)*cos(deltaPhi_flat/2); //y
    ringch_flat_cart[rch][0].z = 0; //z

    //top right:
    ringch_flat_cart[rch][1].x = (Rinner_flat + (rch+1)*deltaR_flat_ring)*sin(-deltaPhi_flat/2); //x
    ringch_flat_cart[rch][1].y = (Rinner_flat + (rch+1)*deltaR_flat_ring)*cos(-deltaPhi_flat/2); //y
    ringch_flat_cart[rch][1].z = 0; //z

    //bottom right:
    ringch_flat_cart[rch][2].x = (Rinner_flat + rch*deltaR_flat_ring)*sin(-deltaPhi_flat/2); //x
    ringch_flat_cart[rch][2].y = (Rinner_flat + rch*deltaR_flat_ring)*cos(-deltaPhi_flat/2); //y
    ringch_flat_cart[rch][2].z = 0; //z

    //bottom left:
    ringch_flat_cart[rch][3].x = (Rinner_flat + rch*deltaR_flat_ring)*sin(deltaPhi_flat/2); //x
    ringch_flat_cart[rch][3].y = (Rinner_flat + rch*deltaR_flat_ring)*cos(deltaPhi_flat/2); //y
    ringch_flat_cart[rch][3].z = 0; //z

  }

  //WEDGES:

  for (int wch=0; wch<NUMWEDGES; wch++) {

    //top left:
    wedgech_flat_cart[wch][0].x = Router_flat*sin(-deltaPhi_flat/2 + (wch+1)*deltaPhi_flat_wedge); //x
    wedgech_flat_cart[wch][0].y = Router_flat*cos(-deltaPhi_flat/2 + (wch+1)*deltaPhi_flat_wedge); //y
    wedgech_flat_cart[wch][0].z = 0; //z

    //top right:
    wedgech_flat_cart[wch][1].x = Router_flat*sin(-deltaPhi_flat/2 + wch*deltaPhi_flat_wedge); //x
    wedgech_flat_cart[wch][1].y = Router_flat*cos(-deltaPhi_flat/2 + wch*deltaPhi_flat_wedge); //y
    wedgech_flat_cart[wch][1].z = 0; //z

    //bottom right:
    wedgech_flat_cart[wch][2].x = Rinner_flat*sin(-deltaPhi_flat/2 + wch*deltaPhi_flat_wedge); //x
    wedgech_flat_cart[wch][2].y = Rinner_flat*cos(-deltaPhi_flat/2 + wch*deltaPhi_flat_wedge); //y
    wedgech_flat_cart[wch][2].z = 0; //z

    //bottom left:
    wedgech_flat_cart[wch][3].x = Rinner_flat*sin(-deltaPhi_flat/2 + (wch+1)*deltaPhi_flat_wedge); //x
    wedgech_flat_cart[wch][3].y = Rinner_flat*cos(-deltaPhi_flat/2 + (wch+1)*deltaPhi_flat_wedge); //y
    wedgech_flat_cart[wch][3].z = 0; //z

  }

  //now tilting, rotating, and translating

  //RINGS:

  for (int rch=0; rch<NUMRINGS; rch++) {
    for (int c=0; c<4; c++) {
      ringch_tilted_cart[rch][c] = TransformVector(ringch_flat_cart[rch][c]);
    }
  }

  //WEDGES:

  for (int wch=0; wch<NUMWEDGES; wch++) {
    for (int c=0; c<4; c++) {
      wedgech_tilted_cart[wch][c] = TransformVector(wedgech_flat_cart[wch][c]);
    }
  }

}

SabreDetGeometry::~SabreDetGeometry() {

  for (int i=0; i<NUMRINGS; i++) {
    delete [] ringch_flat_cart[i];
    delete [] ringch_tilted_cart[i];
  }
  delete [] ringch_flat_cart;
  delete [] ringch_tilted_cart;

  for (int i=0; i<NUMWEDGES; i++) {
    delete [] wedgech_flat_cart[i];
    delete [] wedgech_tilted_cart[i];
  } 
  delete [] wedgech_flat_cart;
  delete [] wedgech_tilted_cart;

  for(int i=0; i<3; i++) {
    delete [] XRotMatrix[i];
    delete [] ZRotMatrix[i];
    delete [] XRotMatrixInv[i];
    delete [] ZRotMatrixInv[i];
  }

  delete [] XRotMatrix;
  delete [] ZRotMatrix;

  delete random;

}

double SabreDetGeometry::Get_Ring_Flat_X(int ch, int corner) {
  return CheckBothRing(ch,corner) ? Get_Cart(0,0,ch,corner,0) : 0.;
}
double SabreDetGeometry::Get_Ring_Flat_Y(int ch, int corner) {
  return CheckBothRing(ch,corner) ? Get_Cart(0,0,ch,corner,1) : 0.;
}
double SabreDetGeometry::Get_Ring_Flat_Z(int ch, int corner) {
  return CheckBothRing(ch,corner) ? Get_Cart(0,0,ch,corner,2) : 0.;
}
double SabreDetGeometry::Get_Ring_Flat_R(int ch, int corner) {
  return CheckBothRing(ch,corner) ? 
    sqrt(pow(Get_Cart(0,0,ch,corner,0),2) +
	 pow(Get_Cart(0,0,ch,corner,1),2) +
	 pow(Get_Cart(0,0,ch,corner,2),2))
    : 0.;
}
double SabreDetGeometry::Get_Ring_Flat_Theta(int ch, int corner) {
  return CheckBothRing(ch,corner) ?
    acos(Get_Ring_Flat_Z(ch,corner)/Get_Ring_Flat_R(ch,corner))
    : 0.;
}
//GWM: Note this now returns in normal phi from x-axis format
double SabreDetGeometry::Get_Ring_Flat_Phi(int ch, int corner) {
  if (!(CheckBothRing(ch,corner))) return 0.;
  double x = Get_Ring_Flat_X(ch,corner);
  double y = Get_Ring_Flat_Y(ch,corner);
  double phi = std::atan2(y,x);
  if (phi<0) phi += M_PI*2.0;
  return phi;
}

double SabreDetGeometry::Get_Wedge_Flat_X(int ch, int corner) {
  return CheckBothWedge(ch,corner) ? Get_Cart(0,1,ch,corner,0) : 0.;
}
double SabreDetGeometry::Get_Wedge_Flat_Y(int ch, int corner) {
  return CheckBothWedge(ch,corner) ? Get_Cart(0,1,ch,corner,1) : 0.;
}
double SabreDetGeometry::Get_Wedge_Flat_Z(int ch, int corner) {
  return CheckBothWedge(ch,corner) ? Get_Cart(0,1,ch,corner,2) : 0.;
}
double SabreDetGeometry::Get_Wedge_Flat_R(int ch, int corner) {
  return CheckBothWedge(ch,corner) ?
    sqrt(pow(Get_Cart(0,1,ch,corner,0),2) +
	 pow(Get_Cart(0,1,ch,corner,1),2) +
	 pow(Get_Cart(0,1,ch,corner,2),2))
    : 0.;
}
double SabreDetGeometry::Get_Wedge_Flat_Theta(int ch, int corner) {
  return CheckBothWedge(ch,corner) ?
    acos(Get_Wedge_Flat_Z(ch,corner)/Get_Wedge_Flat_R(ch,corner))
    : 0.;
}
//GWM: Note this now returns in normal phi from x-axis format
double SabreDetGeometry::Get_Wedge_Flat_Phi(int ch, int corner) {
  if (!(CheckBothWedge(ch,corner))) return false;
  double x = Get_Wedge_Flat_X(ch,corner);
  double y = Get_Wedge_Flat_Y(ch,corner);
  double phi = std::atan2(y,x);
  if(phi<0) phi += M_PI*2.0;
  return phi;
}

double SabreDetGeometry::Get_Ring_Tilted_X(int ch, int corner) {
  return CheckBothRing(ch,corner) ? Get_Cart(1,0,ch,corner,0) : 0.;
}
double SabreDetGeometry::Get_Ring_Tilted_Y(int ch, int corner) {
  return CheckBothRing(ch,corner) ? Get_Cart(1,0,ch,corner,1) : 0.;
}
double SabreDetGeometry::Get_Ring_Tilted_Z(int ch, int corner) {
  return CheckBothRing(ch,corner) ? Get_Cart(1,0,ch,corner,2) : 0.;
}
double SabreDetGeometry::Get_Ring_Tilted_R(int ch, int corner) {
  return CheckBothRing(ch,corner) ?
    sqrt(pow(Get_Cart(1,0,ch,corner,0),2) +
	 pow(Get_Cart(1,0,ch,corner,1),2) +
	 pow(Get_Cart(1,0,ch,corner,2),2))
    : 0.;
}
double SabreDetGeometry::Get_Ring_Tilted_Theta(int ch, int corner) {
  return CheckBothRing(ch,corner) ?
    acos(Get_Ring_Tilted_Z(ch,corner)/Get_Ring_Tilted_R(ch,corner))
    : 0.;
}

//GWM: Note this now returns in normal phi from x-axis format
double SabreDetGeometry::Get_Ring_Tilted_Phi(int ch, int corner) {
  if (!(CheckBothRing(ch,corner))) return false;
  double x = Get_Ring_Tilted_X(ch,corner);
  double y = Get_Ring_Tilted_Y(ch,corner);
  double phi = std::atan2(y,x);
  if(phi<0) phi += M_PI*2.0;
  return phi;
}

double SabreDetGeometry::Get_Wedge_Tilted_X(int ch, int corner) {
  return CheckBothWedge(ch,corner) ? Get_Cart(1,1,ch,corner,0) : 0.;
}
double SabreDetGeometry::Get_Wedge_Tilted_Y(int ch, int corner) {
  return CheckBothWedge(ch,corner) ? Get_Cart(1,1,ch,corner,1) : 0.;
}
double SabreDetGeometry::Get_Wedge_Tilted_Z(int ch, int corner) {
  return CheckBothWedge(ch,corner) ? Get_Cart(1,1,ch,corner,2) : 0.;
}
double SabreDetGeometry::Get_Wedge_Tilted_R(int ch, int corner) {
  return CheckBothWedge(ch,corner) ?
    sqrt(pow(Get_Cart(1,1,ch,corner,0),2) +
	 pow(Get_Cart(1,1,ch,corner,1),2) +
	 pow(Get_Cart(1,1,ch,corner,2),2))
    : 0.;
}
double SabreDetGeometry::Get_Wedge_Tilted_Theta(int ch, int corner) {
  return CheckBothWedge(ch,corner) ?
    acos(Get_Wedge_Tilted_Z(ch,corner)/Get_Wedge_Tilted_R(ch,corner))
    : 0.;
}

//GWM: Note this now returns in normal phi from x-axis format
double SabreDetGeometry::Get_Wedge_Tilted_Phi(int ch, int corner) {
  if (!(CheckBothWedge(ch,corner))) return false;
  double x = Get_Wedge_Tilted_X(ch,corner);
  double y = Get_Wedge_Tilted_Y(ch,corner);
  double phi = std::atan2(y, x);
  if(phi<0) phi += 2.0*M_PI;
  return phi;
}

int SabreDetGeometry::NumRings()  {return NUMRINGS;}
int SabreDetGeometry::NumWedges() {return NUMWEDGES;}

double SabreDetGeometry::Get_Cart(int fot, int row, int ch, int corner, int cart) {

  //fot = flat (0) or tilted (1)
  //row = ring (0) or wedge (1)

  int marker = 10*fot + row;

  switch (marker) {
  case (00): return ringch_flat_cart[ch][corner][cart];
  case (01): return wedgech_flat_cart[ch][corner][cart];
  case (10): return ringch_tilted_cart[ch][corner][cart];
  case (11): return wedgech_tilted_cart[ch][corner][cart];
  default: return 0;
  }

}

bool SabreDetGeometry::CheckRingChannel(int ich) {
  return ((ich >= 0 && ich < NUMRINGS) ? true : false);
}

bool SabreDetGeometry::CheckWedgeChannel(int ich) {
  return ((ich >= 0 && ich < NUMWEDGES) ? true: false);
}

bool SabreDetGeometry::CheckCorner(int icn) {
  return ((icn >= 0 && icn < 4) ? true : false);
}

bool SabreDetGeometry::CheckBothRing(int ich, int icn) {
  return (CheckRingChannel(ich)&&CheckCorner(icn));
}

bool SabreDetGeometry::CheckBothWedge(int ich, int icn) {
  return (CheckWedgeChannel(ich)&&CheckCorner(icn));
}

/*Method by which the coordinates of a hit in a wedge/ring pixel are calculated.
 *Currently, takes the center of the pixel as the value of a pixel hit, could be altered
 *to take a uniformly sampled random point with in the pixel
 */
CartCoords SabreDetGeometry::GetCoordinates(int ringch, int wedgech) {
  if(!CheckRingChannel(ringch) || !CheckWedgeChannel(wedgech)) {
    CartCoords temp;
    temp.x = 0;
    temp.y = 0;
    temp.z = 0;
    return temp;
  }
  //define pixel by its center (half way between top and bottom radius, halfway between left and right phi)
  //EDIT GWM: July 2020 change to randomize the location within the sabre pixel
  double rcenter = Rinner_flat + deltaR_flat_ring*(ringch + random->Uniform(0.0, 1.0));
  double phi_center = deltaPhi_flat/2.0-(wedgech + random->Uniform(0.0, 1.0))*deltaPhi_flat_wedge;

  CartCoords PixelCenter;
  PixelCenter.x = rcenter*sin(phi_center);
  PixelCenter.y = rcenter*cos(phi_center);
  PixelCenter.z = 0;

  //return coords in final orientation
  return TransformVector(PixelCenter);
}

CartCoords SabreDetGeometry::TransformVector(CartCoords vector) {
  CartCoords xrot_vector, xzrot_vector, xzrot_t_vector;
  std::cout<<"Starting coords -- x: "<<vector.x<<" y: "<<vector.y<<" z: "<<vector.z<<" r: "<<vector.GetR()<<std::endl;
  //rotation about x-axis
  xrot_vector.x = vector.x;
  xrot_vector.y = XRotMatrix[1][1]*vector.y+XRotMatrix[1][2]*vector.z;
  xrot_vector.z = XRotMatrix[2][1]*vector.y+XRotMatrix[2][2]*vector.z;

  std::cout<<"after x-rotation -- x: "<<xrot_vector.x<<" y: "<<xrot_vector.y<<" z: "<<xrot_vector.z<<" r: "<<vector.GetR()<<std::endl;

  //rotation about z-axis
  xzrot_vector.x = ZRotMatrix[0][0]*xrot_vector.x+ZRotMatrix[0][1]*xrot_vector.y;
  xzrot_vector.y = ZRotMatrix[1][0]*xrot_vector.x+ZRotMatrix[1][1]*xrot_vector.y;
  xzrot_vector.z = ZRotMatrix[2][2]*xrot_vector.z;

  std::cout<<"after z-rotation -- x: "<<xzrot_vector.x<<" y: "<<xzrot_vector.y<<" z: "<<xzrot_vector.z<<" r: "<<vector.GetR()<<std::endl;

  //translate as needed
  xzrot_t_vector.x = xzrot_vector.x-XdistFromTarget;
  xzrot_t_vector.y = xzrot_vector.y-YdistFromTarget;
  xzrot_t_vector.z = xzrot_vector.z-ZdistFromTarget;

  std::cout<<"after z-translation -- x: "<<xzrot_t_vector.x<<" y: "<<xzrot_t_vector.y<<" z: "<<xzrot_t_vector.z<<std::endl;


  return xzrot_t_vector;
}

CartCoords SabreDetGeometry::InverseTransformVector(CartCoords vector) {
  CartCoords it_vector, izrot_vector, ixzrot_vector;
  //Inverse of translation
  /*it_vector.x = vector.x + XdistFromTarget;
  it_vector.y = vector.y + YdistFromTarget;
  it_vector.z = vector.z + ZdistFromTarget;*/

  //Inverse of rotation about z-axis
  izrot_vector.x = vector.x*ZRotMatrixInv[0][0] + vector.y*ZRotMatrixInv[0][1] + vector.z*ZRotMatrixInv[0][2];
  izrot_vector.y = vector.x*ZRotMatrixInv[1][0] + vector.y*ZRotMatrixInv[1][1] + vector.z*ZRotMatrixInv[1][2];
  izrot_vector.z = vector.x*ZRotMatrixInv[2][0] + vector.y*ZRotMatrixInv[2][1] + vector.z*ZRotMatrixInv[2][2];

  //Inverse of rotation about x-axis
  ixzrot_vector.x = izrot_vector.x*XRotMatrixInv[0][0] + izrot_vector.y*XRotMatrixInv[0][1] + izrot_vector.z*XRotMatrixInv[0][2];
  ixzrot_vector.y = izrot_vector.x*XRotMatrixInv[1][0] + izrot_vector.y*XRotMatrixInv[1][1] + izrot_vector.z*XRotMatrixInv[1][2];
  ixzrot_vector.z = izrot_vector.x*XRotMatrixInv[2][0] + izrot_vector.y*XRotMatrixInv[2][1] + izrot_vector.z*XRotMatrixInv[2][2];

  return ixzrot_vector;

}

void SabreDetGeometry::Recenter(double x, double y) { 
  XdistFromTarget = x;
  YdistFromTarget = y;
}

double SabreDetGeometry::GetScatteringAngle(int ringch, int wedgech) {
  if(!CheckRingChannel(ringch) || !CheckWedgeChannel(wedgech)) return 0;
  CartCoords these_coords = GetCoordinates(ringch, wedgech);
  double r = sqrt(pow(these_coords.x, 2.) + pow(these_coords.y, 2.) + pow(these_coords.z, 2.));
  return acos(these_coords.z/r);
}

/*
  
*/
bool SabreDetGeometry::IsInside(double theta, double phi) {
  //Make unit vector
  CartCoords unit, origin;
  origin.x = 0; origin.y = 0; origin.z = ZdistFromTarget;
  origin = InverseTransformVector(origin);

  unit.x = sin(theta)*cos(phi); unit.y = sin(theta)*sin(phi); unit.z = cos(theta);

  unit = InverseTransformVector(unit);

  double s = origin.y/unit.y;
  double px = origin.x + s*unit.x;
  double py = origin.z + s*unit.z;

  double radius = std::sqrt(px*px + py*py);
  std::cout<<"radius: "<<radius<<std::endl;



  double phi_in_detector = unit.GetPhi();
  double theta_in_detector = unit.GetTheta();
  double xy_radius = std::sqrt(unit.x*unit.x + unit.y*unit.y);
  double xz_radius = std::sqrt(unit.x*unit.x + unit.z*unit.z);

  double val = std::atan2(py, px);
  if(val < 0 ) val += M_PI*2.0;
  std::cout<<"val: "<<val*180.0/M_PI<<std::endl;


  std::cout<<"phi in det: "<<phi_in_detector*180.0/M_PI<<std::endl;
  std::cout<<"theta in det: "<<theta_in_detector*180.0/M_PI<<std::endl;
  std::cout<<"xy radius: "<<xy_radius<<std::endl;
  std::cout<<"xz radius: "<<xz_radius<<std::endl;

  std::cout<<"unit -- x: "<<unit.x<<" y: "<<unit.y<<" z: "<<unit.z<<std::endl;

  std::cout<<" phi range tilt -- min: "<<Get_Wedge_Tilted_Phi(0,1)*180.0/M_PI<<" max: "<<Get_Wedge_Tilted_Phi(7,0)*180.0/M_PI<<std::endl;
  std::cout<<" phi range flat -- min: "<<Get_Wedge_Flat_Phi(0,1)*180.0/M_PI<<" max: "<<Get_Wedge_Flat_Phi(7,0)*180.0/M_PI<<std::endl;
  std::cout<<" wedge 0 1 theta -- "<<Get_Wedge_Tilted_Theta(0,1)*180.0/M_PI<<std::endl;
  std::cout<<" ring 0 0 x: "<<Get_Ring_Flat_X(0,0)<<std::endl;
  std::cout<<" ring 0 1 x: "<<Get_Ring_Flat_X(0,1)<<std::endl;
  std::cout<<" ring 0 0 y: "<<Get_Ring_Flat_Y(0,0)<<std::endl;
  std::cout<<" ring 0 1 y: "<<Get_Ring_Flat_Y(0,1)<<std::endl;

  std::cout<<" t ring 0 0 x: "<<Get_Ring_Tilted_X(0,0)<<std::endl;
  std::cout<<" t ring 0 1 x: "<<Get_Ring_Tilted_X(0,1)<<std::endl;
  std::cout<<" t ring 0 0 y: "<<Get_Ring_Tilted_Y(0,0)<<std::endl;
  std::cout<<" t ring 0 1 y: "<<Get_Ring_Tilted_Y(0,1)<<std::endl;

  std::cout<<"Test inversion"<<std::endl;
  CartCoords test = InverseTransformVector(wedgech_tilted_cart[0][1]);
  std::cout<<"Inverted x: "<<test.x<<" y: "<<test.y<<" z: "<<test.z<<std::endl;
  std::cout<<"Given x: "<<wedgech_flat_cart[0][1].x<<" y: "<<wedgech_flat_cart[0][1].y<<" z: "<<wedgech_flat_cart[0][1].z<<std::endl;

  //Check phi condition
  //double phi_in_detector = phi - beamPhi_central;
  double phi_min = Get_Wedge_Tilted_Phi(0,1);
  double phi_max = Get_Wedge_Tilted_Phi(7,0);
  bool passed;

  if(phi > phi_min && phi < phi_max) {
    passed = true;
  } else {
    return false;
  }

  std::cout<<" passed phi "<<std::endl;

  if(xy_radius>Router_flat || xy_radius<Rinner_flat) passed = false;

  return passed;
  
}

CartCoords SabreDetGeometry::GetTrajectoryCoordinates(double theta, double phi) {
  CartCoords temp;
  temp.x = 0; temp.y = 0; temp.z = 0;

  int ringchan = -1, wedgechan = -1;
  double phi_in_detector = phi - beamPhi_central;
  double phi_min, phi_max;
  for(int i=0; i<NUMWEDGES; i++) {
    phi_min = wedgech_flat_cart[i][0].GetPhi();
    phi_max = wedgech_flat_cart[i][1].GetPhi();
    if(phi_in_detector > phi_min && phi_in_detector < phi_max) {
      wedgechan = i;
      break;
    }
  }
  if(wedgechan == -1) return temp;

  CartCoords innerPosition_flat, innerPosition_tilt;
  CartCoords outerPosition_flat, outerPosition_tilt;
  double r_inner, r_outer;
  double theta_min, theta_max;
  //Since we passed phi, use the phi coordinate to calculate the point at the "top" and
  //bottom of the detector at the given phi, then transform to the lampshade frame.
  for(int i=0; i<NUMRINGS; i++) {
    r_inner = Rinner_flat + deltaR_flat_ring*(i);
    r_outer = Rinner_flat + deltaR_flat_ring*(i+1);

    innerPosition_flat.x = r_inner*sin(phi_in_detector);
    innerPosition_flat.y = r_inner*cos(phi_in_detector);
    innerPosition_flat.z = 0;
    innerPosition_tilt = TransformVector(innerPosition_flat);

    outerPosition_flat.x = r_outer*sin(phi_in_detector);
    outerPosition_flat.y = r_outer*cos(phi_in_detector);
    outerPosition_flat.z = 0;
    outerPosition_tilt = TransformVector(outerPosition_flat);

    //Get the theta values
    theta_max = innerPosition_tilt.GetTheta();
    theta_min = outerPosition_tilt.GetTheta();

    if(theta<=theta_max && theta>=theta_min) {
      ringchan = i;
      break;
    }
  }
  if(ringchan == -1) return temp;

  return GetCoordinates(ringchan, wedgechan);
}

//Coordinate functions
double CartCoords::GetTheta() {
  double r = std::sqrt(std::pow(x, 2.) + std::pow(y, 2.) + std::pow(z, 2.));
  return std::acos(z/r);
}

double CartCoords::GetR() {
  return std::sqrt(std::pow(x, 2.) + std::pow(y, 2.) + std::pow(z, 2.));
}

double CartCoords::GetPhi() {
  /*double r = std::sqrt(std::pow(x, 2.) + std::pow(y, 2.));
  if((x >= 0 && y >= 0) || (x <= 0 && y >= 0)) return std::acos(x/r);
  else if((x <= 0 && y <= 0) || (x >= 0 && y <= 0)) return (2.0*M_PI - std::acos(x/r));
  else return 0.0;*/
  double phi = std::atan2(y, x);
  if(phi<0) phi += M_PI*2.0;
  return phi;
}

double CartCoords::GetDetectorPhi() {
  /*double r = std::sqrt(std::pow(x, 2.) + std::pow(y, 2.));
  if((x >= 0 && y >= 0) || (x <= 0 && y >= 0)) return std::acos(x/r);
  else if((x <= 0 && y <= 0) || (x >= 0 && y <= 0)) return (2.0*M_PI - std::acos(x/r));
  else return 0.0;*/
  double phi = std::atan2(x, y);
  if(phi<0) phi += M_PI*2.0;
  return phi;
}

