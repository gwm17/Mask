#ifndef __SABREDETECTOR_H
#define __SABREDETECTOR_H

#include <TRandom3.h>

struct CartCoords {

  double x;
  double y;
  double z;

  double operator[](int i) { //Overloaded for compatibility with Get_Cart. Only for access
    switch(i) {
      case(0): return this->x;
      case(1): return this->y;
      case(2): return this->z;
      default: return 0;
    }
  }

  double GetTheta();
  double GetR();
  double GetPhi();
  double GetDetectorPhi();
};

class SabreDetGeometry {

 public:

  SabreDetGeometry(double iRinner_flat,
		double iRouter_flat,
		double ideltaPhi_flat,
		double ibeamPhi_central,
		double itiltFromVertical,
    double idistFromTarget,
    double xoffset=0,
    double yoffset=0);
  ~SabreDetGeometry();

  double Get_Ring_Flat_X(int ch, int corner);
  double Get_Ring_Flat_Y(int ch, int corner);
  double Get_Ring_Flat_Z(int ch, int corner);
  double Get_Ring_Flat_R(int ch, int corner);
  double Get_Ring_Flat_Theta(int ch, int corner);
  double Get_Ring_Flat_Phi(int ch, int corner);

  double Get_Wedge_Flat_X(int ch, int corner);
  double Get_Wedge_Flat_Y(int ch, int corner);
  double Get_Wedge_Flat_Z(int ch, int corner);
  double Get_Wedge_Flat_R(int ch, int corner);
  double Get_Wedge_Flat_Theta(int ch, int corner);
  double Get_Wedge_Flat_Phi(int ch, int corner);

  double Get_Ring_Tilted_X(int ch, int corner);
  double Get_Ring_Tilted_Y(int ch, int corner);
  double Get_Ring_Tilted_Z(int ch, int corner);
  double Get_Ring_Tilted_R(int ch, int corner);
  double Get_Ring_Tilted_Theta(int ch, int corner);
  double Get_Ring_Tilted_Phi(int ch, int corner);

  double Get_Wedge_Tilted_X(int ch, int corner);
  double Get_Wedge_Tilted_Y(int ch, int corner);
  double Get_Wedge_Tilted_Z(int ch, int corner);
  double Get_Wedge_Tilted_R(int ch, int corner);
  double Get_Wedge_Tilted_Theta(int ch, int corner);
  double Get_Wedge_Tilted_Phi(int ch, int corner);

  bool IsInside(double theta, double phi);
  CartCoords GetTrajectoryCoordinates(double theta, double phi);

  void Recenter(double x, double y);


  int NumRings();
  int NumWedges();

  /*** Determine coordinates of the hit (ringch, wedgech) ***/
  CartCoords GetCoordinates(int ringch, int wedgech);
  double GetScatteringAngle(int ringch, int wedgech); //shortcut to Scattering angle

 private:

  const int NUMRINGS;
  const int NUMWEDGES;

  //detector corners
  int rbc; //ring bottom channel
  int rtc; //ring top channel
  int wrc; //wedge right channel
  int wlc; //wedge left channel

  double Rinner_flat;
  double Router_flat;
  double deltaR_flat;
  double deltaR_flat_ring;

  double deltaPhi_flat;
  double deltaPhi_flat_wedge;

  double beamPhi_central;
  double tiltFromVertical;
  double ZdistFromTarget;
  double XdistFromTarget;
  double YdistFromTarget;

  TRandom3* random;
  //default storage is cartesian
  //0=x, 1=y, 2=z
  CartCoords **ringch_flat_cart;
  CartCoords **wedgech_flat_cart;

  CartCoords **ringch_tilted_cart;
  CartCoords **wedgech_tilted_cart;

  double Get_Cart(int fot, int row, int ch, int corner, int cart);
  //fot = flat (0) or tilted (1)
  //row = ring (0) or wedge (1)

  bool CheckRingChannel(int);
  bool CheckWedgeChannel(int);
  bool CheckCorner(int);

  bool CheckBothRing(int,int);
  bool CheckBothWedge(int,int);

  /*** Perform transformation for arbitrary point on plane ***/
  double **XRotMatrix; //rotation about x-axis
  double **ZRotMatrix; //rotation about z-axis
  double **XRotMatrixInv; //inverse of x-rotation
  double **ZRotMatrixInv; //inverse of z-rotation
  CartCoords TransformVector(CartCoords vector);
  CartCoords InverseTransformVector(CartCoords vector);

};

#endif
