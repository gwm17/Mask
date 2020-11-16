#include <iostream>
#include <cmath>
#include <fstream>

#include <TGraph.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraph2D.h>
#include <TLatex.h>
#include <TRandom3.h>

#include "SabreDetector.h"

using namespace std;

int main() {

  const int NUMDETS = 5;

  SabreDetGeometry **sabdet = new SabreDetGeometry*[NUMDETS];

  //dimensions common to all dets:
  const double INNER_RADIUS = 0.0326; //meters
  const double OUTER_RADIUS = 0.1351; //meters
  const double TILT_FROM_VERTICAL = 40.0; //degrees
  const double DIST_FROM_TARGET = 0.1245; //meters

  sabdet[0] = new SabreDetGeometry(INNER_RADIUS,
				OUTER_RADIUS,
				54.4*M_PI/180., //delta phi
				36.*M_PI/180.,  //center phi
				TILT_FROM_VERTICAL*M_PI/180.,
				DIST_FROM_TARGET);
  sabdet[1] = new SabreDetGeometry(INNER_RADIUS,
				OUTER_RADIUS,
				54.4*M_PI/180.,
				108.*M_PI/180.,
				TILT_FROM_VERTICAL*M_PI/180.,
				DIST_FROM_TARGET);
  sabdet[2] = new SabreDetGeometry(INNER_RADIUS,
				OUTER_RADIUS,
				54.4*M_PI/180.,
				180.*M_PI/180.,
				TILT_FROM_VERTICAL*M_PI/180.,
				DIST_FROM_TARGET);
  sabdet[3] = new SabreDetGeometry(INNER_RADIUS,
				OUTER_RADIUS,
				54.4*M_PI/180.,
				252.*M_PI/180.,
				TILT_FROM_VERTICAL*M_PI/180.,
				DIST_FROM_TARGET);
  sabdet[4] = new SabreDetGeometry(INNER_RADIUS,
				OUTER_RADIUS,
				54.4*M_PI/180.,
				324.*M_PI/180.,
				TILT_FROM_VERTICAL*M_PI/180.,
				DIST_FROM_TARGET);

  const int NUMCORNERS = 4;
  const int NUMRINGS   = sabdet[0]->NumRings();
  const int NUMWEDGES  = sabdet[0]->NumWedges();

  const int NUMRINGPTS = NUMDETS*NUMRINGS*NUMCORNERS;
  const int NUMWEDGEPTS = NUMDETS*NUMWEDGES*NUMCORNERS;
  
  //const int NUMMCCOUNTS = 1E8;

  double ring_xpts_flat[NUMRINGPTS];
  double ring_ypts_flat[NUMRINGPTS];
  double ring_zpts_flat[NUMRINGPTS];
  double ring_rpts_flat[NUMRINGPTS];
  double ring_tpts_flat[NUMRINGPTS];
  double ring_ppts_flat[NUMRINGPTS];
  double wedge_xpts_flat[NUMWEDGEPTS];
  double wedge_ypts_flat[NUMWEDGEPTS];
  double wedge_zpts_flat[NUMWEDGEPTS];
  double wedge_rpts_flat[NUMWEDGEPTS];
  double wedge_tpts_flat[NUMWEDGEPTS];
  double wedge_ppts_flat[NUMWEDGEPTS];

  double ring_xpts_tilted[NUMRINGPTS];
  double ring_ypts_tilted[NUMRINGPTS];
  double ring_zpts_tilted[NUMRINGPTS];
  double ring_rpts_tilted[NUMRINGPTS];
  double ring_tpts_tilted[NUMRINGPTS];
  double ring_ppts_tilted[NUMRINGPTS];
  double wedge_xpts_tilted[NUMWEDGEPTS];
  double wedge_ypts_tilted[NUMWEDGEPTS];
  double wedge_zpts_tilted[NUMWEDGEPTS];
  double wedge_rpts_tilted[NUMWEDGEPTS];
  double wedge_tpts_tilted[NUMWEDGEPTS];
  double wedge_ppts_tilted[NUMWEDGEPTS];


  ofstream checkfile;
  checkfile.open("checkfile.txt");

  checkfile << "RINGS\n"
	    << "det   ch   corner   x_fl   y_fl   z_fl   r_fl   t_fl   p_fl   x_ti   y_ti   z_ti   r_ti   t_ti   p_ti\n";
  for (int d=0; d<NUMDETS; d++) {
    for (int i=0; i<NUMRINGS; i++) {
      for (int j=0; j<NUMCORNERS; j++) {
	int index = NUMCORNERS*NUMRINGS*d + NUMCORNERS*i + j;
	ring_xpts_flat[index] = sabdet[d]->Get_Ring_Flat_X(i,j)*100;
	ring_ypts_flat[index] = sabdet[d]->Get_Ring_Flat_Y(i,j)*100;
	ring_zpts_flat[index] = sabdet[d]->Get_Ring_Flat_Z(i,j)*100;
	ring_rpts_flat[index] = sabdet[d]->Get_Ring_Flat_R(i,j)*100;
	ring_tpts_flat[index] = sabdet[d]->Get_Ring_Flat_Theta(i,j)*180./M_PI;
	ring_ppts_flat[index] = sabdet[d]->Get_Ring_Flat_Phi(i,j)*180./M_PI;
	ring_xpts_tilted[index] = sabdet[d]->Get_Ring_Tilted_X(i,j)*100;
	ring_ypts_tilted[index] = sabdet[d]->Get_Ring_Tilted_Y(i,j)*100;
	ring_zpts_tilted[index] = sabdet[d]->Get_Ring_Tilted_Z(i,j)*100;
	ring_rpts_tilted[index] = sabdet[d]->Get_Ring_Tilted_R(i,j)*100;
	ring_tpts_tilted[index] = sabdet[d]->Get_Ring_Tilted_Theta(i,j)*180./M_PI;
	ring_ppts_tilted[index] = sabdet[d]->Get_Ring_Tilted_Phi(i,j)*180./M_PI;

	checkfile << d << "   " << i << "    " << j << "        "
		  << ring_xpts_flat[index] << " "
		  << ring_ypts_flat[index] << " "
		  << ring_zpts_flat[index] << " "
		  << ring_rpts_flat[index] << " "
		  << ring_tpts_flat[index] << " "
		  << ring_ppts_flat[index] << " "
		  << ring_xpts_tilted[index] << " "
		  << ring_ypts_tilted[index] << " "
		  << ring_zpts_tilted[index] << " "
		  << ring_rpts_tilted[index] << " "
		  << ring_tpts_tilted[index] << " "
		  << ring_ppts_tilted[index] << endl;
      }
    }
  }

  checkfile << "WEDGES\n"
	    << "det   ch   corner   x_fl   y_fl   z_fl   r_fl   t_fl   p_fl   x_ti   y_ti   z_ti   r_ti   t_ti   p_ti\n";
  for (int d=0; d<NUMDETS; d++) {
    for (int i=0; i<NUMWEDGES; i++) {
      for (int j=0; j<NUMCORNERS; j++) {
	int index = NUMCORNERS*NUMWEDGES*d + NUMCORNERS*i + j;
	wedge_xpts_flat[index] = sabdet[d]->Get_Wedge_Flat_X(i,j)*100;
	wedge_ypts_flat[index] = sabdet[d]->Get_Wedge_Flat_Y(i,j)*100;
	wedge_zpts_flat[index] = sabdet[d]->Get_Wedge_Flat_Z(i,j)*100;
	wedge_rpts_flat[index] = sabdet[d]->Get_Wedge_Flat_R(i,j)*100;
	wedge_tpts_flat[index] = sabdet[d]->Get_Wedge_Flat_Theta(i,j)*180./M_PI;
	wedge_ppts_flat[index] = sabdet[d]->Get_Wedge_Flat_Phi(i,j)*180./M_PI;
	wedge_xpts_tilted[index] = sabdet[d]->Get_Wedge_Tilted_X(i,j)*100;
	wedge_ypts_tilted[index] = sabdet[d]->Get_Wedge_Tilted_Y(i,j)*100;
	wedge_zpts_tilted[index] = sabdet[d]->Get_Wedge_Tilted_Z(i,j)*100;
	wedge_rpts_tilted[index] = sabdet[d]->Get_Wedge_Tilted_R(i,j)*100;
	wedge_tpts_tilted[index] = sabdet[d]->Get_Wedge_Tilted_Theta(i,j)*180./M_PI;
	wedge_ppts_tilted[index] = sabdet[d]->Get_Wedge_Tilted_Phi(i,j)*180./M_PI;

	checkfile << d << "   " << i << "    " << j << "        "
		  << wedge_xpts_flat[index] << " "
		  << wedge_ypts_flat[index] << " "
		  << wedge_zpts_flat[index] << " "
		  << wedge_rpts_flat[index] << " "
		  << wedge_tpts_flat[index] << " "
		  << wedge_ppts_flat[index] << " "
		  << wedge_xpts_tilted[index] << " "
		  << wedge_ypts_tilted[index] << " "
		  << wedge_zpts_tilted[index] << " "
		  << wedge_rpts_tilted[index] << " "
		  << wedge_tpts_tilted[index] << " "
		  << wedge_ppts_tilted[index] << endl;
      }
    }
  }

  checkfile.close();

  /*
    UNDER CONSTRUCTION: MC efficiency integration

  TRandom3 *rand = new TRandom3();
  rand->SetSeed();
  int hits=0;
  bool didhit=false;
  double randx=0;
  double randy=0;
  double randz=0;
  double randr=0;
  double randt=0;
  double randp=0;
  for (int mc=0; mc<NUMMCCOUNTS; mc++) {
    didhit = false;
    randx = rand->Uniform(0,1);
    randy = rand->Uniform(0,1);
    randz = rand->Uniform(0,1);
    randr = sqrt(pow(randx,2)+pow(randy,2)+pow(randz,2));
    randt = acos(randz/randr)*180./M_PI;
    randp = atan(randy/randx)*180./M_PI;
    if (randx <= 0 && randy >= 0)      randp += 90;
    else if (randx <= 0 && randy <= 0) randp += 90;
    else if (randx >= 0 && randy <= 0) randp += 270;
    else if (randx >= 0 && randy >= 0) randp += 270;
    for (int d=0; d<NUMDETS; d++) {
      if (didhit) break;
      for (int f=0; f<NUMRINGS; f++) {
	if (didhit) break;
	for (int b=0; b<NUMWEDGES; b++) {
	  if (randt <= sabdet[d]->Get_Ring_Tilted_Theta(f,)) {
	    hit++;
	    didhit=true;
	    break;
	  }
	}
      }
    }
  }
  cout << "Geom eff = " << 1.*hits/NUMMCCOUNTS*100  << "%" << endl;
  */
  TGraph ring_flat_xy_gr(NUMRINGPTS,ring_xpts_flat,ring_ypts_flat);
  ring_flat_xy_gr.SetMarkerStyle(2);
  ring_flat_xy_gr.SetName("ring_flat_xy_gr");
  TGraph wedge_flat_xy_gr(NUMWEDGEPTS,wedge_xpts_flat,wedge_ypts_flat);
  wedge_flat_xy_gr.SetMarkerStyle(2);
  wedge_flat_xy_gr.SetName("wedge_flat_xy_gr");
  wedge_flat_xy_gr.SetTitle("Flat (Untilted);x (cm);y (cm)");
  TGraph ring_flat_zy_gr(NUMRINGPTS,ring_zpts_flat,ring_ypts_flat);
  ring_flat_zy_gr.SetMarkerStyle(2);
  ring_flat_zy_gr.SetName("ring_flat_zy_gr");
  TGraph wedge_flat_zy_gr(NUMWEDGEPTS,wedge_zpts_flat,wedge_ypts_flat);
  wedge_flat_zy_gr.SetMarkerStyle(2);
  wedge_flat_zy_gr.SetName("wedge_flat_zy_gr");
  wedge_flat_zy_gr.SetTitle("Flat (Untilted);z (cm);y (cm)");
  TGraph ring_flat_pt_gr(NUMRINGPTS,ring_ppts_flat,ring_tpts_flat);
  ring_flat_pt_gr.SetMarkerStyle(2);
  ring_flat_pt_gr.SetName("ring_flat_pt_gr");
  TGraph wedge_flat_pt_gr(NUMWEDGEPTS,wedge_ppts_flat,wedge_tpts_flat);
  wedge_flat_pt_gr.SetMarkerStyle(2);
  wedge_flat_pt_gr.SetName("wedge_flat_pt_gr");
  wedge_flat_pt_gr.SetTitle("Flat (Untilted);#phi (deg.);#theta (deg.)");
  TGraph2D ring_flat_2d_gr(NUMRINGPTS,ring_xpts_flat,ring_zpts_flat,ring_ypts_flat);
  ring_flat_2d_gr.SetMarkerStyle(2);
  ring_flat_2d_gr.SetName("ring_flat_2d_gr");
  TGraph2D wedge_flat_2d_gr(NUMWEDGEPTS,wedge_xpts_flat,wedge_zpts_flat,wedge_ypts_flat);
  wedge_flat_2d_gr.SetMarkerStyle(2);
  wedge_flat_2d_gr.SetName("wedge_flat_2d_gr");
  wedge_flat_2d_gr.SetTitle("Flat (Untilted);x (cm); z (cm); y (cm)");

  TGraph ring_tilted_xy_gr(NUMRINGPTS,ring_xpts_tilted,ring_ypts_tilted);
  ring_tilted_xy_gr.SetMarkerStyle(2);
  ring_tilted_xy_gr.SetName("ring_tilted_xy_gr");
  TGraph wedge_tilted_xy_gr(NUMWEDGEPTS,wedge_xpts_tilted,wedge_ypts_tilted);
  wedge_tilted_xy_gr.SetMarkerStyle(2);
  wedge_tilted_xy_gr.SetName("wedge_tilted_xy_gr");
  wedge_tilted_xy_gr.SetTitle("Tilted;x (cm);y (cm)");
  TGraph ring_tilted_zy_gr(NUMRINGPTS,ring_zpts_tilted,ring_ypts_tilted);
  ring_tilted_zy_gr.SetMarkerStyle(2);
  ring_tilted_zy_gr.SetName("ring_tilted_zy_gr");
  TGraph wedge_tilted_zy_gr(NUMWEDGEPTS,wedge_zpts_tilted,wedge_ypts_tilted);
  wedge_tilted_zy_gr.SetMarkerStyle(2);
  wedge_tilted_zy_gr.SetName("wedge_tilted_zy_gr");
  wedge_tilted_zy_gr.SetTitle("Tilted;z (cm);y (cm)");
  TGraph ring_tilted_pt_gr(NUMRINGPTS,ring_ppts_tilted,ring_tpts_tilted);
  ring_tilted_pt_gr.SetMarkerStyle(2);
  ring_tilted_pt_gr.SetName("ring_tilted_pt_gr");
  TGraph wedge_tilted_pt_gr(NUMWEDGEPTS,wedge_ppts_tilted,wedge_tpts_tilted);
  wedge_tilted_pt_gr.SetMarkerStyle(2);
  wedge_tilted_pt_gr.SetName("wedge_tilted_pt_gr");
  wedge_tilted_pt_gr.SetTitle("Tilted;#phi (deg.);#theta (deg.)");
  TGraph2D ring_tilted_2d_gr(NUMRINGPTS,ring_xpts_tilted,ring_zpts_tilted,ring_ypts_tilted);
  ring_tilted_2d_gr.SetMarkerStyle(2);
  ring_tilted_2d_gr.SetName("ring_tilted_2d_gr");
  TGraph2D wedge_tilted_2d_gr(NUMWEDGEPTS,wedge_xpts_tilted,wedge_zpts_tilted,wedge_ypts_tilted);
  wedge_tilted_2d_gr.SetMarkerStyle(2);
  wedge_tilted_2d_gr.SetName("wedge_tilted_2d_gr");
  wedge_tilted_2d_gr.SetTitle("Tilted;x (cm); z (cm); y (cm)");
 
  int nCenters = 5*16*8;
  int nCent_det = 16*8;
  double centerx[nCenters], centery[nCenters], centerz[nCenters];
  double **det_centx, **det_centy, **det_centz;
  det_centx = new double*[5]; det_centy = new double*[5]; det_centz = new double*[5];
  for(int i=0; i<5; i++) {
    det_centx[i] = new double[nCent_det];
    det_centy[i] = new double[nCent_det];
    det_centz[i] = new double[nCent_det];
  }
  CartCoords point;
  int index = 0;
  int small_index;
  for(int d=0; d<5; d++) {
    small_index = 0;
    for(int r=0; r<16; r++) {
      for(int w=0; w<8; w++) {
        point = sabdet[d]->GetCoordinates(r,w);
        centerx[index] = point.x*100;   
        centery[index] = point.y*100;   
        centerz[index] = point.z*100; 

        det_centx[d][small_index] = point.x*100;  
        det_centy[d][small_index] = point.y*100;  
        det_centz[d][small_index] = point.z*100;  
        index++; small_index++;
      }
    }
  }
  TGraph2D centerpoint(nCenters, &centerx[0], &centerz[0], &centery[0]);
  centerpoint.SetName("centerpoint");
  centerpoint.SetTitle("centerpoint;x;z;y");
  centerpoint.SetMarkerStyle(2);
  centerpoint.SetMarkerColor(2);

  TGraph2D centerpoint_d1(nCent_det, &(det_centx[0][0]), &(det_centz[0][0]), &(det_centy[0][0]));
  centerpoint_d1.SetName("centerpoint_d1");
  centerpoint_d1.SetMarkerStyle(2);
  centerpoint_d1.SetMarkerColor(1);
  TGraph2D centerpoint_d2(nCent_det, &(det_centx[1][0]), &(det_centz[1][0]), &(det_centy[1][0]));
  centerpoint_d2.SetName("centerpoint_d2");
  centerpoint_d2.SetMarkerStyle(2);
  centerpoint_d2.SetMarkerColor(2);
  TGraph2D centerpoint_d3(nCent_det, &(det_centx[2][0]), &(det_centz[2][0]), &(det_centy[2][0]));
  centerpoint_d3.SetName("centerpoint_d3");
  centerpoint_d3.SetMarkerStyle(2);
  centerpoint_d3.SetMarkerColor(3);
  TGraph2D centerpoint_d4(nCent_det, &(det_centx[3][0]), &(det_centz[3][0]), &(det_centy[3][0]));
  centerpoint_d4.SetName("centerpoint_d4");
  centerpoint_d4.SetMarkerStyle(2);
  centerpoint_d4.SetMarkerColor(4);
  TGraph2D centerpoint_d5(nCent_det, &(det_centx[4][0]), &(det_centz[4][0]), &(det_centy[4][0]));
  centerpoint_d5.SetName("centerpoint_d5");
  centerpoint_d5.SetMarkerStyle(2);
  centerpoint_d5.SetMarkerColor(5);

  TCanvas flat_xy_canv("flat_xy_canv");
  wedge_flat_xy_gr.Draw("AP");
  ring_flat_xy_gr.Draw("same P");
  TCanvas flat_zy_canv("flat_zy_canv");
  wedge_flat_zy_gr.Draw("AP");
  ring_flat_zy_gr.Draw("same P");
  TCanvas flat_pt_canv("flat_pt_canv");
  wedge_flat_pt_gr.Draw("AP");
  ring_flat_pt_gr.Draw("same P");
  TCanvas flat_2d_canv("flat_2d_canv");
  wedge_flat_2d_gr.Draw("P");
  ring_flat_2d_gr.Draw("same P");
  TCanvas tilted_xy_canv("tilted_xy_canv");
  wedge_tilted_xy_gr.Draw("AP");
  ring_tilted_xy_gr.Draw("same P");
  TCanvas tilted_zy_canv("tilted_zy_canv");
  wedge_tilted_zy_gr.Draw("AP");
  ring_tilted_zy_gr.Draw("same P");
  TCanvas tilted_pt_canv("tilted_pt_canv");
  wedge_tilted_pt_gr.Draw("AP");
  ring_tilted_pt_gr.Draw("same P");
  TCanvas tilted_2d_canv("tilted_2d_canv");
  wedge_tilted_2d_gr.Draw("P");
  ring_tilted_2d_gr.Draw("same P");
  centerpoint.Draw("same P");
  TCanvas tilted_point("tilted_point");
  centerpoint.Draw("P");
  TCanvas tilted_point_det("tilted_point_det");
  wedge_tilted_2d_gr.Draw("P");
  ring_tilted_2d_gr.Draw("same P");
  centerpoint_d1.Draw("P same");
  centerpoint_d2.Draw("P same");
  centerpoint_d3.Draw("P same");
  centerpoint_d4.Draw("P same");
  centerpoint_d5.Draw("P same");

  TFile outfile("outfile.root","RECREATE");
  ring_flat_xy_gr.Write();
  wedge_flat_xy_gr.Write();
  flat_xy_canv.Write();
  ring_tilted_xy_gr.Write();
  wedge_tilted_xy_gr.Write();
  tilted_xy_canv.Write();
  ring_flat_zy_gr.Write();
  wedge_flat_zy_gr.Write();
  flat_zy_canv.Write();
  ring_tilted_zy_gr.Write();
  wedge_tilted_zy_gr.Write();
  tilted_zy_canv.Write();
  ring_flat_pt_gr.Write();
  wedge_flat_pt_gr.Write();
  flat_pt_canv.Write();
  ring_tilted_pt_gr.Write();
  wedge_tilted_pt_gr.Write();
  tilted_pt_canv.Write();
  ring_flat_2d_gr.Write();
  wedge_flat_2d_gr.Write();
  flat_2d_canv.Write();
  ring_tilted_2d_gr.Write();
  wedge_tilted_2d_gr.Write();
  tilted_2d_canv.Write();
  tilted_point.Write();
  tilted_point_det.Write();
  outfile.Close();

  delete [] sabdet;

  return 0;

}
