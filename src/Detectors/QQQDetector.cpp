#include "QQQDetector.h"

QQQDetector::QQQDetector(double R_in, double R_out, double deltaPhi, double phiCentral, double z, double x, double y) :
	m_Rinner(R_in), m_Router(R_out), m_deltaPhi(deltaPhi), m_phiCentral(phiCentral), m_translation(x,y,z)
{
	m_deltaR = (m_Router - m_Rinner)/nrings;
	m_deltaPhi_per_wedge = m_deltaPhi/nwedges;
	m_ZRot.SetAngle(m_phiCentral);
	m_ringCoords.resize(nrings);
	m_wedgeCoords.resize(nwedges);
	for(auto& ring : m_ringCoords) {
		ring.resize(4);
	}
	for(auto& wedge : m_wedgeCoords) {
		wedge.resize(4);
	}

	CalculateCorners();
}

QQQDetector::~QQQDetector() {}

void QQQDetector::CalculateCorners() {
	double x0, x1, x2, x3;
	double y0, y1, y2, y3;
	double z0, z1, z2, z3;

	//Generate flat ring corner coordinates
	for(int i=0; i<nrings; i++) {
		x0 = (m_Rinner + m_deltaR*(i+1))*std::cos(-m_deltaPhi/2.0);
		y0 = (m_Rinner + m_deltaR*(i+1))*std::sin(-m_deltaPhi/2.0);
		z0 = 0.0;
		m_ringCoords[i][0].SetVectorCartesian(x0, y0, z0);

		x1 = (m_Rinner + m_deltaR*(i))*std::cos(-m_deltaPhi/2.0);
		y1 = (m_Rinner + m_deltaR*(i))*std::sin(-m_deltaPhi/2.0);
		z1 = 0.0;
		m_ringCoords[i][1].SetVectorCartesian(x1, y1, z1);

		x2 = (m_Rinner + m_deltaR*(i))*std::cos(m_deltaPhi/2.0);
		y2 = (m_Rinner + m_deltaR*(i))*std::sin(m_deltaPhi/2.0);
		z2 = 0.0;
		m_ringCoords[i][2].SetVectorCartesian(x2, y2, z2);

		x3 = (m_Rinner + m_deltaR*(i+1))*std::cos(m_deltaPhi/2.0);
		y3 = (m_Rinner + m_deltaR*(i+1))*std::sin(m_deltaPhi/2.0);
		z3 = 0.0;
		m_ringCoords[i][3].SetVectorCartesian(x3, y3, z3);
	}

	//Generate flat wedge corner coordinates
	for(int i=0; i<nwedges; i++) {
		x0 = m_Router * std::cos(-m_deltaPhi/2.0 + i*m_deltaPhi_per_wedge);
		y0 = m_Router * std::sin(-m_deltaPhi/2.0 + i*m_deltaPhi_per_wedge);
		z0 = 0.0;
		m_wedgeCoords[i][0].SetVectorCartesian(x0, y0, z0);

		x1 = m_Rinner * std::cos(-m_deltaPhi/2.0 + i*m_deltaPhi_per_wedge);
		y1 = m_Rinner * std::sin(-m_deltaPhi/2.0 + i*m_deltaPhi_per_wedge);
		z1 = 0.0;
		m_wedgeCoords[i][1].SetVectorCartesian(x1, y1, z1);

		x2 = m_Rinner * std::cos(-m_deltaPhi/2.0 + (i+1)*m_deltaPhi_per_wedge);
		y2 = m_Rinner * std::sin(-m_deltaPhi/2.0 + (i+1)*m_deltaPhi_per_wedge);
		z2 = 0.0;
		m_wedgeCoords[i][2].SetVectorCartesian(x2, y2, z2);

		x3 = m_Router * std::cos(-m_deltaPhi/2.0 + (i+1)*m_deltaPhi_per_wedge);
		y3 = m_Router * std::sin(-m_deltaPhi/2.0 + (i+1)*m_deltaPhi_per_wedge);
		z3 = 0.0;
		m_wedgeCoords[i][3].SetVectorCartesian(x3, y3, z3);
	}

	for(int i=0; i<nrings; i++) {
		for(int j=0; j<4; j++) {
			m_ringCoords[i][j] = TransformCoordinates(m_ringCoords[i][j]);
		}
	}

	for(int i=0; i<nwedges; i++) {
		for(int j=0; j<4; j++) {
			m_wedgeCoords[i][j] = TransformCoordinates(m_wedgeCoords[i][j]);
		}
	}

}

Mask::Vec3 QQQDetector::GetTrajectoryCoordinates(double theta, double phi) {
	double z_to_detector = m_translation.GetZ();
	double rho_traj = z_to_detector*std::tan(theta);
	double r_traj = std::sqrt(rho_traj*rho_traj + z_to_detector*z_to_detector);
	double min_rho, max_rho, min_phi, max_phi;

	Mask::Vec3 result;

	for(auto& ring : m_ringCoords) {
		min_rho = ring[1].GetRho();
		max_rho = ring[0].GetRho();
		if(rho_traj < max_rho && rho_traj > min_rho) {
			for(auto& wedge : m_wedgeCoords) {
				min_phi = wedge[0].GetPhi();
				max_phi = wedge[3].GetPhi();
				if(phi < min_phi && phi < max_phi) {
					result.SetVectorSpherical(r_traj, theta, phi);
					break;
				}
			}
		}
	}
	
	return result;

}

std::pair<int,int> QQQDetector::GetTrajectoryRingWedge(double theta, double phi) {
	double z_to_detector = m_translation.GetZ();
	double rho_traj = z_to_detector*std::tan(theta);
	double min_rho, max_rho, min_phi, max_phi;


	for(int r=0; r<nrings; r++) {
		auto& ring = m_ringCoords[r];
		min_rho = ring[1].GetRho();
		max_rho = ring[0].GetRho();
		if(rho_traj < max_rho && rho_traj > min_rho) {
			for(int w=0; w<nwedges; w++) {
				auto& wedge = m_wedgeCoords[w];
				min_phi = wedge[0].GetPhi();
				max_phi = wedge[3].GetPhi();
				if(phi > min_phi && phi < max_phi) {
					return std::make_pair(r, w);
				}
			}
		}
	}
	
	return std::make_pair(-1, -1);

}

Mask::Vec3 QQQDetector::GetHitCoordinates(int ringch, int wedgech) {
	if(!CheckChannel(ringch) || !CheckChannel(wedgech)) {
		return Mask::Vec3();
	}

	double r_center  = m_Rinner + (0.5+ringch)*m_deltaR;
	double phi_center = -m_deltaPhi/2.0 + (0.5+wedgech)*m_deltaPhi_per_wedge;
	double x = r_center*std::cos(phi_center);
	double y = r_center*std::sin(phi_center);
	double z = 0;

	Mask::Vec3 hit(x, y, z);

	return TransformCoordinates(hit);
}