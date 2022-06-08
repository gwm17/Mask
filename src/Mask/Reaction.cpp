/*
	Reaction.cpp
	Reaction is a class which implements either a decay or scattering reaction. As such it requires either
	3 (decay) or 4 (scattering) nuclei to perform any calcualtions. I also links together the target, which provides
	energy loss calculations, with the kinematics. Note that Reaction does not own the LayeredTarget.

	--GWM Jan. 2021
*/
#include "Reaction.h"
#include "KinematicsExceptions.h"

namespace Mask {

	Reaction::Reaction() :
		target(nullptr), m_bke(0), m_theta(0), m_phi(0), m_ex(0), rxnLayer(0), m_eject_theta_type(lab), nuc_initFlag(false), resid_elossFlag(false)
	{
	}
	
	Reaction::Reaction(int zt, int at, int zp, int ap, int ze, int ae) :
		target(nullptr), m_bke(0), m_theta(0), m_phi(0), m_ex(0), rxnLayer(0), m_eject_theta_type(lab), resid_elossFlag(false)
	{
		SetNuclei(zt, at, zp, ap, ze, ae);
	}
	
	Reaction::~Reaction()
	{
	}
	
	bool Reaction::Calculate() {
	
		if(!nuc_initFlag) 
			return false;
	
		if(decayFlag) {
			CalculateDecay();
			return true;
		} else {
			CalculateReaction();
			return true;
		}
	}
	
	//Deep copy of nucleus array
	void Reaction::SetNuclei(const Nucleus* nucs) {
		reactants[0] = nucs[0];
		reactants[1] = nucs[1];
		reactants[2] = nucs[2];
		reactants[3] = nucs[3];
		nuc_initFlag = true;
	}
	
	void Reaction::SetNuclei(int zt, int at, int zp, int ap, int ze, int ae) {
		int zr, ar;
		reactants[0] = Nucleus(zt, at);
		reactants[2] = Nucleus(ze, ae);
		if(ap == 0) {
			decayFlag = true;
			zr = zt - ze;
			ar = at - ae;
		} else {
			reactants[1] = Nucleus(zp, ap);
			decayFlag = false;
			zr = zt + zp - ze;
			ar = at + ap - ae;
		}
	
		if(zr < 0 || ar <= 0) {
			nuc_initFlag = false;
		} else {
			reactants[3] = Nucleus(zr, ar);
			nuc_initFlag = true;
		}
	}
	
	void Reaction::SetBeamKE(double bke) {
		if(!nuc_initFlag || decayFlag) 
			return;
	
		m_bke = bke - target->GetProjectileEnergyLoss(reactants[1].GetZ(), reactants[1].GetA(), bke, rxnLayer, 0);
	}
	
	void Reaction::SetEjectileThetaType(int type) {
		if(decayFlag) return;
		if(type != center_of_mass && type != lab) return;
	
		m_eject_theta_type = type;
	}
	
	//Methods given by Iliadis in Nuclear Physics of Stars, Appendix C
	//For use with lab frame restricted angles. May not give appropriate disribution for ejectile
	void Reaction::CalculateReactionThetaLab() {
		reactants[0].SetVectorCartesian(0.,0.,0.,reactants[0].GetGroundStateMass());
		double beam_pz = std::sqrt(m_bke*(m_bke + 2.0 * reactants[1].GetGroundStateMass()));
		double beam_E = m_bke + reactants[1].GetGroundStateMass();
		reactants[1].SetVectorCartesian(0.,0.,beam_pz,beam_E);
	
		double Q = reactants[0].GetGroundStateMass() + reactants[1].GetGroundStateMass() - (reactants[2].GetGroundStateMass() + reactants[3].GetGroundStateMass() + m_ex);
	
		double Ethresh = -Q*(reactants[2].GetGroundStateMass()+reactants[3].GetGroundStateMass())/(reactants[2].GetGroundStateMass() + reactants[3].GetGroundStateMass() - reactants[1].GetGroundStateMass());
		if(m_bke < Ethresh) {
			throw EnergyThresholdException();
		}
	
		double term1 = sqrt(reactants[1].GetGroundStateMass()*reactants[2].GetGroundStateMass()*m_bke)/(reactants[2].GetGroundStateMass()+reactants[3].GetGroundStateMass())*cos(m_theta);
		double term2 = (m_bke*(reactants[3].GetGroundStateMass() - reactants[1].GetGroundStateMass()) + reactants[3].GetGroundStateMass()*Q)/(reactants[3].GetGroundStateMass() + reactants[2].GetGroundStateMass());
		double sqrt_pos_ejectKE = term1 + sqrt(term1*term1 + term2);
		double sqrt_neg_ejectKE = term1 - sqrt(term1*term1 + term2);
		double ejectKE;
		if(sqrt_pos_ejectKE > 0) {
			ejectKE = sqrt_pos_ejectKE*sqrt_pos_ejectKE;
		} else {
			ejectKE = sqrt_neg_ejectKE*sqrt_neg_ejectKE;
		}
		double ejectP = std::sqrt(ejectKE*(ejectKE + 2.0 * reactants[2].GetGroundStateMass()));
		double ejectE = ejectKE + reactants[2].GetGroundStateMass();
	
		reactants[2].SetVectorSpherical(m_theta, m_phi, ejectP, ejectE);
	
		reactants[3] = reactants[0] + reactants[1] - reactants[2];
	
		ejectKE -= target->GetEjectileEnergyLoss(reactants[2].GetZ(), reactants[2].GetA(), ejectKE, rxnLayer, m_theta);
		ejectP = std::sqrt(ejectKE*(ejectKE + 2.0 * reactants[2].GetGroundStateMass()));
		ejectE = ejectKE + reactants[2].GetGroundStateMass();
		reactants[2].SetVectorSpherical(m_theta, m_phi, ejectP, ejectE);
	
		if(resid_elossFlag) {
			double residKE = reactants[3].GetKE() - target->GetEjectileEnergyLoss(reactants[3].GetZ(), reactants[3].GetA(), reactants[3].GetKE(), rxnLayer, reactants[3].GetTheta());
			double residP = std::sqrt(residKE*(residKE + 2.0*reactants[3].GetInvMass()));
			double residE = residKE + reactants[3].GetInvMass();
			reactants[3].SetVectorSpherical(reactants[3].GetTheta(), reactants[3].GetPhi(), residP, residE);
		}
	}
	
	//Methods from original ANASEN. Gives proper distribution for inverse kinematics.
	void Reaction::CalculateReactionThetaCM() {
		//Target assumed at rest, with 0 excitation energy
		reactants[0].SetVectorCartesian(0.,0.,0.,reactants[0].GetGroundStateMass());
		double beam_pz = std::sqrt(m_bke*(m_bke + 2.0 * reactants[1].GetGroundStateMass()));
		double beam_E = m_bke + reactants[1].GetGroundStateMass();
		reactants[1].SetVectorCartesian(0.,0.,beam_pz,beam_E);
	
	
		double Q = reactants[0].GetGroundStateMass() + reactants[1].GetGroundStateMass() - (reactants[2].GetGroundStateMass() + reactants[3].GetGroundStateMass() + m_ex);
	
		double Ethresh = -Q*(reactants[2].GetGroundStateMass()+reactants[3].GetGroundStateMass())/(reactants[2].GetGroundStateMass() + reactants[3].GetGroundStateMass() - reactants[1].GetGroundStateMass());
		if(m_bke < Ethresh) {
			throw EnergyThresholdException();
		}
		
		auto parent = reactants[0] + reactants[1];
		double boost2lab[3];
		double boost2cm[3];
		const double* boost = parent.GetBoost();
		for(int i=0; i<3; i++) {
			boost2lab[i] = boost[i];
			boost2cm[i] = boost[i]*(-1.0);
		}
		parent.ApplyBoost(boost2cm);
		double ejectE_cm = (std::pow(reactants[2].GetGroundStateMass(), 2.0) - std::pow(reactants[3].GetGroundStateMass() + m_ex, 2.0) + std::pow(parent.GetE(),2.0))/(2.0*parent.GetE());
		double ejectP_cm = std::sqrt(ejectE_cm*ejectE_cm - std::pow(reactants[2].GetGroundStateMass(), 2.0));
		reactants[2].SetVectorSpherical(m_theta, m_phi, ejectP_cm, ejectE_cm);
		reactants[2].ApplyBoost(boost2lab);
		reactants[3] = reactants[0] + reactants[1] - reactants[2];
	
		double ejectKE = reactants[2].GetKE();
		double ejectP = reactants[2].GetP();
		double ejectE = reactants[2].GetE();
		//energy loss for ejectile (after reaction!)
		ejectKE -= target->GetEjectileEnergyLoss(reactants[2].GetZ(), reactants[2].GetA(), ejectKE, rxnLayer, reactants[2].GetTheta());
		ejectP = std::sqrt(ejectKE*(ejectKE + 2.0 * reactants[2].GetGroundStateMass()));
		ejectE = ejectKE + reactants[2].GetGroundStateMass();
		reactants[2].SetVectorSpherical(reactants[2].GetTheta(), reactants[2].GetPhi(), ejectP, ejectE);
	
		//if on, get eloss for residual (after reaction!)
		if(resid_elossFlag) {
			double residKE = reactants[3].GetKE() - target->GetEjectileEnergyLoss(reactants[3].GetZ(), reactants[3].GetA(), reactants[3].GetKE(), rxnLayer, reactants[3].GetTheta());
			double residP = std::sqrt(residKE*(residKE + 2.0*reactants[3].GetInvMass()));
			double residE = residKE + reactants[3].GetInvMass();
			reactants[3].SetVectorSpherical(reactants[3].GetTheta(), reactants[3].GetPhi(), residP, residE);
		}
	}
	
	void Reaction::CalculateReaction() {
		switch(m_eject_theta_type) {
			case center_of_mass:
			{
				CalculateReactionThetaCM();
				break;
			}
			case lab:
			{
				CalculateReactionThetaLab();
				break;
			}
		}
	}
	
	//Calculate in CM, where decay is isotropic
	void Reaction::CalculateDecay() {
	
		double Q = reactants[0].GetInvMass() - reactants[2].GetGroundStateMass() - reactants[3].GetGroundStateMass();
		if(Q < 0) {
			throw QValueException();
		}
	
		const double* boost = reactants[0].GetBoost();
		double boost2cm[3];
		double boost2lab[3];
		for(int i=0; i<3; i++) {
			boost2lab[i] = boost[i];
			boost2cm[i] = boost[i]*(-1.0);
		}
	
		reactants[0].ApplyBoost(&(boost2cm[0]));
		double ejectE_cm = (reactants[2].GetGroundStateMass()*reactants[2].GetGroundStateMass() - reactants[3].GetGroundStateMass()*reactants[3].GetGroundStateMass() + reactants[0].GetE()*reactants[0].GetE())/
					(2.0*reactants[0].GetE());
		double ejectP_cm = std::sqrt(ejectE_cm*ejectE_cm - reactants[2].GetGroundStateMass()*reactants[2].GetGroundStateMass());
	
		reactants[2].SetVectorSpherical(m_theta, m_phi, ejectP_cm, ejectE_cm);
		reactants[2].SetThetaCM(m_theta);
	
		reactants[0].ApplyBoost(boost2lab);
		reactants[2].ApplyBoost(boost2lab);
	
		reactants[3] = reactants[0] - reactants[2];
	
		//energy loss for the *light* break up nucleus
		double ejectKE = reactants[2].GetKE() - target->GetEjectileEnergyLoss(reactants[2].GetZ(), reactants[2].GetA(), reactants[2].GetKE(), rxnLayer, reactants[2].GetTheta());
		double ejectP = std::sqrt(ejectKE*(ejectKE + 2.0*reactants[2].GetGroundStateMass()));
		double ejectE = ejectKE + reactants[2].GetGroundStateMass();
		reactants[2].SetVectorSpherical(reactants[2].GetTheta(), reactants[2].GetPhi(), ejectP, ejectE);
	
		//if on, get eloss for *heavy* break up nucleus
		if(resid_elossFlag) {
	
			double residKE = reactants[3].GetKE() - target->GetEjectileEnergyLoss(reactants[3].GetZ(), reactants[3].GetA(), reactants[3].GetKE(), rxnLayer, reactants[3].GetTheta());
			double residP = std::sqrt(residKE*(residKE + 2.0*reactants[3].GetInvMass()));
			double residE = residKE + reactants[3].GetInvMass();
			reactants[3].SetVectorSpherical(reactants[3].GetTheta(), reactants[3].GetPhi(), residP, residE);
		}
	}

}



