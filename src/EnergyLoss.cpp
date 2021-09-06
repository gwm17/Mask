/*

EnergyLoss.cpp
Code for calculating the energy loss of a charged, massive particle through an arbitrary medium.
Based on code written by D.W. Visser at Yale for the original SPANC. Uses energy loss calulations
described by Ziegler in various SRIM textbooks.

Written by G.W. McCann Aug. 2020

*/

#include "Eloss_Tables.h"
#include "EnergyLoss.h"
#include "KinematicsExceptions.h"
#include <iostream>

namespace Mask {

	EnergyLoss::EnergyLoss() :
		ZP(-1), AP(-1), MP(-1.0), comp_denom(0)
	{
	}
	
	EnergyLoss::~EnergyLoss() {}
	
	/*Targets are defined by their atomic number, total number of nucleons, and their stoichiometry within the target compound*/
	void EnergyLoss::SetTargetComponents(std::vector<int>& Zt, std::vector<int>& At, std::vector<int>& Stoich) {
		comp_denom = 0;
		ZT = Zt;
		AT = At;
		for(unsigned int i=0; i<Stoich.size(); i++) {
			comp_denom += Stoich[i];
			if(ZT[i] > MAX_Z)
				throw ELossException("Maximum allowed target Z exceeded");
		}
		targ_composition.resize(Stoich.size());
	
		for(unsigned int i=0; i<Stoich.size(); i++)
			targ_composition[i] = Stoich[i]/comp_denom;
	}
	
	/*
	  Returns units of MeV; thickness in ug/cm^2; e_initial in units of MeV
	  Energy loss going through the target
	 */
	double EnergyLoss::GetEnergyLoss(int zp, int ap, double e_initial, double thickness) {
		if( ZP != zp) { 
			ZP = zp;
			AP = ap;
			MP = MassLookup::GetInstance()->FindMass(ZP, AP)*MEV2U;
		}
	
		double e_final = e_initial;
		double x_traversed = 0;
		double x_step = 0.25*thickness; //initial step in x
		double e_step = GetTotalStoppingPower(e_final)*x_step/1000.0; //initial step in e
		double e_threshold = 0.05*e_initial;
	
		int depth=0;
	
		
		if(thickness == 0.0 || e_initial == 0.0) 
			return 0;
	
		bool go = true;
		while(go) {
			//If intial guess of step size is too large, shrink until in range
			if(e_step/e_final > MAX_FRACTIONAL_STEP && depth < MAX_DEPTH) {
				depth++;
				x_step *= 0.5;
				e_step = GetTotalStoppingPower(e_final)*x_step/1000.0;
			} else if((x_step + x_traversed) >= thickness) { //last chunk
			  	go = false;
			  	x_step = thickness - x_traversed; //get valid portion of last chunk
			  	e_final -= GetTotalStoppingPower(e_final)*x_step/1000.0;
				if(e_final <= e_threshold)
					return e_initial;
			} else if(depth == MAX_DEPTH) {
				return e_initial;
			} else {
				x_traversed += x_step;
				e_step = GetTotalStoppingPower(e_final)*x_step/1000.0;
				e_final -= e_step;
				if(e_final <= e_threshold)
					return e_initial;
			}
		}
		return e_initial - e_final;
	}
	
	/*
	  Returns units of MeV; thickness in ug/cm^2; e_final in units of MeV
	  Energy loss going through the target using energy of particle after traversal
	*/
	double EnergyLoss::GetReverseEnergyLoss(int zp, int ap, double e_final, double thickness) {
		if( ZP != zp) { 
			ZP = zp;
			AP = ap;
			MP = MassLookup::GetInstance()->FindMass(ZP, AP)*MEV2U;
		}
	
		double e_initial = e_final;
		double x_traversed = 0.0;
		double x_step = 0.25*thickness; //initial step in x
		double e_step = GetTotalStoppingPower(e_initial)*x_step/1000.0; //initial step in E
	
		bool go = true;
		while(go) {
			if(e_step/e_initial > MAX_FRACTIONAL_STEP) {
				x_step *= 0.5;
				e_step = GetTotalStoppingPower(e_initial)*x_step/1000.0;
			} else if (x_traversed+x_step > thickness) {
				go = false;
				x_step = thickness - x_traversed;
				e_initial += GetTotalStoppingPower(e_initial)*x_step/1000.0;
			} else {
				x_traversed += x_step;
				e_step = GetTotalStoppingPower(e_initial)*x_step/1000.0;
				e_initial += e_step;
			}
		}
	
		return e_initial-e_final;
	}
	
	/*
	  Returns units of keV/(ug/cm^2)
	  Calculates Electronic stopping power by first calculating SE for hydrogen through the target and then using
	  corrections to calculate SE for projectile of interest
	*/
	double EnergyLoss::GetElectronicStoppingPower(double energy) {
		//Wants in units of keV
		energy *= 1000.0;
		double e_per_u = energy/MP;
		std::vector<double> values;
		if(e_per_u > MAX_H_E_PER_U) {
			throw ELossException("Exceeded maximum allowed energy per nucleon");
		} else if (e_per_u > 1000.0) {
			for(auto& z: ZT)
				values.push_back(Hydrogen_dEdx_High(e_per_u, energy, z));
		} else if (e_per_u > 10.0) {
			for(auto& z: ZT)
				values.push_back(Hydrogen_dEdx_Med(e_per_u, z));
		} else if (e_per_u > 0.0) {
			for(auto& z: ZT)
				values.push_back(Hydrogen_dEdx_Low(e_per_u, z));
		} else {
			throw ELossException("Negative energy per nucleon");
		}
	
		if(values.size() == 0)
			throw ELossException("Size of value array is 0. Unable to iterate over target components");
	
		if(ZP > 1) { //not hydrogen, need to account for effective charge
			for(unsigned int i=0; i<values.size(); i++)
				values[i] *= CalculateEffectiveChargeRatio(e_per_u, ZT[i]);
		}
		
		double stopping_total = 0;
		double conversion_factor = 0;
		for(unsigned int i=0; i<ZT.size(); i++) {
			conversion_factor += targ_composition[i]*NATURAL_MASS[ZT[i]];
			stopping_total += values[i]*targ_composition[i];
		}
		stopping_total *= AVOGADRO/conversion_factor;
	
		return stopping_total;
	}
	
	//Returns units of keV/(ug/cm^2)
	double EnergyLoss::GetNuclearStoppingPower(double energy) {
		energy *= 1000.0;
		double stopping_total = 0.0;
		double sn, x, epsilon, conversion_factor;
		for(unsigned int i=0; i<ZT.size(); i++) {
			x = (MP + NATURAL_MASS[ZT[i]]) * std::sqrt(std::pow(ZP, 2.0/3.0) + std::pow(ZT[i], 2.0/3.0));
			epsilon = 32.53*NATURAL_MASS[ZT[i]]*energy/(ZP*ZT[i]*x);
			sn = 8.462*(0.5*std::log(1.0+epsilon)/(epsilon+0.10718*std::pow(epsilon, 0.37544)))*ZP*ZT[i]*MP/x;
			conversion_factor = AVOGADRO/NATURAL_MASS[ZT[i]];
			stopping_total += sn*conversion_factor*targ_composition[i];
		}
	
		return stopping_total;
	}
	
	/*Wrapper function for aquiring total stopping (elec + nuc)*/
	double EnergyLoss::GetTotalStoppingPower(double energy) {
		if(ZP == 0)
			return GetNuclearStoppingPower(energy);
	
		return GetElectronicStoppingPower(energy)+GetNuclearStoppingPower(energy);
	}
	
	/*Charge rel to H*/
	double EnergyLoss::CalculateEffectiveChargeRatio(double e_per_u, int z) {
		double z_ratio;
		if(ZP == 2) {
			double ln_epu = std::log(e_per_u);
			double gamma = 1.0+(0.007+0.00005*z)*std::exp(-std::pow(7.6-ln_epu,2.0));
			double alpha = 0.7446 + 0.1429*ln_epu + 0.01562*std::pow(ln_epu, 2.0) - 0.00267*std::pow(ln_epu,3.0)
							+ 1.338E-6*std::pow(ln_epu,8.0);
			z_ratio = gamma*(1.0-std::exp(-alpha))*2.0;
		} else if (ZP == 3) {
			double ln_epu = std::log(e_per_u);
			double gamma = 1.0+(0.007+0.00005*z)*std::exp(-std::pow(7.6-ln_epu,2.0));
			double alpha = 0.7138+0.002797*e_per_u+1.348E-6*std::pow(e_per_u, 2.0);
			z_ratio = gamma*(1-std::exp(-alpha))*3.0;
		} else {
			double B = 0.886*std::pow(e_per_u/25.0, 0.5)/std::pow(ZP, 2.0/3.0);
			double A = B + 0.0378*std::sin(M_PI/2.0*B);
			z_ratio = 1.0 - std::exp(-A)*(1.034-0.1777*std::exp(-0.08114*ZP))*z;
		}
		return z_ratio*z_ratio; //for stopping power uses ratio sq. 
	}
	
	double EnergyLoss::Hydrogen_dEdx_Low(double e_per_u, int z) {
		return std::sqrt(e_per_u)*HYDROGEN_COEFF[z][0];
	}
	
	double EnergyLoss::Hydrogen_dEdx_Med(double e_per_u, int z) {
		double x = HYDROGEN_COEFF[z][1]*std::pow(e_per_u, 0.45);
		double y = HYDROGEN_COEFF[z][2]/e_per_u * std::log(1.0+HYDROGEN_COEFF[z][3]/e_per_u+HYDROGEN_COEFF[z][4]*e_per_u);
		return x*y/(x+y);
	}
	
	double EnergyLoss::Hydrogen_dEdx_High(double e_per_u, double energy, int z) {
		energy /= 1000.0; //back to MeV for ease of beta calc
		double beta_sq = energy * (energy+2.0*MP/MEV2U)/std::pow(energy+MP/MEV2U, 2.0);
		double alpha = HYDROGEN_COEFF[z][5]/beta_sq;
		double epsilon = HYDROGEN_COEFF[z][6]*beta_sq/(1.0-beta_sq) - beta_sq - HYDROGEN_COEFF[z][7];
		for(int i=1; i<5; i++)
			epsilon += HYDROGEN_COEFF[z][7+i]*std::pow(std::log(e_per_u), i);
	
		return alpha * std::log(epsilon);
	}
	
	//unimplemented
	double EnergyLoss::GetRange(double energy) {
		std::cerr<<"EnergyLoss::GetRange is not implemented! Returning 0.0"<<std::endl;
		return 0.0;
	}

}
