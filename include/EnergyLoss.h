/*

EnergyLoss.h
Code for calculating the energy loss of a charged, massive particle through an arbitrary medium.
Based on code written by D.W. Visser at Yale for the original SPANC. Uses energy loss calulations
described by Ziegler in various SRIM textbooks.

Written by G.W. McCann Aug. 2020

*/

#ifndef ENERGYLOSS_H
#define ENERGYLOSS_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include "MassLookup.h"

class EnergyLoss {
  
  public:
    EnergyLoss();
    ~EnergyLoss();
    double GetEnergyLoss(int zp, int ap, double e_initial, double thickness);
    double GetReverseEnergyLoss(int zp, int ap, double e_final, double thickness);
    double GetRange(double energy);
    void SetTargetComponents(std::vector<int>& Zt, std::vector<int>& At, std::vector<int>& Stoich);

  private:
    double GetElectronicStoppingPower(double energy);
    double GetNuclearStoppingPower(double energy);
    double GetTotalStoppingPower(double energy);
    double Hydrogen_dEdx_Low(double e_per_u, int z);
    double Hydrogen_dEdx_Med(double e_per_u, int z);
    double Hydrogen_dEdx_High(double e_per_u, double energy, int z);
    double CalculateEffectiveChargeRatio(double e_per_u, int z);

    int ZP, AP;
    double MP; //units of u, isotopic
    double comp_denom;
    std::vector<int> ZT, AT;
    std::vector<double> targ_composition;

    //constants for calculations
    static constexpr double MAX_FRACTIONAL_STEP = 0.001;
    static constexpr double MAX_H_E_PER_U = 100000.0;
    static constexpr double AVOGADRO = 0.60221367; //N_A times 10^(-24) for converting
    static constexpr double MEV2U = 1.0/931.4940954;
    static constexpr double PI = 3.14159265358979323846;
    static constexpr double H_RESTMASS = 938.27231; //MeV, for beta calc
};

#endif
