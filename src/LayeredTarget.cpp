/*

LayeredTarget.h
Functional unit for targets in the SPANCRedux environment. Contains a
set (read: vector) of Targets for use in reaction calculations. In this
way handles situations such as carbon backed targets

Based on code by D.W. Visser written at Yale for the original SPANC

Written by G.W. McCann Aug. 2020

*/
#include "LayeredTarget.h"

LayeredTarget::LayeredTarget() {}

LayeredTarget::~LayeredTarget() {}

/*Add in a Target made of a compound defined by a set of Zs, As, Ss, and a thickness*/
void LayeredTarget::AddLayer(std::vector<int>& Z, std::vector<int>& A, std::vector<int>& stoich, double thickness) {
  Target t(thickness);
  t.SetElements(Z, A, stoich);
  layers.push_back(t);
}

/*
  Here projectile refers to the incoming reactant particle (i.e. the beam)
  Calculates energy loss assuming that the reaction occurs in the middle of the target
  Note that the layer order can matter!
*/
double LayeredTarget::GetProjectileEnergyLoss(int zp, int ap, double startEnergy, int rxnLayer, double angle) {

  if(rxnLayer < 0 || rxnLayer > layers.size()) {
    std::cerr<<"Reaction layer in eloss calculation is not in range! Returning 0"<<std::endl;
    return 0.0;
  }

  double eloss = 0.0;
  double newEnergy = startEnergy;
  for(int i=0; i<=rxnLayer; i++) {
    if(i == rxnLayer) {
      eloss += layers[i].getEnergyLossHalf(zp, ap, newEnergy, angle);
      newEnergy = startEnergy - eloss;
    } else {
      eloss += layers[i].getEnergyLossTotal(zp, ap, newEnergy, angle);
      newEnergy = startEnergy-eloss;
    }
  }

  return eloss;
}

/*
  Here ejectile refers to the outgoing reactant particle
  Calculates energy loss assuming that the reaction occurs in the middle of the target
  Note that the layer order can matter!
*/
double LayeredTarget::GetEjectileEnergyLoss(int ze, int ae, double startEnergy, int rxnLayer, double angle) {

  if(rxnLayer < 0 || rxnLayer > layers.size()) {
    std::cerr<<"Reaction layer in eloss calculation is not in range! Returning 0"<<std::endl;
    return 0.0;
  }

  double eloss = 0.0;
  double newEnergy = startEnergy;
  for(unsigned int i=rxnLayer; i<layers.size(); i++) {
    if(i == rxnLayer) {
      eloss += layers[i].getEnergyLossHalf(ze, ae, newEnergy, angle);
      newEnergy = startEnergy - eloss;
    } else {
      eloss += layers[i].getEnergyLossTotal(ze, ae, newEnergy, angle);
      newEnergy = startEnergy - eloss;
    }
  }

  return eloss;
}

/*ReverseEnergyLoss version of GetEjectileEnergyLoss*/
double LayeredTarget::GetEjectileReverseEnergyLoss(int ze, int ae, double startEnergy, int rxnLayer, double angle) {

  if(rxnLayer < 0 || rxnLayer > layers.size()) {
    std::cerr<<"Reaction layer in eloss calculation is not in range! Returning 0"<<std::endl;
    return 0.0;
  }

  double eloss = 0.0;
  double newEnergy = startEnergy;
  for(int i=(layers.size()-1); i>=rxnLayer; i--) {
    if(i == rxnLayer) {
      eloss += layers[i].getReverseEnergyLossHalf(ze, ae, newEnergy, angle);
      newEnergy = startEnergy + eloss;
    } else {
      eloss += layers[i].getReverseEnergyLossTotal(ze, ae, newEnergy, angle);
      newEnergy = startEnergy + eloss;
    }
  }

  return eloss;
}

/*Getters and Setters*/

int LayeredTarget::GetNumberOfLayers() {
  return layers.size();
}

int LayeredTarget::FindLayerContaining(int Z, int A) {
  for(unsigned int i=0; i<layers.size(); i++) {
    if(layers[i].ContainsElement(Z, A)) return i;
  }
  return -1;
}

void LayeredTarget::SetName(std::string& n) {
  name = n;
}

Target& LayeredTarget::GetLayerInfo(int index) {
  return layers[index];
}

std::string& LayeredTarget::GetName() {
  return name;
}

