/*

LayeredTarget.h
Functional unit for targets in the SPANCRedux environment. Contains a
set (read: vector) of Targets for use in reaction calculations. In this
way handles situations such as carbon backed targets

Based on code by D.W. Visser written at Yale for the original SPANC

Written by G.W. McCann Aug. 2020

*/
#ifndef LAYEREDTARGET_H
#define LAYEREDTARGET_H

#include <vector>
#include <string>
#include <iostream>
#include "Target.h"

class LayeredTarget {
  
  public:
   LayeredTarget();
   ~LayeredTarget();
   void AddLayer(std::vector<int>& Z, std::vector<int>& A, std::vector<int>& stoich, double thickness);
   double GetProjectileEnergyLoss(int zp, int ap, double startEnergy, int rxnLayer, double angle);
   double GetEjectileEnergyLoss(int ze, int ae, double startEnergy, int rxnLayer, double angle);
   double GetEjectileReverseEnergyLoss(int ze, int ae, double startEnergy, int rxnLayer, double angle);
   int GetNumberOfLayers();
   int FindLayerContaining(int Z, int A);
   void SetName(std::string& n);
   Target& GetLayerInfo(int index);
   std::string& GetName();

  private:
   std::vector<Target> layers;
   std::string name;
};

#endif
