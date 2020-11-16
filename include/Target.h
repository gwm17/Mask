/*

Target.h
A basic target unit for use in the SPANCRedux environment. A target
is defined as a single compound with elements Z,A of a given stoichiometry
Holds an energy loss class

Based on code by D.W. Visser written at Yale for the original SPANC

Written by G.W. McCann Aug. 2020

*/
#ifndef TARGET_H
#define TARGET_H

#include <string>
#include <vector>
#include "EnergyLoss.h"

class Target {

  public:
    Target(double thick);
    ~Target();
    void SetElements(std::vector<int>& z, std::vector<int>& a, std::vector<int>& stoich);
    bool ContainsElement(int z, int a);
    double getEnergyLossTotal(int zp, int ap, double startEnergy, double angle);
    double getEnergyLossHalf(int zp, int ap, double startEnergy, double angle);
    double getReverseEnergyLossTotal(int zp, int ap, double finalEnergy, double angle);
    double getReverseEnergyLossHalf(int zp, int ap, double finalEnergy, double angle);
    double& GetThickness();
    int GetNumberOfElements();
    int GetElementZ(int index);
    int GetElementA(int index);
    int GetElementStoich(int index);
   

  private:
    EnergyLoss eloss;
    double thickness;
    std::vector<int> Z, A, Stoich;
    static constexpr double PI = 3.14159265358979323846;

};

#endif 
