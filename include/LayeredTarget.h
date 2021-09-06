/*

LayeredTarget.h
Functional unit for targets in the Mask environment. Contains a
set (read: vector) of Targets for use in reaction calculations. In this
way handles situations such as carbon backed targets

Based on code by D.W. Visser written at Yale for the original SPANC

Written by G.W. McCann Aug. 2020

*/
#ifndef LAYEREDTARGET_H
#define LAYEREDTARGET_H

#include <vector>
#include <string>
#include "Target.h"

namespace Mask {

	class LayeredTarget {
	  
	public:
		LayeredTarget();
		~LayeredTarget();
		void AddLayer(std::vector<int>& Z, std::vector<int>& A, std::vector<int>& stoich, double thickness);
		double GetProjectileEnergyLoss(int zp, int ap, double startEnergy, int rxnLayer, double angle);
		double GetEjectileEnergyLoss(int ze, int ae, double startEnergy, int rxnLayer, double angle);
		double GetEjectileReverseEnergyLoss(int ze, int ae, double startEnergy, int rxnLayer, double angle);
		int FindLayerContaining(int Z, int A);
		inline int GetNumberOfLayers() { return layers.size(); }
		inline void SetName(std::string& n) { name = n; }
		inline const Target& GetLayerInfo(int index) { return layers[index]; }
		inline const std::string& GetName() { return name; }
	
	private:
		std::vector<Target> layers;
		std::string name;
	};

}

#endif
