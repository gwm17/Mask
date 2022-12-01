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
#include <random>
#include "RandomGenerator.h"
#include "Target.h"

namespace Mask {

	class LayeredTarget
	{
	public:
		LayeredTarget();
		~LayeredTarget();
		void AddLayer(const std::vector<int>& Z, const std::vector<int>& A, const std::vector<int>& stoich, double thickness);
		double GetProjectileEnergyLoss(int zp, int ap, double startEnergy, std::size_t rxnLayer, double angle, double rxnDepth);
		double GetEjectileEnergyLoss(int ze, int ae, double startEnergy, std::size_t rxnLayer, double angle, double rxnDepth);
		double GetEjectileReverseEnergyLoss(int ze, int ae, double startEnergy, std::size_t rxnLayer, double angle, double rxnDepth);
		std::size_t FindLayerContaining(int Z, int A);
		std::size_t GetNumberOfLayers() { return m_layers.size(); }
		void SetName(std::string& n) { m_name = n; }
		const Target& GetLayerInfo(int index) { return m_layers[index]; }
		const std::string& GetName() { return m_name; }
	
	private:
		std::vector<Target> m_layers;
		std::string m_name;
	};

}

#endif
