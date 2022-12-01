/*
	Reaction.h
	Reaction is a class which implements either a decay or scattering reaction. As such it requires either
	3 (decay) or 4 (scattering) nuclei to perform any calcualtions. I also links together the target, which provides
	energy loss calculations, with the kinematics. Note that Reaction does not own the LayeredTarget.

	--GWM Jan. 2021
*/
#ifndef REACTION_H
#define REACTION_H

#include "Nucleus.h"
#include "LayeredTarget.h"
#include "RxnType.h"

namespace Mask {

	class Reaction
	{
	public:
		Reaction();
		Reaction(Nucleus* target, Nucleus* projectile, Nucleus* ejectile, Nucleus* residual);
		~Reaction();
		bool Calculate(); //do sim

		void BindNuclei(Nucleus* target, Nucleus* projectile, Nucleus* ejectile, Nucleus* residual);
		void SetBeamKE(double bke);
		void SetEjectileThetaType(RxnThetaType type);
	
		void SetLayeredTarget(LayeredTarget* targ) { m_layeredTarget = targ; };

		void SetPolarRxnAngle(double theta) { m_theta = theta; };
		void SetAzimRxnAngle(double phi) { m_phi = phi; };
		void SetExcitation(double ex) { m_ex = ex; };
		void SetReactionDepth(double depth) { m_rxnDepth = depth; }

		void BindTarget(Nucleus* nuc) { m_target = nuc; };
		void BindProjectile(Nucleus* nuc) { m_projectile = nuc; };
		void BindEjectile(Nucleus* nuc) { m_ejectile = nuc; };
		void BindResidual(Nucleus* nuc) { m_residual = nuc; };

		void SetRxnLayer(std::size_t layer) { m_rxnLayer = layer; };
		void SetResidualEnergyLoss(bool isEloss) { m_isResidEloss = isEloss; };

		bool IsDecay() const { return m_isDecay; };

		std::size_t GetRxnLayer() const { return m_rxnLayer; };
	
	private:
		void CalculateDecay(); //target -> light_decay (eject) + heavy_decay(resid)
		void CalculateReaction(); //target + project -> eject + resid
		void CalculateReactionThetaLab();
		void CalculateReactionThetaCM();
	
		//Reactants -> NOT OWNED BY RXN
		Nucleus* m_target;
		Nucleus* m_projectile;
		Nucleus* m_ejectile;
		Nucleus* m_residual;

		LayeredTarget* m_layeredTarget; //not owned by Reaction
	
		double m_bke, m_theta, m_phi, m_ex, m_rxnDepth;
	
		int m_rxnLayer;
		RxnThetaType m_ejectThetaType; 
	
		bool m_isDecay, m_isInit, m_isResidEloss;
	};

}

#endif