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

namespace Mask {

class Reaction {
public:
	Reaction();
	Reaction(int zt, int at, int zp, int ap, int ze, int ae);
	~Reaction();
	bool Calculate();
	void SetNuclei(int zt, int at, int zp, int ap, int ze, int ae);
	void SetNuclei(const Nucleus* nucs);
	void SetBeamKE(double bke);

	/*Setters and getters*/
	inline void SetLayeredTarget(LayeredTarget* targ) { target = targ; };
	inline void SetPolarRxnAngle(double theta) { m_theta = theta; };
	inline void SetAzimRxnAngle(double phi) { m_phi = phi; };
	inline void SetExcitation(double ex) { m_ex = ex; };
	inline void SetTarget(const Nucleus& nuc) { reactants[0] = nuc; };
	inline void SetTarget(int z, int a) { reactants[0] = Nucleus(z, a); };
	inline void SetProjectile(const Nucleus& nuc) { reactants[1] = nuc; };
	inline void SetProjectile(int z, int a) { reactants[1] = Nucleus(z, a); };
	inline void SetEjectile(const Nucleus& nuc) { reactants[2] = nuc; };
	inline void SetEjectile(int z, int a) { reactants[2] = Nucleus(z, a); };
	inline void SetResidual(const Nucleus& nuc) { reactants[3] = nuc; };
	inline void SetResidual(int z, int a) { reactants[3] = Nucleus(z, a); };
	inline void SetRxnLayer(int layer) { rxnLayer = layer; };
	inline void TurnOffResidualEloss() { resid_elossFlag = false; };
	inline void TurnOnResidualEloss() { resid_elossFlag = true; };
	inline bool IsDecay() { return decayFlag; };
	inline const Nucleus* GetNuclei() const { return &(reactants[0]); };
	inline const Nucleus& GetProjectile() const { return reactants[1]; };
	inline const Nucleus& GetTarget() const { return reactants[0]; };
	inline const Nucleus& GetEjectile() const { return reactants[2]; };
	inline const Nucleus& GetResidual() const { return reactants[3]; };
	inline int GetRxnLayer() { return rxnLayer; };

private:
	void CalculateReaction(); //target + project -> eject + resid
	void CalculateDecay(); //target -> light_decay (eject) + heavy_decay(resid)

	Nucleus reactants[4]; //0=target, 1=projectile, 2=ejectile, 3=residual
	LayeredTarget* target; //not owned by Reaction

	double m_bke, m_theta, m_phi, m_ex;

	int  rxnLayer;

	bool decayFlag, nuc_initFlag, resid_elossFlag;

};

};

#endif