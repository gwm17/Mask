# MASK: Monte cArlo Simulation of Kinematics
MASK is a Monte Carlo simulation of reaction kinematics intended for use with the Super-Enge Split-pole Spectrograph (SESPS) at Florida State University.
MASK is capable of simulating multi-step reaction-decay sequences, however in a purely kinematic sense, as it currently has no quantum mechanical input (this
is planned to be added in the next version). It is also capable of testing detector efficiency; this version contains the methods necessary to simulate the efficiency
of a reaction into the Silicion Array for Branching Ratio Detectors (SABRE).

## Building MASK
Download the repository from github using your favorite method. To build simply run

`make`

in the MASK directory, and the executable should be built and found in the bin directory.

## Running MASK
By default MASK is capable of simulating reactions of up to three steps. Here is a brief outline of each type:

0. A reaction of type 0 is a pure decay. It is assumed isotropic by default; any other case will require the modification of the code.
1. A reaction of type 1 is a pure reaction. It can incorporate all of the input file sampling parameters.
2. A reaction of type 2 is a reaction followed by a subsequent decay of the residual nucleus. Again, all sampling is allowed.
3. A reaction of type 3 is a reaction followed by a subsequent decay of the residual, followed by a decay of one of the products. Again, all sampling is allowed

For decays, a specific angular momentum L can be assumed. These are given in the input file as Decay1_AngularMomentum, and Decay2_AngularMomentum. This essentially modifies the center-of-mass angular distribution (as well as the lab frame). It is assumed that the decays in the center-of-mass frame are isotropic in phi (i.e. m=0). Decay1 corresponds to the first decay, if there are multiple steps, Decay2 to the second. If there are no decays, these parameters are not used (or if only one decay, Decay2_AngularMomentum is not used). The input file requires that the user include target information, which will be used to calculate energy loss for all of the reactants and reaction products. The target can contain layers, and each layer can be composed of a compound of elements with a given stoichiometry. If the user wishes to not include energy loss in the kinematics, simply give all target layers a thickness of 0. Note that more layers and more thickness = more time spent calculating energy loss. These energy loss methods are only applicable for solid targets, and should not be applied to gas or liquid targets. Energy loss calculations have a stated uncertainty of approximately five percent.

The default MASK program includes a calculation of SABRE efficiency, whose methods are contained in the SabreEfficiency and SabreDetector classes. This can be disabled by modifying the
main.cpp file appropriately.

In the input file the user also has the option to select to save the ROOT tree of the simulated data and the default plots. The options are yes or no. Yes saves them, no doesn't.

To run MASK simply do the following from the MASK directory:

`./bin/mask input.txt`

Input.txt can be replaced by any text file with the correct format.

## Requirements
MASK requires that ROOT is installed for data writting and visualization, as well as for random number generation. It also requires gsl to calculate Legendre Polynomials. Testing has been done only on ROOT 6. Mileage on all other ROOT versions will vary.