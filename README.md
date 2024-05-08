# Mask: Monte cArlo Simulation of Kinematics

Mask is a Monte Carlo simulation of reaction kinematics for use detector systems at Florida State University.
Mask is capable of simulating multi-step kinematic reaction-decay sequences, storing data in ROOT Trees, after which the kinematic data can be fed to a detector geometry for efficiency testing. Currently geometries for ANASEN and SABRE are included in the code.

## Building Mask

Dowload the repository from GitHub using `git clone --recursive https://github.com/gwm17/Mask.git`. CMake is use to build the project; in most environments you can build Mask using the following methods:

- `mkdir build`
- `cd build`
- `cmake ..`
- `make`

Executables will be installed to the repostiory's `bin` directory. Libraries will be installed to the `lib` directory.

By default Mask builds for release. To build for debug replace `cmake ..` with `cmake -DCMAKE_BUILD_TYPE=Debug ..`. Mask uses CMake to find the installed ROOT libraries and headers.

## Using the kinematics simulation

By default Mask is capable of simulating reactions of up to three steps. In the configuration file, the reaction is specified by the `ReactionChain`, which is a list of reaction specifications. Each reaction specification has a `Type` which is either Reaction or Decay, and a list of `Reactants` and sampling parameters. To run Mask simply do the following from the Mask repository:

`./bin/Kinematics <your_config>.yaml`

`<your_config.yaml>` is a YAML configuration file. An example is given in the repository named `kinematics.yaml` and can be replaced by any yaml file with the correct format.

### Reaction

To specify a reaction you need 3 reactants, each specified by a `Z` (proton number) and `A` (mass number). The first reactant is the target, the second the projectile (beam), and the third is the ejectile. The residual is calculated for you assuming conservation of proton and mass number (no weak decays). Reactions require a `ThetaType` which is either Lab or CenterOfMass. This tells the simulation to sample the reaction angle in the Lab (useful for detector constraints) or CenterOfMass (useful for the correct distribution in open systems). `ThetaMin`, `ThetaMax`, `PhiMin`, `PhiMax` all specify the limits in angles for the simulation (here Theta is the polar (reaction) angle, and Phi is the azimuthal angle). `BeamEnergyMean` and `BeamEnergySigma` specify the beam energy distribution to be used, assuming the energy is a gaussian distribution. `ResidualExcitationMean` and `ResidualExcitationSigma` specify the excited state energy of the residual nucleus assuming a gaussian distribution.

### Decay

To specify a decay you need 2 reactants, each specified by a `Z` (proton number) and `A` (mass number). The first is the parent nucleus, and the second is one of the decay products. The other product will be calculated for you assuming conservation of proton and mass number (no weak decays). Decays only have a phi limit for sampling, and are only sampled in the CenterOfMass frame. The calculated decay product (residual) can still have an excitation distribution specified. Additionally, Decays can have an angular distribution file specified. The file contains the weights for a Legendre Polynomial series description of an angular distribution. An example file of an isotropic distribution is included with the repository in the `etc` directory.

### Limitations

Mask can only accept certain types of chains. That is, a Decay only chain is allowed, and the a Reaction + up to 2 subsequent Decays are allowed. Any other types of chains are not supported at this time. Mask will check your chain to make sure it complies with these requirements.

## Using the detector geometry simulation

Detector geometry is encoded using ROOT math libraries in the `src/Detectors` folder. Two different detector geometries are already present: SPS-SABRE and ANASEN. To add a new geometry, follow the guidelines outlined by each of these cases.

To choose which detector scheme is run, modify the main function in `src/Detectors/main.cpp`. The included geometries also have options to do an internal geometry consistency check and print out coordinates for drawing the detector arrays, which can be useful for testing.

To run the geometry code, one needs to provide an input file containing the following: the path of a Mask kinematics data file, the path to which data should be written, the path to a file containing a list of dead channels (optional, if not used, write None for the path), the number of threads to be used by the thread pool, and a keyword for the array type (current options are Sabre or Anasen)

To run Detectors use the format

`./bin/Detectors <your_config>.yaml`

`<your_config>.yaml` is a YAML configuration file. An example, `detector.yaml` is included in the repository.

## Data visualization

All data is saved as ROOT trees of std::vectors of Mask::Nucleus classes. To enable this, a ROOT dictionary is generated and linked into a shared library found in the `lib` directory of the repository. This allows the user to link to the shared library for accessing and analyzing the data generated by Mask.

Mask also provides a default visualization tool called RootPlot. RootPlot is run as

`./bin/RootPlot <datafile> <outputfile>`

where the datafile can be either the datafile from Mask or the datafile from DetEff. The outputfile is saved in the ROOT file format.

## Requirements

ROOT version 6.22 or greater is required
CMake version 3.0 or greater is required
