#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from Nucleus import Nucleus
from MaskFile import MaskFile, MaskFileData
from NucData import Masses
import pickle
import sys

def PlotData(inputname, outputname):
	rad2deg = 180.0/np.pi

	datafile = MaskFile(inputname)
	datafile.ReadHeader()

	print("MaskFile opened -- rxntype:", datafile.rxntype, "number of samples:", datafile.nsamples)

	data = MaskFileData(datafile.N_nuclei)

	ke = np.zeros((datafile.N_nuclei, datafile.nsamples))
	ke_d = np.zeros((datafile.N_nuclei, datafile.nsamples))
	theta = np.zeros((datafile.N_nuclei, datafile.nsamples))
	theta_d = np.zeros((datafile.N_nuclei, datafile.nsamples))
	phi = np.zeros((datafile.N_nuclei, datafile.nsamples))
	phi_d = np.zeros((datafile.N_nuclei, datafile.nsamples))
	detect_mask = np.ones((datafile.N_nuclei, datafile.nsamples), dtype=bool)
	names = []
	for i in range(datafile.N_nuclei):
		names.append(" ")

	nuc = Nucleus()

	for i in range(datafile.nsamples):
		data = datafile.ReadData()

		if i == 0:
			for j in range(datafile.N_nuclei):
				names[j] = Masses.GetSymbol(data.Z[j], data.A[j])

		for j in range(datafile.N_nuclei):
			nuc.SetIsotope(data.Z[j], data.A[j])
			nuc.SetVectorSpher(data.theta[j], data.phi[j], data.p[j], data.E[j])

			ke[j][i] = nuc.GetKE()
			theta[j][i] = data.theta[j]*rad2deg
			phi[j][i] = data.phi[j]*rad2deg
			if data.dFlag[j] == True:
				ke_d[j][i] = data.KE[j]
				theta_d[j][i] = data.theta[j]*rad2deg
				phi_d[j][i] = data.phi[j]*rad2deg
			else:
				detect_mask[j][i] = False

	datafile.Close()

	#Remove empty values from detection arrays
	final_theta_d = theta_d[detect_mask]
	final_phi_d = phi_d[detect_mask]
	final_ke_d = ke_d[detect_mask]

	#figs = {}
	#axes = {}

	'''
	for i in range(len(names)):
		figs[i], axes[i] = plt.subplots(2,2)
	'''
	fig, axes = plt.subplots(len(names)-1,4)
	fig.set_size_inches(12, 12)

	for i in range(1, len(names)):
		'''
		axes[i][0][0].plot(theta[i], ke[i], marker=',', linestyle='None')
		axes[i][0][0].set_title(names[i]+" KE vs. Theta")
		axes[i][0][0].set_xlabel(r"$\theta$ (degrees)")
		axes[i][0][0].set_ylabel("KE (MeV)")

		axes[i][0][1].plot(phi[i], ke[i], marker=",", linestyle='None')
		axes[i][0][1].set_title(names[i]+" KE vs. Phi")
		axes[i][0][1].set_xlabel(r"$\phi$ (degrees)")
		axes[i][0][1].set_ylabel("KE (MeV)")

		axes[i][1][0].plot(theta_d[i], ke_d[i], marker=',', linestyle='None')
		axes[i][1][0].set_title(names[i]+" KE vs. Theta -- Detected")
		axes[i][1][0].set_xlabel(r"$\theta$ (degrees)")
		axes[i][1][0].set_ylabel("KE (MeV)")

		axes[i][1][1].plot(phi_d[i], ke_d[i], marker=",", linestyle='None')
		axes[i][1][1].set_title(names[i]+" KE vs. Phi -- Detected")
		axes[i][1][1].set_xlabel(r"$\phi$ (degrees)")
		axes[i][1][1].set_ylabel("KE (MeV)")
		'''

		axes[i-1][0].plot(theta[i], ke[i], marker=',', linestyle='None')
		axes[i-1][0].set_title(names[i]+" KE vs. Theta")
		axes[i-1][0].set_xlabel(r"$\theta$ (degrees)")
		axes[i-1][0].set_ylabel("KE (MeV)")

		axes[i-1][1].plot(phi[i], ke[i], marker=",", linestyle='None')
		axes[i-1][1].set_title(names[i]+" KE vs. Phi")
		axes[i-1][1].set_xlabel(r"$\phi$ (degrees)")
		axes[i-1][1].set_ylabel("KE (MeV)")

		axes[i-1][2].plot(theta_d[i], ke_d[i], marker=',', linestyle='None')
		axes[i-1][2].set_title(names[i]+" KE vs. Theta -- Detected")
		axes[i-1][2].set_xlabel(r"$\theta$ (degrees)")
		axes[i-1][2].set_ylabel("KE (MeV)")

		axes[i-1][3].plot(phi_d[i], ke_d[i], marker=",", linestyle='None')
		axes[i-1][3].set_title(names[i]+" KE vs. Phi -- Detected")
		axes[i-1][3].set_xlabel(r"$\phi$ (degrees)")
		axes[i-1][3].set_ylabel("KE (MeV)")

	plt.tight_layout()
	plt.show(block=True)

	print("Writing figure to file:", outputname)
	with open(outputname, "wb") as outfile:
		pickle.dump(fig, outfile)
		outfile.close()

	print("Finished.")

def main():
	if len(sys.argv) == 3:
		PlotData(sys.argv[1], sys.argv[2])
	else:
		print("Unable to run PyPlotter, incorrect number of arguments! Need an input datafile name, and an output plot file name")

if __name__ == '__main__':
	main()



