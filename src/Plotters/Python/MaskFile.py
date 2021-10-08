#!/usr/bin/env python3

import numpy as np
import struct

class MaskFileData :
	def __init__(self, n):
		self.Z = np.zeros(n, dtype=int)
		self.A = np.zeros(n, dtype=int)
		self.dFlag = np.zeros(n, dtype=bool)
		self.E = np.zeros(n)
		self.KE = np.zeros(n)
		self.p = np.zeros(n)
		self.theta = np.zeros(n)
		self.phi = np.zeros(n)

class MaskFile:
	int_size = 4
	double_size = 8
	bool_size = 1

	def __init__(self, filename=""):
		self.eofFlag = False
		self.openFlag = False
		if filename != "" :
			self.Open(filename)

	def Open(self, filename):
		self.filename = filename
		self.file = open(self.filename, mode="rb")
		if self.file.closed  :
			self.openFlag = False
		else:
			self.openFlag = True

	def ReadHeader(self):
		data = self.file.read(2*self.int_size)

		(self.nsamples, self.rxntype) = struct.unpack("ii", data)
		self.datasize = (5*self.double_size+2*self.int_size+self.bool_size)
		self.datastr = "=ii?ddddd"

		if self.rxntype == 0:
			self.N_nuclei = 3
		elif self.rxntype == 1:
			self.N_nuclei = 4
		elif self.rxntype == 2:
			self.N_nuclei = 6
		elif self.rxntype == 3:
			self.N_nuclei = 8

	def ReadData(self):
		data = MaskFileData(self.N_nuclei)
		for i in range(self.N_nuclei):
			buffer = self.file.read(self.datasize)
			(data.Z[i], data.A[i], data.dFlag[i], data.E[i], data.KE[i], data.p[i], data.theta[i], data.phi[i]) = struct.unpack(self.datastr, buffer)
		if buffer == "":
			self.eofFlag = True

		return data

	def Close(self):
		self.file.close()

def main() :
	file = MaskFile(filename="/data1/gwm17/mask_tests/7Bedp_870keV_beam_50CD2.mask")
	file.ReadHeader()
	print("samples: ", file.nsamples, "rxntype:", file.rxntype, "datasize:", file.datasize)
	count=0
	for i in range(file.nsamples):
		file.ReadData()
		count += 1

	print("count:",count)
	print("eofFlag:",file.eofFlag)

	file.Close()

if __name__ == '__main__':
	main()