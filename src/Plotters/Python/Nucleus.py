#!/usr/bin/env python3

import numpy as np
from NucData import Masses

class Nucleus:
	deg2rad = np.pi/180.0
	def __init__(self, Z=0, A=0):
		self.Z = Z
		self.A = A
		if Z != 0 and A != 0:
			self.gsMass = Masses.GetMass(self.Z, self.A)
			self.vec4 = np.zeros(4)
			self.symbol = Masses.GetSymbol(self.Z, self.A)
			self.vec4[3] = self.gsMass

	def SetIsotope(self, Z, A):
		self.gsMass = Masses.GetMass(Z, A)
		self.symbol = Masses.GetSymbol(Z, A)
		self.Z = Z
		self.A = A
		self.vec4 = np.zeros(4)
		self.vec4[3] = self.gsMass

	def SetVectorCart(self, px, py, pz, E):
		self.vec4[0] = px
		self.vec4[1] = py
		self.vec4[2] = pz
		self.vec4[3] = E

	def SetVectorSpher(self, theta, phi, p, E):
		self.vec4[0] = p*np.sin(theta)*np.cos(phi)
		self.vec4[1] = p*np.sin(theta)*np.sin(phi)
		self.vec4[2] = p*np.cos(theta)
		self.vec4[3] = E

	def __add__(self, other):
		vec4 = self.vec4 + other.vec4
		newNuc = Nucleus(self.Z + other.Z, self.A + other.A)
		newNuc.SetVectorCart(vec4[0], vec4[1], vec4[2], vec4[3])
		return newNuc

	def __sub__(self, other):
		vec4 = self.vec4 - other.vec4
		newNuc = Nucleus(self.Z - other.Z, self.A - other.A)
		newNuc.SetVectorCart(vec4[0], vec4[1], vec4[2], vec4[3])
		return newNuc

	def __str__(self):
		return "Nucleus({0},{1}) with 4-vector({2})".format(self.Z, self.A, self.vec4)

	def GetP(self):
		return np.sqrt(self.vec4[0]**2.0 + self.vec4[1]**2.0 + self.vec4[2]**2.0)

	def GetInvMass(self):
		return np.sqrt(self.vec4[3]**2.0 - self.GetP()**2.0)

	def GetKE(self):
		return self.vec4[3] - self.GetInvMass()

	def GetTheta(self):
		return np.arccos(self.vec4[2]/self.GetP())

	def GetPhi(self):
		result = np.arctan2(self.vec4[1], self.vec4[0])
		if result < 0.0:
			result += 2.0*np.pi
		return result

	def GetExcitation(self):
		return self.GetInvMass() - self.gsMass

	def GetBoostToCMFrame(self):
		boost_vec = np.zeros(3)
		boost_vec[0] = self.vec4[0]/self.vec4[3]
		boost_vec[1] = self.vec4[1]/self.vec4[3]
		boost_vec[2] = self.vec4[2]/self.vec4[3]
		return boost_vec

	def ApplyBoost(self, boost_vec):
		beta2 = np.linalg.norm(boost_vec)**2.0
		gamma  = 1.0/np.sqrt(1.0 - beta2)
		bdotp = boost_vec[0]*self.vec4[0] + boost_vec[1]*self.vec4[1] + boost_vec[2]*self.vec4[2]
		gfactor = (gamma-1.0)/beta2 if beta2 > 0.0 else 0.0

		px = self.vec4[0]+gfactor*bdotp*boost_vec[0]+gamma*boost_vec[0]*self.vec4[3]
		py = self.vec4[1]+gfactor*bdotp*boost_vec[1]+gamma*boost_vec[1]*self.vec4[3]
		pz = self.vec4[2]+gfactor*bdotp*boost_vec[2]+gamma*boost_vec[2]*self.vec4[3]
		E = gamma*(self.vec4[3] + bdotp)
	
		self.SetVectorCart(px, py, pz, E);


def main():
	nuc = Nucleus(1,1)
	print("First",nuc)
	nuc2 = Nucleus(2,4)
	print("Second", nuc2)
	result = nuc + nuc2
	print("Addition", result)
	result2 = nuc2 - nuc
	print("Subtraction", result2)


if __name__ == '__main__':
	main()