# Neutrino_Oscillation_Calculator
# stephen_greenberg@berkeley.edu

from globaldef import *

import numpy as np
import scipy.linalg as la
from numpy import sin, arcsin
from numpy import cos, arccos
from numpy import exp
from numpy import sqrt
import matplotlib.pyplot as plt
from scipy.stats import norm


class Neutrino_Oscillation_Calculator:

	def __init__(self, nu, mat=None):
		self.nu = nu
	#-----------------------------------------------------------------------
	#construct as a product of 3 rotations (NOTE: excluding mayorana phases)
	#-----------------------------------------------------------------------
		self.U_23 = np.array([[1,			 0,				0],
						   [0, cos(nu.theta_23),	sin(nu.theta_23)],
	   					[0, -1*sin(nu.theta_23), cos(nu.theta_23)]
	
		])

		self.U_13 = np.array([[				 cos(nu.theta_13),		0,  sin(nu.theta_13)*exp(nu.d_cp * -1j)],
						   [							 0,		1,							 0],
	   					[-1*sin(nu.theta_13)*exp(nu.d_cp * 1j),		0,				 cos(nu.theta_13)]
		])

		self.U_12 = np.array([[   cos(nu.theta_12),	sin(nu.theta_12),	0],
						   [-1*sin(nu.theta_12),	cos(nu.theta_12),	0],
							[			   0,				0,	1]
	
		])

		self.U = self.U_23 @ self.U_13 @ self.U_12
		self.U_inv = la.inv(self.U)

		if not mat is None:
			A = mat.coupling_strength()*nu.mode
			self.V_f = np.array([[A, 0, 0], # electron-neutrino forward scattering potential in 
			 				 [0, 0, 0], # flavor basis
								 [0, 0, 0]
			])
		else:
			self.V_f = np.array([[0, 0, 0],
			 				 [0, 0, 0], 
								 [0, 0, 0]
			])

		self.Set_Hamiltonian___Mass__E_independent()

	def Set_Hamiltonian___Mass__E_independent(self): 
	#hamiltonian in mass basis!!!
	#hamiltonian as a function of the parameter E = (E_1 + E_2)/2
	#in relativistic limit where energies are parametrized by mass-squared differences
	#overall energy is neglected, only a determines overall phase
		E_21 = self.nu.dM2_21/(2)
		E_32 = self.nu.dM2_32/(2)
		E_1 = -1*E_21/2
		E_2 = E_21 + E_1
		E_3 = E_32 + E_2
	
		self.H = np.array([[E_1,   0,   0],
	   				[  0, E_2,   0],
		  			[  0,   0, E_3]
			   
		])

	def Hamiltonian___Mass(self, E):
		return 1/E * self.H

	def Hamiltonian___Flavor(self, E):
		return self.U @ self.Hamiltonian___Mass(E) @ self.U_inv

	def Full_H___Flavor(self, E):
		return self.Hamiltonian___Flavor(E) + self.V_f

	def Time_Evolution___Flavor(self, E, L): ## L is baseline in m
		return la.expm( (-1j)*L/(c*hbar)*self.Full_H___Flavor(E))

	def Hamiltonian___Mass_Wave_Packet(self, E, E_width): 
	#hamiltonian in mass basis!!!
	#hamiltonian as a function of the parameter E = (E_1 + E_2)/2
	#in relativistic limit where energies are parametrized by mass-squared differences
	#overall energy is neglected, only a determines overall phase
		H = np.array([[0,0,0],
	   			[0,0,0],
		  		 [0,0,0]
	   
   		 ])
		nsample = 100
		_sum = 0
		energy_samples = np.linspace(E-5*E_width, E+5*E_width, nsample)
		sample_width = 10.*float(E_width)/float(nsample)
		for esample in energy_samples:
	   		weight = norm.pdf((esample-E)/E_width)*sample_width
	   		_sum += weight
	   		E_21 = self.nu.dM2_21/(2*esample)
	   		E_32 = self.nu.dM2_32/(2*esample)
	   		E_1 = (esample-E) - 1*E_21/2
	   		E_2 = E_21 + E_1
	   		E_3 = E_32 + E_2
	   		H_inc = np.array([[E_1,   0,   0],
		  					  [  0, E_2,   0],
			 		    	[  0,   0, E_3]
			   
	   			])*weight
	   		H = H + H_inc
		return H/_sum

	def Hamiltonian___Flavor_Wave_Packet(self, E, E_width):
		return self.U @ self.Hamiltonian___Mass_Wave_Packet(E, E_width) @ self.U_inv

	def Full_H___Flavor_Wave_Packet(self, E, E_width):
		return self.Hamiltonian___Flavor_Wave_Packet(E, E_width) + self.V_f

	def Time_Evolution___Flavor_Wave_Packet(self, E, E_width, L): ## L is baseline in m
		return la.expm( (-1j)*L/(c*hbar)*self.Full_H___Flavor_Wave_Packet(E, E_width))


class Neutrino_Oscillation_Calculator_Sterile:

	def __init__(self, nu, mat=None):
		self.nu = nu
	#-----------------------------------------------------------------------
	#construct as a product of 3 rotations (NOTE: excluding mayorana phases)
	#-----------------------------------------------------------------------
		self.U_23 = np.array([[1,			 0,				0, 0],
						   [0, cos(nu.theta_23),	sin(nu.theta_23), 0],
	   					[0, -1*sin(nu.theta_23), cos(nu.theta_23), 0],
	   					[			   0,				0,	0,  1]
	
		])

		self.U_13 = np.array([[				 cos(nu.theta_13),		0,  sin(nu.theta_13)*exp(nu.d_cp * -1j), 0],
						   [							 0,		1,							 0, 0],
	   					[-1*sin(nu.theta_13)*exp(nu.d_cp * 1j),		0,				 cos(nu.theta_13), 0],
	   					[			   0,				0,	0,  1]
		])

		self.U_12 = np.array([[   cos(nu.theta_12),	sin(nu.theta_12),	0, 0],
						   [-1*sin(nu.theta_12),	cos(nu.theta_12),	0, 0],
							[			   0,				0,	1,  0],
							[			   0,				0,	0,  1]
	
		])

		self.C = np.array([[1,  0, 0, 0],
						   [0, 1,	0, 0],
							[ 0,0,	cos(nu.c),  sin(nu.c)],
							[0, 0,	-sin(nu.c),  cos(nu.c)]
	
		])

		self.B = np.array([[1,  0, 0, 0],
						   [0, cos(nu.b),	0, sin(nu.b)],
							[ 0,0,	1,  0],
							[0,-sin(nu.b),	0,  cos(nu.b)]
	
		])

		self.A = np.array([[cos(nu.a), 0, 0, sin(nu.a)],
						   [0, 1,	0, 0],
							[ 0,0,	1,  0],
							[-sin(nu.a),0,	0,  cos(nu.a)]
	
		])



		self.U = self.U_23 @ self.U_13 @ self.U_12 @ self.A @ self.B @ self.C
		self.U_inv = la.inv(self.U)

		if not mat is None:
			A = mat.coupling_strength()*nu.mode
			self.V_f = np.array([[A, 0, 0, 0], # electron-neutrino forward scattering potential in 
			 				     [0, 0, 0, 0], # flavor basis
			 				     [0, 0, 0, 0],
								 [0, 0, 0, 0]
			])
		else:
			self.V_f = np.array([[0, 0, 0, 0],
			 				     [0, 0, 0, 0],
			 				     [0, 0, 0, 0],  
								 [0, 0, 0, 0]
			])

		self.Set_Hamiltonian___Mass__E_independent()

	def Set_Hamiltonian___Mass__E_independent(self): 
	#hamiltonian in mass basis!!!
	#hamiltonian as a function of the parameter E = (E_1 + E_2)/2
	#in relativistic limit where energies are parametrized by mass-squared differences
	#overall energy is neglected, only a determines overall phase
		E_21 = self.nu.dM2_21/(2)
		E_32 = self.nu.dM2_32/(2)
		E_1 = -1*E_21/2
		E_2 = E_21 + E_1
		E_3 = E_32 + E_2
		E_4 = self.nu.OM4
	
		self.H = np.array([[E_1, 0,   0, 0],
	   					   [0, E_2,   0, 0],
		  			 	   [0,   0, E_3, 0],
		  			 	   [0, 0, 0, E_4]
			   
		])

	def Hamiltonian___Mass(self, E):
		return 1/E * self.H

	def Hamiltonian___Flavor(self, E):
		return self.U @ self.Hamiltonian___Mass(E) @ self.U_inv

	def Full_H___Flavor(self, E):
		return self.Hamiltonian___Flavor(E) + self.V_f

	def Time_Evolution___Flavor(self, E, L): ## L is baseline in m
		return la.expm( (-1j)*L/(c*hbar)*self.Full_H___Flavor(E))











