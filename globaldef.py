import numpy as np
from numpy import sin, arcsin
from numpy import cos, arccos
from numpy import exp
from numpy import sqrt

# numerical and mathematical constants go here

#UNITS
eV = 1.000000000
KeV = 10**3*eV
MeV = 10**6*eV
GeV = 10**9*eV

m = 1.000000000
cm = 10**(-2)*m
km = 10**3*m
s = 1.000000000
g = 5.61*(10**32)*eV


#PHYSICAL CONSTANTS
c = 299792458.0*(m/s)
hbar = 6.582119569*10**(-16)*eV*s
m_p = 938.272088*MeV
G_F = (1.1663787*(10**(-5)))/GeV**2 #electroweak coupling constant

#MATH CONSTANTS
Pi = np.pi


class NuNO_nufit_4:
	def __init__(self):
		self.theta_12 = 0.5903       
		self.theta_13 = 0.150
		self.theta_23 = 0.866        #sin2(θ23) = 0.546 ± 0.021 (normal ordering!!! different for inverted order)
		self._d_cp = 0.5*Pi                          	#1.36 (+ 0.20 − 0.16) π rad
		self.dM2_21 = 7.39*10**(-5)*(eV**2)          #(7.53 ± 0.18) × 10−5 eV2
		self.dM2_32 = 2.451*10**(-3)*(eV**2)         #(2.453 ± 0.033) × 10−3 eV2 (normal ord_ering!!!)
		self.dM2_31 = (self.dM2_32 + self.dM2_21)
		self.mode = 1                                # +1 neutrino, -1 anti-neutrino
		self.name = 'NO'

	@property
	def d_cp(self):
		return self._d_cp*self.mode

	@d_cp.setter
	def d_cp(self, val):
		self._d_cp = val
		
class NuIO_nufit_4:
	def __init__(self):
		self.theta_12 = 0.5903       
		self.theta_13 = 0.151   
		self.theta_23 = 0.869         #sin2(θ23) = 0.546 ± 0.021 (normal ordering!!! different for inverted order)
		self._d_cp = 0.5*Pi                          	#1.36 (+ 0.20 − 0.16) π rad
		self.dM2_21 = 7.39*10**(-5)*(eV**2)          #(7.53 ± 0.18) × 10−5 eV2
		self.dM2_32 = -2.512*10**(-3)*(eV**2)         #(2.453 ± 0.033) × 10−3 eV2 (normal ordering!!!)
		self.dM2_31 = (self.dM2_32 + self.dM2_21)
		self.mode = 1  
		self.name = 'IO'

	@property
	def d_cp(self):
		return self._d_cp*self.mode

	@d_cp.setter
	def d_cp(self, val):
		self._d_cp = val                              # +1 neutrino, -1 anti-neutrino


class NuNO_nufit_2:
	def __init__(self):
		self.theta_12 = arcsin(sqrt(0.308))         
		self.theta_13 = arcsin(sqrt(2.163e-2))   
		self.theta_23 = arcsin(sqrt(0.440))          #sin2(θ23) = 0.546 ± 0.021 (normal ordering!!! different for inverted order)
		self._d_cp = 289*Pi/180                          	#1.36 (+ 0.20 − 0.16) π rad
		self.dM2_21 = 7.49e-5*(eV**2)         
		self.dM2_31 = 2.526*10**(-3)*(eV**2)         #(2.453 ± 0.033) × 10−3 eV2 (normal ordering!!!)
		self.dM2_32 = (self.dM2_31 - self.dM2_21)
		self.mode = 1    

	@property
	def d_cp(self):
		return self._d_cp*self.mode

	@d_cp.setter
	def d_cp(self, val):
		self._d_cp = val                              # +1 neutrino, -1 anti-neutrino
  

class NuIO_nufit_2:
	def __init__(self):
		self.theta_12 = arcsin(sqrt(0.308))         
		self.theta_13 = arcsin(sqrt(2.175e-2))   
		self.theta_23 = arcsin(sqrt(0.584))          #sin2(θ23) = 0.546 ± 0.021 (normal ordering!!! different for inverted order)
		self._d_cp = 269*Pi/180                          	#1.36 (+ 0.20 − 0.16) π rad
		self.dM2_21 = 7.49e-5*(eV**2)         
		self.dM2_32 = -2.518*10**(-3)*(eV**2)         #(2.453 ± 0.033) × 10−3 eV2 (normal ordering!!!)
		self.dM2_31 = (self.dM2_32 + self.dM2_21)
		self.mode = 1

	@property
	def d_cp(self):
		return self._d_cp*self.mode

	@d_cp.setter
	def d_cp(self, val):
		self._d_cp = val                              # +1 neutrino, -1 anti-neutrino

class NuSterile:
	def __init__(self):
		self.theta_12 = 0.5903       
		self.theta_13 = 0.150
		self.theta_23 = 0.866        #sin2(θ23) = 0.546 ± 0.021 (normal ordering!!! different for inverted order)
		self._d_cp = 0.5*Pi                          	#1.36 (+ 0.20 − 0.16) π rad
		self.dM2_21 = 7.39*10**(-5)*(eV**2)          #(7.53 ± 0.18) × 10−5 eV2
		self.dM2_32 = 2.451*10**(-3)*(eV**2)         #(2.453 ± 0.033) × 10−3 eV2 (normal ord_ering!!!)
		self.dM2_31 = (self.dM2_32 + self.dM2_21)
		self.mode = 1                                # +1 neutrino, -1 anti-neutrino
		self.name = 'NO'
		self.OM4 = 1.
		self.a = 0.1
		self.b = 0.1
		self.c = 0.1


	@property
	def d_cp(self):
		return self._d_cp*self.mode

	@d_cp.setter
	def d_cp(self, val):
		self._d_cp = val
        


A_nuc = 12. #mass number for carbon
Z = 6. #atomic number for carbon
rho_earth_crust = 2.848*(g/(cm**3))

class Material:
	def __init__(self, A_nuc, Z, rho):
		self.A_nuc = A_nuc
		self.Z = Z
		self.rho = rho

	def electron_density(self):
		return self.rho*self.Z/(self.A_nuc*m_p)

	def coupling_strength(self):
		return sqrt(2)*G_F*self.electron_density()
Carbon = Material(A_nuc, Z, rho_earth_crust)

class EarthCrust:
	def __init__(self):
		self.N_e = 2.848 * 2.9805e23 
		print(self.N_e)
	#	self.N_e =  8488460000.0 # https://arxiv.org/pdf/1707.02322.pdf
		'''
		self.element_Z_weight_ppm = { \
		'O' : (8, 462700.),
		'Si' : (14, 276300.),
		'Al' : (13, 81667.),
		'Fe' : (26, 56600.),
		'Ca' : (20, 42667.),
		'Na' : (11, 24700.),
		'K' : (19, 20567.),
		'Mg' : (12, 29000.),
		'Ti' : (22, 6133.)
		}
		'''

	def coupling_strength(self):
		return G_F*self.N_e*sqrt(2)

class FakeMaterial:
	def __init__(self):
		self.N_e = 10**31 #number density of electrons in earth crust
	def coupling_strength(self):
		#print(self.N_e*sqrt(2)*G_F)
		return 0
		#return self.N_e*sqrt(2)*G_F
		