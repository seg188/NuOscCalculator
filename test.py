from Neutrino_Oscillation_Calculator import Neutrino_Oscillation_Calculator
import matplotlib.pyplot as plt
from scipy.stats import norm
from globaldef import *

import sys
sys.path.append('../NuOscProbExact/src')

import oscprob3nu
import hamiltonians3nu
from globaldefs import *

if __name__=='__main__':

	E = 1*GeV
	baselines = np.linspace(0, 5000*km, 700)

	osc_prob_ee_calc = []
	osc_prob_me_calc = []
	osc_prob_te_calc = []
	osc_prob_em_calc = []
	osc_prob_mm_calc = []
	osc_prob_tm_calc = []
	osc_prob_et_calc = []
	osc_prob_mt_calc = []
	osc_prob_tt_calc = []

	nu = NuNO_nufit_2()
	crust = EarthCrust()
	fm = FakeMaterial()
	#print(fm.coupling_strength())
	h_vacuum_energy_indep = \
		hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(  sin(nu.theta_12),
																sin(nu.theta_23),
																sin(nu.theta_13),
																nu.d_cp,
																nu.dM2_21,
																nu.dM2_31)

	# Units of VCC_EARTH_CRUST: [eV]
	h_matter = hamiltonians3nu.hamiltonian_3nu_matter(  h_vacuum_energy_indep,
													E,
													crust.coupling_strength()*nu.mode)
	for baseline in baselines:  # Baseline [km]
		Pee, Pem, Pet, Pme, Pmm, Pmt, Pte, Ptm, Ptt = \
			oscprob3nu.probabilities_3nu(h_matter, (baseline/km)*CONV_KM_TO_INV_EV)

		osc_prob_ee_calc.append(Pee)
		osc_prob_me_calc.append(Pme)
		osc_prob_te_calc.append(Pte)
		osc_prob_em_calc.append(Pem)
		osc_prob_mm_calc.append(Pmm)
		osc_prob_tm_calc.append(Ptm)
		osc_prob_et_calc.append(Pet)
		osc_prob_mt_calc.append(Pmt)
		osc_prob_tt_calc.append(Ptt)

		#print(Pme + Pte, Pee)
	
	osc_prob_ee = []
	osc_prob_me = []
	osc_prob_te = []
	osc_prob_em = []
	osc_prob_mm = []
	osc_prob_tm = []
	osc_prob_et = []
	osc_prob_mt = []
	osc_prob_tt = []


	for b in baselines:
		calc = Neutrino_Oscillation_Calculator(nu, crust)
		TEV = calc.Time_Evolution___Flavor(E, b)
		osc_prob_ee.append(np.absolute(TEV[0][0])**2)
		osc_prob_em.append(np.absolute(TEV[1][0])**2)
		osc_prob_et.append(np.absolute(TEV[2][0])**2)
		osc_prob_me.append(np.absolute(TEV[0][1])**2)
		osc_prob_mm.append(np.absolute(TEV[1][1])**2)
		osc_prob_mt.append(np.absolute(TEV[2][1])**2)
		osc_prob_te.append(np.absolute(TEV[0][2])**2)
		osc_prob_tm.append(np.absolute(TEV[1][2])**2)
		osc_prob_tt.append(np.absolute(TEV[2][2])**2)

		#print(np.absolute(TEV[1][0])**2 + np.absolute(TEV[2][0])**2, np.absolute(TEV[0][0])**2 )



	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.plot((baselines/km)/(E/GeV), osc_prob_ee_calc, label='analytic', linestyle='dashed', alpha=0.5)
	plt.plot((baselines/km)/(E/GeV), osc_prob_ee, label='calculator', alpha=0.5)
	ax.set_title('P_ee - Plane Wave in Earth Crust')
	ax.set_xlabel('L/E [km/GeV]')
	ax.set_ylabel('Probability')
	plt.legend()

	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.plot((baselines/km)/(E/GeV), osc_prob_em_calc, label='analytic', linestyle='dashed', alpha=0.5)
	plt.plot((baselines/km)/(E/GeV), osc_prob_em, label='calculator', alpha=0.5)
	
	ax.set_title('P_em - Plane Wave in Earth Crust')
	ax.set_xlabel('L/E [km/GeV]')
	ax.set_ylabel('Probability')
	plt.legend()

	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.plot((baselines/km)/(E/GeV), osc_prob_et_calc, label='analytic', linestyle='dashed', alpha=0.5)
	plt.plot((baselines/km)/(E/GeV), osc_prob_et, label='calculator', alpha=0.5)

	ax.set_title('P_et - Plane Wave in Earth Crust')
	ax.set_xlabel('L/E [km/GeV]')
	ax.set_ylabel('Probability')
	plt.legend()

	diff1 = [0]
	diff2 = [0]
	diff3 = [0]

	diff1 += [ (osc_prob_ee_calc[i] - osc_prob_ee[i])/osc_prob_ee_calc[i] for i in range(1, len(osc_prob_em_calc))]
	diff2 += [ (osc_prob_em_calc[i] - osc_prob_em[i])/osc_prob_em_calc[i] for i in range(1, len(osc_prob_em_calc))]
	diff3 += [ (osc_prob_et_calc[i] - osc_prob_et[i])/osc_prob_et_calc[i] for i in range(1, len(osc_prob_em_calc))]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.plot((baselines/km)/(E/GeV), diff1, label='ee')
	plt.plot((baselines/km)/(E/GeV), diff2, label='em')
	plt.plot((baselines/km)/(E/GeV), diff3, label='et')
	ax.set_title('Diff P_ee')
	ax.set_xlabel('L/E [km/GeV]')
	ax.set_ylabel('Probability')
	plt.legend()

	plt.show()