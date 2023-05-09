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
	energies = np.linspace(5e8, 8e9*eV, 1000)

	osc_prob_ee_calc = []
	osc_prob_me_calc = []
	osc_prob_te_calc = []
	osc_prob_em_calc = []
	osc_prob_mm_calc = []
	osc_prob_tm_calc = []
	osc_prob_et_calc = []
	osc_prob_mt_calc = []
	osc_prob_tt_calc = []

	nu = NuNO_nufit_4()
	nu.d_cp = -1*Pi/2
	nu.mode=1
	plot_title = 'Nu Oscillations Probs. in Earth\'s Crust - ' + nu.name
	if nu.mode < 0: plot_title = 'anti-'+plot_title

	savename = 'nu_'+nu.name 
	if nu.mode < 0: savename = 'anti'+savename

	crust = EarthCrust()
	baseline = 1284.9*km


	h_vacuum_energy_indep = \
		hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(  sin(nu.theta_12),
																sin(nu.theta_23),
																sin(nu.theta_13),
																nu.d_cp,
																nu.dM2_21,
																nu.dM2_31)

	
	for E in energies:  # Baseline [km]
		# Units of VCC_EARTH_CRUST: [eV]
		h_matter = hamiltonians3nu.hamiltonian_3nu_matter(  h_vacuum_energy_indep,
														E,
														crust.coupling_strength()*nu.mode)
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


	### my calculator
	
	osc_prob_ee = []
	osc_prob_me = []
	osc_prob_te = []
	osc_prob_em = []
	osc_prob_mm = []
	osc_prob_tm = []
	osc_prob_et = []
	osc_prob_mt = []
	osc_prob_tt = []


	osc_prob_me_posdcp = []
	osc_prob_mm_posdcp = []
	osc_prob_mt_posdcp = []

	osc_prob_me_negdcp = []
	osc_prob_mm_negdcp = []
	osc_prob_mt_negdcp = []

	for E in energies:
		nu.d_cp = 0
		calc = Neutrino_Oscillation_Calculator(nu, crust)
		TEV = calc.Time_Evolution___Flavor(E, baseline)
		osc_prob_ee.append(np.absolute(TEV[0][0])**2)
		osc_prob_em.append(np.absolute(TEV[1][0])**2)
		osc_prob_et.append(np.absolute(TEV[2][0])**2)
		osc_prob_me.append(np.absolute(TEV[0][1])**2)
		osc_prob_mm.append(np.absolute(TEV[1][1])**2)
		osc_prob_mt.append(np.absolute(TEV[2][1])**2)
		osc_prob_te.append(np.absolute(TEV[0][2])**2)
		osc_prob_tm.append(np.absolute(TEV[1][2])**2)
		osc_prob_tt.append(np.absolute(TEV[2][2])**2)

		nu.d_cp= -1*Pi/2
		calc = Neutrino_Oscillation_Calculator(nu, crust)
		TEV = calc.Time_Evolution___Flavor(E, baseline)
		osc_prob_me_negdcp.append(np.absolute(TEV[0][1])**2)
		osc_prob_mm_negdcp.append(np.absolute(TEV[1][1])**2)
		osc_prob_mt_negdcp.append(np.absolute(TEV[2][1])**2)

		nu.d_cp= Pi/2
		calc = Neutrino_Oscillation_Calculator(nu, crust)
		TEV = calc.Time_Evolution___Flavor(E, baseline)
		osc_prob_me_posdcp.append(np.absolute(TEV[0][1])**2)
		osc_prob_mm_posdcp.append(np.absolute(TEV[1][1])**2)
		osc_prob_mt_posdcp.append(np.absolute(TEV[2][1])**2)

		#print(np.absolute(TEV[1][0])**2 + np.absolute(TEV[2][0])**2, np.absolute(TEV[0][0])**2 )

	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.plot(energies/GeV, osc_prob_me_calc, label='analytic dcp=-pi/2', linestyle='dashed', c='red')
	if nu.mode > 0:
		plt.fill_between(energies/GeV,osc_prob_me_negdcp, label='dcp=-pi/2', alpha=0.75, color='green')
		plt.fill_between(energies/GeV,osc_prob_me,  label='dcp=0', alpha=0.75, color='blue')
		plt.fill_between(energies/GeV,osc_prob_me_posdcp, label='dcp=pi/2', alpha=0.75, color='orange')
	else:
		plt.fill_between(energies/GeV,osc_prob_me_posdcp, label='dcp=pi/2', alpha=0.75, color='orange')
		plt.fill_between(energies/GeV,osc_prob_me,  label='dcp=0', alpha=0.75, color='blue')
		plt.fill_between(energies/GeV,osc_prob_me_negdcp, label='dcp=-pi/2', alpha=0.75, color='green')

	ax.set_title(plot_title)
	ax.set_xlabel('E [GeV]')
	ax.set_ylabel('Prob.')
	ax.set_xscale('log')
	plt.legend()
	plt.savefig('plots/'+savename+'.png')

	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.plot(energies/GeV, osc_prob_mm_calc, label='analytic', linestyle='dashed', alpha=0.5)
	plt.plot(energies/GeV, osc_prob_mm, label='calculator', alpha=0.5)
	
	ax.set_title('P_mm - Plane Wave in Earth Crust')
	ax.set_xlabel('E [GeV]')
	ax.set_ylabel('Probability')
	ax.set_xscale('log')
	plt.legend()

	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.plot(energies/GeV, osc_prob_mt_calc, label='analytic', linestyle='dashed', alpha=0.5)
	plt.plot(energies/GeV, osc_prob_mt, label='calculator', alpha=0.5)

	ax.set_title('P_mt - Plane Wave in Earth Crust')
	ax.set_xlabel('E [GeV]')
	ax.set_ylabel('Probability')
	ax.set_xscale('log')
	plt.legend()

	diff1 = [0]
	diff2 = [0]
	diff3 = [0]

	diff1 += [ (osc_prob_me_calc[i] - osc_prob_me[i])/osc_prob_me_calc[i] for i in range(1, len(osc_prob_em_calc))]
	diff2 += [ (osc_prob_mm_calc[i] - osc_prob_mm[i])/osc_prob_mm_calc[i] for i in range(1, len(osc_prob_em_calc))]
	diff3 += [ (osc_prob_mt_calc[i] - osc_prob_mt[i])/osc_prob_mt_calc[i] for i in range(1, len(osc_prob_em_calc))]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.plot(energies/GeV, diff1, label='me')
	plt.plot(energies/GeV, diff2, label='mm')
	plt.plot(energies/GeV, diff3, label='mt')
	ax.set_title('Analytic to Calc. Difference')
	ax.set_xlabel('E [GeV]')
	ax.set_xscale('log')
	ax.set_ylabel('Probability')
	plt.legend()

	plt.show()