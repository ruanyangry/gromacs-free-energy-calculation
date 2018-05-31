# _*_ coding: utf-8 _*_

'''
@author: Ruan Yang
Created on 2018.5.24
Free energy calculation.
Reference: http://www.bevanlab.biochem.vt.edu/Pages/Personal/justin/gmx-tutorials/free_energy
1. Steepest descents minimization
2. L-BFGS minimization
3. NVT equilibration
4. NPT equilibration
5. Data collection under an NPT ensemble
'''
import os

def Femdp(names,methods,nsteps,lambda_state,vdw_lambdas,coul_lambdas,bonded_lambdas,\
restraint_lambdas,mass_lambdas,temperature_lambdas):
	'''
	Simulation steps:
	1. Steepest descents minimization
	2. L-BFGS minimization
	3. NVT equilibration
	4. NPT equilibration
	5. Data collection under an NPT ensemble
	
	names(list): default: ["steep","lbfgs","nvt","npt","pd"]
	methods(list):steep,l-bfgs,sd,sd,sd, default ["steep","l-bfgs","sd","sd","sd"]
	nsteps(list): default [5000,5000,50000,50000,500000]
	lambda_state(list):
	vdw_lambdas(list):
	coul_lambdas(list):
	bonded_lambdas(list):
	restraint_lambdas(list):
	mass_lambdas(list):
	temperature_lambdas(list):
	'''
	for j in range(len(lambda_state)):
		os.mkdir("lambda-%d"%(j))
		os.chdir("lambda-%d"%(j))
		for i in range(5):
			with open("%s-%d.mdp"%(names[i],j),"w") as f:
				f.write("; Run control\n")
				f.write("integrator               = %s\n"%(methods[i]))
				if names[i] in ["steep","lbfgs"]:
					f.write("nsteps                   = %d\n"%(nsteps[i]))
					if names[i] == "lbfgs":
						f.write("define                   = -DFLEXIBLE\n")
					f.write("; EM criteria and other stuff\n")
					f.write("emtol                    = 100\n")
					f.write("emstep                   = 0.01\n")
					f.write("niter                    = 20\n")
					f.write("nbfgscorr                = 10\n")
					f.write("; Output control\n")
					f.write("nstlog                   = 1\n")
					f.write("nstenergy                = 1\n")
				else:
					f.write("tinit                    = 0\n")
					f.write("dt                       = 0.002\n")
					f.write("nsteps                   = %d\n"%(nsteps[i]))
					f.write("nstcomm                  = 100\n")
					f.write("; Output control\n")
					f.write("nstxout                  = 500\n")
					f.write("nstvout                  = 500\n")
					f.write("nstfout                  = 0\n")
					f.write("nstlog                   = 500\n")
					f.write("nstenergy                = 500\n")
					f.write("nstxout-compressed       = 0\n")
				f.write("; Neighborsearching and short-range nonbonded interactions\n")
				f.write("cutoff-scheme            = verlet\n")
				if names[i] in ["steep","lbfgs"]:
					f.write("nstlist                  = 1\n")
				else:
					f.write("nstlist                  = 20\n")
				f.write("ns_type                  = grid\n")
				f.write("pbc                      = xyz\n")
				f.write("rlist                    = 1.2\n")
				f.write("; Electrostatics\n")
				f.write("coulombtype              = PME\n")
				f.write("rcoulomb                 = 1.2\n")
				f.write("; van der Waals\n")
				f.write("vdwtype                  = cutoff\n")
				f.write("vdw-modifier             = potential-switch\n")
				f.write("rvdw-switch              = 1.0\n")
				f.write("rvdw                     = 1.2\n")
				f.write("; Apply long range dispersion corrections for Energy and Pressure\n")
				f.write("DispCorr                  = EnerPres\n")
				f.write("; Spacing for the PME/PPPM FFT grid\n")
				f.write("fourierspacing           = 0.12\n")
				f.write("; EWALD/PME/PPPM parameters\n")
				f.write("pme_order                = 6\n")
				f.write("ewald_rtol               = 1e-06\n")
				f.write("epsilon_surface          = 0\n")
				if names[i] in ["steep","lbfgs"]:
					f.write("; Temperature and pressure coupling are off during EM\n")
					f.write("tcoupl                   = no\n")
					f.write("pcoupl                   = no\n")
				else:
					f.write("; Temperature coupling\n")
					f.write("; tcoupl is implicitly handled by the sd integrator\n")
					f.write("tc_grps                  = system\n")
					f.write("tau_t                    = 1.0\n")
					f.write("ref_t                    = 300\n")
				if names[i] in ["nvt","npt","pd"]:
					if names[i] == "nvt":
						f.write("; Pressure coupling is off for NVT\n")
						f.write("Pcoupl                   = No\n")
						f.write("tau_p                    = 0.5\n")
						f.write("compressibility          = 4.5e-05\n")
						f.write("ref_p                    = 1.0\n")
					else:
						f.write("; Pressure coupling is on for NPT\n")
						f.write("Pcoupl                   = Parrinello-Rahman\n")
						f.write("tau_p                    = 1.0\n")
						f.write("compressibility          = 4.5e-05\n")
						f.write("ref_p                    = 1.0\n")
				f.write("; Free energy control stuff\n")
				f.write("free_energy              = yes\n")
				f.write("init_lambda_state        = %d\n"%(lambda_state[j]))
				f.write("delta_lambda             = 0\n")
				f.write("calc_lambda_neighbors    = 1\n")
				f.write("; Vectors of lambda specified here\n")
				f.write("; Each combination is an index that is retrieved from init_lambda_state for each simulation\n")
				#f.write("; init_lambda_state \n".join(str(lambda_state[j])))
				f.write("; init_lambda_state\n")
				f.write("vdw_lambdas              = "+"%s\n"%(" ".join([str(vdw_lambdas[i]) for i in range(len(vdw_lambdas))])))
				f.write("coul_lambdas             = "+"%s\n"%(" ".join([str(coul_lambdas[i]) for i in range(len(coul_lambdas))])))
				f.write("; We are not transforming any bonded or restrained interactions\n")
				f.write("bonded_lambdas           = "+"%s\n"%(" ".join([str(bonded_lambdas[i]) for i in range(len(bonded_lambdas))])))
				f.write("restraint_lambdas        = "+"%s\n"%(" ".join([str(restraint_lambdas[i]) for i in range(len(restraint_lambdas))])))
				f.write("; Masses are not changing (particle identities are the same at lambda = 0 and lambda = 1)\n")
				f.write("mass_lambdas             = "+"%s\n"%(" ".join([str(mass_lambdas[i]) for i in range(len(mass_lambdas))])))
				f.write("; Not doing simulated temperting here\n")
				f.write("temperature_lambdas      = "+"%s\n"%(" ".join([str(temperature_lambdas[i]) for i in range(len(temperature_lambdas))])))
				f.write("; Options for the decoupling\n")
				f.write("sc-alpha                 = 0.5\n")
				f.write("sc-coul                  = no\n")
				f.write("sc-power                 = 1.0\n")
				f.write("sc-sigma                 = 0.3\n")
				f.write("; name of moleculetype to decouple\n")
				f.write("couple-moltype           = TST\n")
				f.write("; only van der Waals interactions\n")
				f.write("couple-lambda0           = vdw\n")
				f.write("; turn off everything, in this case only vdW\n")
				f.write("couple-lambda1           = none\n")
				f.write("couple-intramol          = no\n")
				f.write("nstdhdl                  = 10\n")
				if names[i] in ["steep","lbfgs"]:
					f.write("; No velocities during EM\n")
					f.write("gen_vel                  = no\n")
				if names[i] in ["nvt","npt","pd"]:
					if names[i] == "nvt":
						f.write("; Generate velocities to start\n")
						f.write("gen_vel                  = yes\n")
						f.write("gen_temp                 = 300\n")
						f.write("gen_seed                 = -1\n")
					else:
						f.write("; Do not generate velocities\n")
						f.write("gen_vel                  = no\n")
				f.write("; options for bonds\n")
				if names[i] == "lbfgs":
					f.write("constraints              = none\n")
				else:
					f.write("constraints              = h-bonds\n")
				f.write("; Type of constraint algorithm\n")
				f.write("constraint-algorithm     = lincs\n")
				f.write("; Constrain the starting configuration\n")
				f.write("; since we are continuing from NVT\n")
				if names[i] == "nvt":
					f.write("continuation             = no\n")
				else:
					f.write("continuation             = yes\n")
				f.write("; Highest order in the expansion of the constraint coupling matrix\n")
				f.write("lincs-order              = 12\n")
		os.chdir("../")