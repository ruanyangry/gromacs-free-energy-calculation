# _*_ coding:utf-8 _*_

'''
@author: Ruan Yang
Created on 2018.5.28
Free energy calculation.
Reference: http://www.bevanlab.biochem.vt.edu/Pages/Personal/justin/gmx-tutorials/free_energy
'''

import GenerateFF as GF

def fsolsh(nstates,name,nt):
	'''
	Generated free energy workflow shell file.
	
	nstates: how many states defined in .mdp file.
	name(list):organic names
	nt: nodes used to parallel computing 
	'''
	for i in range(len(name)):
		#os.chdir("/molStructure/%s"%(name[i]))
		#os.mkdir("FE")
		with open("FE-%s.sh"%(name[i]),"w") as f:
			f.write("#!/bin/bash\n")
			f.write('#$ -S /bin/sh\n')
			f.write('#$ -N test\n')
			f.write('#$ -j y\n')
			f.write('#$ -o .\n')
			f.write('#$ -e .\n')
			f.write('#$ -cwd\n')
			f.write('#$ -q all.q\n')
			f.write('#$ -masterq all.q@node9\n')
			f.write('#$ -pe thread %d-%d\n'%(nt,nt))
			f.write('source ~/.bashrc\n')
			f.write('hash -r\n')
			f.write('export PATH=/export/home/ry/gromacs-5.0/bin:PATH\n')
			f.write('\n')
			f.write("\n")
			f.write("# Author: Ruan Yang\n")
			f.write("# Email: ruanyang_njut@163.com\n")
			f.write("\n")
			for j in range(nstates):
				f.write("echo 'Staring state %d'\n"%(j))
				if j==0:
					f.write("cd ../molStructure/%s\n"%(name[i]))
					f.write("mkdir FE\n")
					f.write("cd FE\n")
				else:
					f.write("cd molStructure/%s/FE\n"%(name[i]))
				f.write("mkdir steep-%d\n"%(j))
				f.write("cd steep-%d\n"%(j))
				f.write("echo 'staring steep minimization'\n")
				f.write("\n")
				f.write("gmx grompp -f ../../../../Freesol/lambda-%d/steep-%d.mdp -c ../../%s_pack.pdb\
				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2\n"%(j,j,name[i]))
				#f.write("mpirun -nt %d gmx mdrun -v -deffnm steep\n"%(nt))
				f.write("gmx mdrun -v -nt $NSLOTS -deffnm steep\n")
				f.write("wait\n")
				f.write("cd ../\n")
				f.write("\n")
				f.write("mkdir lbfg-%d\n"%(j))
				f.write("cd lbfg-%d\n"%(j))
				f.write("echo 'starting L-BFGS minimization'\n")
				f.write("\n")
				f.write("gmx grompp -f ../../../../Freesol/lambda-%d/lbfgs-%d.mdp -c ../steep-%d/steep.gro\
				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2\n"%(j,j,j))
				#f.write("mpirun -nt 1 gmx mdrun -v -deffnm lbfgs\n")
				f.write("gmx mdrun -nt 1 -v -deffnm lbfgs\n")
				f.write("wait\n")
				f.write("cd ../\n")
				f.write("\n")
				f.write("mkdir nvt-%d\n"%(j))
				f.write("cd nvt-%d\n"%(j))
				f.write("echo 'Starting constant volume equilibration'\n")
				f.write("\n")
				f.write("gmx grompp -f ../../../../Freesol/lambda-%d/nvt-%d.mdp -c ../lbfg-%d/lbfgs.gro\
				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2\n"%(j,j,j))
				#f.write("mpirun -nt %d gmx mdrun -v -deffnm nvt\n" %(nt))
				f.write("gmx mdrun -v -nt $NSLOTS -deffnm nvt\n")
				f.write("wait\n")
				f.write("cd ../\n")
				f.write("\n")
				f.write("mkdir npt-%d\n"%(j))
				f.write("cd npt-%d\n"%(j))
				f.write("echo 'Starting constant pressure equilibration'\n")
				f.write("\n")
				f.write("gmx grompp -f ../../../../Freesol/lambda-%d/npt-%d.mdp -c ../nvt-%d/nvt.gro\
				-t ../nvt-%d/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2\n"%(j,j,j,j))
				#f.write("mpirun -nt %d gmx mdrun -v -deffnm npt\n" %(nt))
				f.write("gmx mdrun -v -nt $NSLOTS -deffnm npt\n")
				f.write("wait\n")
				f.write("cd ../\n")
				f.write("\n")
				f.write("mkdir pd-%d\n"%(j))
				f.write("cd pd-%d\n"%(j))
				f.write("echo 'Starting production MD simulation'\n")
				f.write("\n")
				f.write("gmx grompp -f ../../../../Freesol/lambda-%d/pd-%d.mdp -c ../npt-%d/npt.gro\
				-t ../npt-%d/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2\n"%(j,j,j,j))
				#f.write("mpirun -nt %d gmx mdrun -v -deffnm pd\n" %(nt))
				f.write("gmx mdrun -v -nt $NSLOTS -deffnm pd\n")
				f.write("wait\n")
				f.write("cd ../\n")
				f.write("\n")
				f.write("cd ../../../\n")
				f.write("\n")
