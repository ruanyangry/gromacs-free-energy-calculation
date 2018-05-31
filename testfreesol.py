# _*_ coding: utf-8 _*_

import GenerateFF as GF
import freesolsh as FSH
import freesolmdp as FSP
import gmxtop     as GMP
import subprocess
import numpy as np
from scipy.interpolate import interp1d
import os

# mol2 file location
name,smi=GF.loadData('molecules.smi')
name=np.array(name)
smi=np.array(smi)

os.mkdir('molStructure')
os.chdir('molStructure')

print('#-------------------------------------------------------#')
print('                        Start                            ')
print('#-------------------------------------------------------#')
print(' \n')



print('#-------------------------------------------------------#')
print('       Generated Gaussian and PACKMOL input file         ')
print('#-------------------------------------------------------#')
print(' \n')

for i in range(len(name)):
	# mkdir cd 
	os.mkdir('%s'%(name[i]))
	os.chdir('%s'%(name[i]))
	GF.Smi2mol2(name[i],smi[i])
	Molmass,tcharge,element,natoms=GF.readmol2(name[i]+'.mol2')
	with open('%s-information.txt'%(name[i]),'w') as f:
		f.write('%s molecule charge = %.2f'%(name[i],tcharge))
		f.write('\n')
		f.write('Unique element in %s\n'%(name[i]))
		f.write('%s'%(element))
	box=[0.,0.,0.,40.,40.,40.]
	nmol=GF.nmolecule(Molmass,box,freesol="on")
	methods=["acpypeantechamber","GaussianAntechamber"]
	method="GaussianAntechamber"
	if(method==methods[0]):
		GF.Packmolinputfile(name[i],'pdb',nmol,box,method,freesol="on")
		GF.acpypeantechamber(name[i],nmol)
	else:
		GF.Packmolinputfile(name[i],'pdb',nmol,box,freesol="on")
		GF.GaussianAntechamber(name[i],nmol,freesol="on")
	
	os.chdir("../../")
	os.system("cp water.pdb molStructure/%s"%(name[i]))
	os.chdir("cd molStructure/%s"%(name[i]))
	os.system("dos2unix *.sh *.inp")
	os.system("chmod 700 *.sh")
	os.system("./*.sh")
	
	GMP.topfreesol("TST_GMX.top",natoms)
	
	os.chdir('../')
	
os.chdir('../')
	
print('#-------------------------------------------------------#')
print('            Generate batch commit task scripts           ')
print('#-------------------------------------------------------#')
print(' \n')

os.mkdir("sh")
os.chdir("sh")
FSH.fsolsh(21,name,16)	
os.chdir("../")	

print('#-------------------------------------------------------#')
print('   Generate GROMACS Free energy calculated mdp file      ')
print('#-------------------------------------------------------#')
print(' \n')

names=["steep","lbfgs","nvt","npt","pd"]
methods=["steep","l-bfgs","sd","sd","sd"]
nsteps=[5000,5000,50000,50000,500000]
lambda_state=[i for i in range(20)]
vdw_lambdas=[0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00]
coul_lambdas=[0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00]
bonded_lambdas=[0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00]
restraint_lambdas=[0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00]
mass_lambdas=[0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00]
temperature_lambdas=[0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00]

os.mkdir("Freesol")
path=os.getcwd()
os.chdir(path+"/Freesol")

FSP.Femdp(names,methods,nsteps,lambda_state,vdw_lambdas,coul_lambdas,bonded_lambdas,\
restraint_lambdas,mass_lambdas,temperature_lambdas)

os.chdir("../")


