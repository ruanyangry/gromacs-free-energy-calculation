# _*_ coding: utf-8 _*_

'''
@author: Ruan Yang
Created on 2018.5.7
Anaconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/archive/
Openbabel-2.4.1: http://openbabel.org/wiki/Main_Page
Antechamber: http://ambermd.org/antechamber/ac.html
Gaussian: http://gaussian.com/
Acpype.py: https://github.com/t-/acpype
amber2lammpsry.py: Can be find in lammps/tools/ directory. I modified the input/output file name.
PACKMOL: http://www.ime.unicamp.br/~martinez/packmol/userguide.shtml
'''

import os
import numpy as np       

# Ellement Mol mass
# Reference: https://ptable.com/

molmass={'H':1.008,
         'HE':4.003,
         'LI':6.941,
         'BE':9.012,
         'B':10.81,
         'C':12.01,
         'N':14.01,
         'O':16.00,
         'F':19.00,
         'NE':20.18,
         'NA':22.99,
         'MG':24.31,
         'AL':26.98,
         'SI':28.09,
         'P':30.97,
         'S':32.06,
         'CL':35.45,
         'AR':39.95,
         'K':39.10,
         'CA':40.08,
         'SC':44.96,
         'TI':47.87,
         'V':50.94,
         'CR':52.00,
         'MN':54.94,
         'FE':55.85,
         'CO':58.93,
         'NI':58.69,
         'CU':63.55,
         'ZN':65.58,
         'GA':69.72,
         'GE':72.63,
         'AS':74.92,
         'SE':78.96,
         'BR':79.90,
         'KR':83.80,
         'RB':85.47,
         'SR':87.62,
         'Y':88.91,
         'ZR':91.22,
         'NB':92.91,
         'MO':95.96,
         'TC':98.00,
         'RU':101.1,
         'RH':102.9,
         'PD':106.4,
         'AG':107.9,
         'CD':112.4,
         'IN':114.8,
         'SN':118.7,
         'SB':121.8,
         'TE':127.6,
         'I':126.9,
         'XE':131.3,
         'CS':133,
         'BA':137.3,
         'LA':139.0,
         'CE':140.0,
         'PR':141.0,
         'ND':144.0,
         'PM':145,
         'SM':150.5,
         'EU':152,
         'GD':157,
         'TB':159,
         'DY':162.5,
         'HO':165,
         'ER':167,
         'TM':169,
         'YB':173,
         'LU':175,
         'HF':178.5,
         'TA':181,
         'W':184,
         'RE':186,
         'OS':190.0,
         'IR':192,
         'PT':195,
         'AU':197,
         'HG':200.6,
         'TL':204.5,
         'PB':207,
         'BI':209,
         'PO':209,
         'AT':210,
         'RN':222,
         'FR':223,
         'RA':226,
         'AC':227,
         'TH':232,
         'PA':231,
         'U':238,
         'NP':237,
         'PU':244,
         'AM':243,
         'CM':247,
         'BK':247,
         'CF':251,
         'ES':252,
         'FM':257,
         'MD':258,
         'NO':259,
         'LR':260,
         'RF':261,
         'DB':262,
         'SG':263,
         'BH':264,
         'HS':265,
         'MT':266,
         'DS':269,
         'RG':272,
         'CN':277,
         'NH':286,
         'FL':289,
         'MC':289,
         'LV':294,
         'TS':294,
         'OG':294} 
  
# Output Element table to check
    
#for k,v in molmass.items():
#	print("%s %.2f"%(k,v))

# Define load molecule smi format file function
# File format: #Compound ID   SMILES
# Reference:Machine Learning of Partial Charges Derived From High-Quality 
# Quantum-Mechanical Calculations.
# DOI: 10.1021/acs.jcim.7b00663

def loadData(filename):
	'''
	filename: input molecule smi
	'''
	name=[]
	smi=[]
	with open(filename,'r') as f:
		lines=f.readlines()
		for line in lines:
			words=line.strip().split()
			if words[1] == 'ID':
				continue
			else:
				name.append(words[0])
				smi.append(words[1])
	return name,smi

# Use openbabel python library pybel convert smi to mol2

def Smi2mol2(molname,molsmi,molformat='mol2',platform='windows'):
	'''
	molname : molecule name
	molsmi : molecule smiles format
	molformat default value : 'mol2'
	pltform : windows or linux, default windows
	'''
	if platform=='windows':
		import pybel
		mol=pybel.readstring("smi","%s"%(molsmi))
		mol.make3D()
		mol.write("%s"%(molformat),"%s%s%s"%(molname,'.',molformat))
	if platform=='linux':
		import openbabel
		pass
	
# Through read mol2 file get the molecule mol mass

def readmol2(filename):
	'''
	filename : molecule_name.mol2
	return value: 
	M: molecule mol mass
	tcharge: total charge of molecules (charge value rstore in .mol2)
	element: get the unique element in molecules
	atoms: total atoms in one molecule
	'''
	M=0.0
	tcharge=0.0
	element=[]
	with open(filename,'r') as f:
		lines=f.readlines()
		for line in lines:
			words=line.split()
			if(len(words)==5):
				atoms=int(words[0])
				bonds=int(words[1])
		for line in lines[7:7+atoms:1]:
			words=line.split()
			for k,v in molmass.items():
				if words[1]==k:
					M += v
					tcharge += float(words[-1])
					element.append(words[1])
		element=set(element)
	return M,tcharge,element,atoms
	
# function nmolecule used to calculated the number of molecules
# packed into the box. The default size of the box 0. 0. 0. 40. 40. 40.
# Ameldeo Avogadro constant = 6.02×10^23/mol
# density=mass/volume
# mass = n*M
# number = (density*volume)/(M*Av)

def nmolecule(Molmass,box,freesol="off"):
	'''
	Molmass : single molecule molar mass
	box : box size (xlo,ylo,zho,xhi,yhi,zhi)
	return value: nmol the number of molecules in pack box, density = 1000 kg/m3
	freesol: generated 1 organic + nmol water system. default='off'
	'''
	Molwater=16.00+1.008*2
	Avc=6.02*(10**(23))
	density=1000.0   #g/m3
	volume=(box[3]-box[0])*(box[4]-box[1])*(box[5]-box[2])*(10**(-27)) # m3
	if freesol=="on":
		nmol=round(((density*volume*Avc)-Molmass)/Molwater)
	else:
		nmol=round((density*volume)/(Molmass*(1/Avc)))-1
	return nmol
	
# First write PACKMOL input file, just pack single molecule
# into big box

def Packmolinputfile(inpname,filetype,number,box,freesol="off",method='GaussianAntechamber'):
	'''
	inpname : The name of the .inp file ,default = molecule name defined in .smi
	filetype : .xyz or .pdb
	number : number of molecules packed into simulation box
	box : box size (xlo,ylo,zho,xhi,yhi,zhi)
	freesol: True ,defaule=false
	method : Generated molecule topology file method, default = GaussianAntechamber
	'''
	with open('%s.inp'%(inpname),'w') as f:
		f.write('# This file generated by writeinputfile.py. Author: Ruan Yang\n')
		f.write(' \n')
		f.write('tolerance 2.0\n')
		f.write('filetype %s\n'%(filetype))
		f.write('output %s%s%s%s\n'%(inpname,'_pack','.',filetype))
		f.write('add_box_sides 1.5\n')
		f.write(' \n')
		if(method=="acpypeantechamber"):
			f.write('structure %s%s%s%s\n'%(inpname,'_NEW','.',filetype))
		else:
			f.write('structure %s%s%s\n'%(inpname,'.',filetype))
		if freesol=="on":
			f.write('  number %d\n'%(1))
			f.write('  fixed %.2f %.2f %.2f 0. 0. 0.\n'%(box[3],box[4],box[5]))
			f.write('end structure\n')
			f.write('\n')
			f.write('structure %s%s%s\n'%("water",'.',filetype))
			f.write('  number %d\n'%(number))
			f.write('  inside box %.2f %.2f %.2f %.2f %.2f %.2f\n'%(box[0],box[1],box[2],\
			box[3],box[4],box[5]))
			f.write('end structure\n')
		else:
			f.write('  number %d\n'%(number))
			f.write('  inside box %.2f %.2f %.2f %.2f %.2f %.2f\n'%(box[0],box[1],box[2],\
			box[3],box[4],box[5]))
			f.write('end structure\n')
		
# Just used acpype.py generated gromacs,AMBER,CHARMM input file.
# for gromacs there have two force field: AMBER and OPLS/AA
# usage: python acpype.py -i mobley_5857.mol2 -b mobley_5857 -c user \
# -n 0 -o gmx -a gaff2 -d
# -i: default input molecule file format .mol2
# mobley_5857: molecule name
# -c: use charge define in .mol2
# -n: net charge 0
# -o: write gromacs topology file
# -a: atom type (default=gaff, I want use gaff2)
# -d: for debugging purposes, keep any temporary file created

def acpypeantechamber(molname,nmol,data='off'):
	'''
	molname : molecule name
	nmol : number of molecules packed into simulation box
	data : generated lammps data file, default='off'
	'''
	with open('acpypeantechamber.sh','w') as f:
		f.write('#!/usr/bin/env python3 \n')
		f.write(' \n')
		f.write('python3 acpype.py -i %s.mol2 -b %s -c user -n 0 -o gmx \
		-a gaff2 -d\n'%(molname,molname))
		if data=='on':
			f.write('python2 amber2lammpsry.py %s\n'%(molname))
		f.write('packmol < *.inp\n')
		print('sed -i " s/TST              1 /TST  %d/g " TST_GMX.top\n'%(nmol),file=f)

# Define function write Gaussian and antechamber input file.
# Please reference: http://ambermd.org/tutorials/ForceField.php
# We fitting Energy landspace to get the RESP charges.
# We used .prmtop and inpcrd file, only output GAFF topology.

def GaussianAntechamber(molname,nmol,basisset='#HF/6-31G*',data='off',freesol="off"):
	'''
	molname : molecule name
	nmol : number of molecules packed into simulation box
	basisset : the basis set in Gaussian, default='#HF/6-31G*',(B3LYP/6-311G(d,p))
	data : generated lammps data file, default='off'
	'''
	with open('GaussianAntechamber.sh','w') as f:
		f.write('#!/bin/bash\n')
		f.write('echo "#-----------------------------------------------------------#"\n')
		f.write('# Reference: http://ambermd.org/tutorials/ForceField.php\n')
		f.write('echo "#-----------------------------------------------------------#"\n')
		f.write(' \n')
		# ''%%'' == %
		f.write('echo "#-----------------------------------------------------------#"\n')
		f.write('echo "Step 1: prepare Gaussian input file .com"\n')
		f.write('echo "#-----------------------------------------------------------#"\n')
		f.write('antechamber -i %s.mol2 -fi mol2 -o %s.com -fo \
		gcrt -pf y -gm "%%mem=4096MB" -gn "%%nproc=8" -nc 0 -gk \
		"%s SCF=tight Test Pop=MK iop(6/33=2) iop(6/42=6)"\n'%(molname,molname,basisset))
		f.write(' \n')
		f.write('echo "#-----------------------------------------------------------#"\n')
		f.write('echo "Step 2: run Gaussian jobs."\n')
		f.write('echo "#-----------------------------------------------------------#"\n')
		f.write('nohup g09 %s.com > %s.log &\n'%(molname,molname))
		f.write('wait \n')
		f.write(' \n')
		#f.write('mv %s.log %s.out \n'%(molname,molname))
		f.write('echo "#-----------------------------------------------------------#"\n')
		f.write('echo "Step 3: fitting Energy landspace get the RESP charge."\n')
		f.write('echo "#-----------------------------------------------------------#"\n')
		f.write('antechamber -i %s.log -fi gout -o %s.mol2 -fo mol2 -c resp -rn TST \n'%(molname,molname))
		f.write(' \n')
		f.write('echo "#-----------------------------------------------------------#"\n')
		f.write('echo "Step 4: Get the atom deficiency force field parameters."\n')
		f.write('echo "#-----------------------------------------------------------#"\n')
		f.write('parmchk2 -i %s.mol2 -f mol2 -o %s.frcmod \n'%(molname,molname))
		f.write(' \n')
		f.write('echo "#-----------------------------------------------------------#"\n')
		f.write('echo "Step 5: Using tleap get amber input file .inpcrd and .prmtop."\n')
		f.write('echo "#-----------------------------------------------------------#"\n')
		f.write('tleap -f leap.in \n')
		f.write(' \n')
		# acpype.py used python2 interpreter
		f.write('echo "#-----------------------------------------------------------#"\n')
		f.write('echo "Step 6: Using acpype.py get GROMACS topology and coordinate file."\n')
		f.write('echo "#-----------------------------------------------------------#"\n')
		f.write('python2 ../../acpype.py -p %s.prmtop -x %s.inpcrd -d\n'%(molname,molname)) 
		f.write(' \n')
		f.write('echo "#-----------------------------------------------------------#"\n')
		f.write('echo "Step 7: Using packmol get the packed coordinate."\n')
		f.write('echo "#-----------------------------------------------------------#"\n')
		f.write('packmol < *.inp\n')
		f.write(' \n')
		f.write('echo "#-----------------------------------------------------------#"\n')
		f.write('echo "Step 8: Adjust the .top file molecules number information."\n')
		f.write('echo "#-----------------------------------------------------------#"\n')
		if freesol=="on":
			print('sed -i "s/TST              1 /TST  %d/g" TST_GMX.top\n'%(1),file=f)
			f.write('echo  SOL    %d >> TST_GMX.top\n'%(nmol))
			print('sed -i \'/defaults/i\#include "amber03.ff/forcefield.itp"\' TST_GMX.top',file=f) 
			#print('sed -i \'/moleculetype/i\#include "amber03.ff/tip3p.itp"\' TST_GMX.top',file=f)
			print('sed -i \'/moleculetype/i\#include "amber03.ff/tip3p.itp"\\n\' TST_GMX.top',file=f)
			print('sed -i \'/defaults/d\' TST_GMX.top',file=f)
			print('sed -i \'/; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ/d\' TST_GMX.top',file=f)
			print('sed -i \'/1               2               yes             0.5     0.8333/d\' TST_GMX.top',file=f)
		else:
			print('sed -i "s/TST              1 /TST  %d/g" TST_GMX.top\n'%(nmol),file=f)
		f.write(' \n')
		if data=='on':
			f.write('echo "#-----------------------------------------------------------#"\n')
			f.write('echo "Step 9: Generated lammps data file."\n')
			f.write('echo "#-----------------------------------------------------------#"\n')
			f.write('python2 ../../amber2lammpsry.py %s \n'%(molname[i]))
			f.write(' \n')
	with open('leap.in','w') as f:
		f.write('source leaprc.ff14SB\n')
		f.write('source leaprc.gaff\n')
		f.write('loadamberparams %s.frcmod\n'%(molname))
		f.write('%s=loadmol2 %s.mol2\n'%(molname,molname))
		f.write('check %s\n'%(molname))
		f.write('saveamberparm %s %s.prmtop %s.inpcrd\n'%(molname,molname,molname))
		f.write('savepdb %s %s.pdb\n'%(molname,molname))
		f.write('quit\n')




