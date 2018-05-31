# _*_ coding: utf-8 _*_

'''
@author: Ruan Yang
Created on 2018.5.27
This code used to prepare gromacs itp or top file.
'''

def topfreesol(filename,natoms):
	'''
	itpfreesol: used to generated gromacs free energy calculated .itp file.
	            make each atom chrage == 0.0
	filename: The name of .itp or .top file.
	natoms: total atoms in .itp file.
	reference: http://www.bevanlab.biochem.vt.edu/Pages/Personal/justin/gmx-tutorials/free_energy
	'''
	with open(filename,"r") as f:
		lines=f.readlines()
	names=filename.split('.')
	with open("%ssol.%s"%(names[0],names[1]),"w") as f:
		i=0
		for line in lines:
			f.write(line)
			i +=1
			words=line.strip().split()
			if len(words) == 3 and words[1] == "atoms":
				#print("i= %d"%(i))
				break
		f.write("%s"%(lines[i]))
		
		for line in lines[i+1:(i+1+natoms)]:
			words=line.strip().split()
			words[6]="0.000"
			line="   ".join([word for word in words])
			f.write("%s\n"%(line))
			
		for line in lines[(i+natoms+1):]:
			f.write(line)
