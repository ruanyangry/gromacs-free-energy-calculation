#!/bin/bash
#$ -S /bin/sh
#$ -N test
#$ -j y
#$ -o .
#$ -e .
#$ -cwd
#$ -q all.q
#$ -masterq all.q@node9
#$ -pe thread 16-16
source ~/.bashrc
hash -r
export PATH=/export/home/ry/gromacs-5.0/bin:PATH


# Author: Ruan Yang
# Email: ruanyang_njut@163.com

echo 'Staring state 0'
cd ../molStructure/mobley_210639
mkdir FE
cd FE
mkdir steep-0
cd steep-0
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-0/steep-0.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-0
cd lbfg-0
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-0/lbfgs-0.mdp -c ../steep-0/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-0
cd nvt-0
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-0/nvt-0.mdp -c ../lbfg-0/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-0
cd npt-0
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-0/npt-0.mdp -c ../nvt-0/nvt.gro				-t ../nvt-0/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-0
cd pd-0
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-0/pd-0.mdp -c ../npt-0/npt.gro				-t ../npt-0/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 1'
cd molStructure/mobley_210639/FE
mkdir steep-1
cd steep-1
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-1/steep-1.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-1
cd lbfg-1
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-1/lbfgs-1.mdp -c ../steep-1/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-1
cd nvt-1
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-1/nvt-1.mdp -c ../lbfg-1/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-1
cd npt-1
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-1/npt-1.mdp -c ../nvt-1/nvt.gro				-t ../nvt-1/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-1
cd pd-1
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-1/pd-1.mdp -c ../npt-1/npt.gro				-t ../npt-1/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 2'
cd molStructure/mobley_210639/FE
mkdir steep-2
cd steep-2
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-2/steep-2.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-2
cd lbfg-2
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-2/lbfgs-2.mdp -c ../steep-2/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-2
cd nvt-2
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-2/nvt-2.mdp -c ../lbfg-2/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-2
cd npt-2
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-2/npt-2.mdp -c ../nvt-2/nvt.gro				-t ../nvt-2/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-2
cd pd-2
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-2/pd-2.mdp -c ../npt-2/npt.gro				-t ../npt-2/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 3'
cd molStructure/mobley_210639/FE
mkdir steep-3
cd steep-3
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-3/steep-3.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-3
cd lbfg-3
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-3/lbfgs-3.mdp -c ../steep-3/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-3
cd nvt-3
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-3/nvt-3.mdp -c ../lbfg-3/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-3
cd npt-3
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-3/npt-3.mdp -c ../nvt-3/nvt.gro				-t ../nvt-3/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-3
cd pd-3
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-3/pd-3.mdp -c ../npt-3/npt.gro				-t ../npt-3/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 4'
cd molStructure/mobley_210639/FE
mkdir steep-4
cd steep-4
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-4/steep-4.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-4
cd lbfg-4
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-4/lbfgs-4.mdp -c ../steep-4/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-4
cd nvt-4
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-4/nvt-4.mdp -c ../lbfg-4/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-4
cd npt-4
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-4/npt-4.mdp -c ../nvt-4/nvt.gro				-t ../nvt-4/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-4
cd pd-4
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-4/pd-4.mdp -c ../npt-4/npt.gro				-t ../npt-4/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 5'
cd molStructure/mobley_210639/FE
mkdir steep-5
cd steep-5
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-5/steep-5.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-5
cd lbfg-5
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-5/lbfgs-5.mdp -c ../steep-5/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-5
cd nvt-5
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-5/nvt-5.mdp -c ../lbfg-5/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-5
cd npt-5
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-5/npt-5.mdp -c ../nvt-5/nvt.gro				-t ../nvt-5/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-5
cd pd-5
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-5/pd-5.mdp -c ../npt-5/npt.gro				-t ../npt-5/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 6'
cd molStructure/mobley_210639/FE
mkdir steep-6
cd steep-6
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-6/steep-6.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-6
cd lbfg-6
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-6/lbfgs-6.mdp -c ../steep-6/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-6
cd nvt-6
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-6/nvt-6.mdp -c ../lbfg-6/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-6
cd npt-6
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-6/npt-6.mdp -c ../nvt-6/nvt.gro				-t ../nvt-6/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-6
cd pd-6
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-6/pd-6.mdp -c ../npt-6/npt.gro				-t ../npt-6/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 7'
cd molStructure/mobley_210639/FE
mkdir steep-7
cd steep-7
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-7/steep-7.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-7
cd lbfg-7
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-7/lbfgs-7.mdp -c ../steep-7/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-7
cd nvt-7
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-7/nvt-7.mdp -c ../lbfg-7/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-7
cd npt-7
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-7/npt-7.mdp -c ../nvt-7/nvt.gro				-t ../nvt-7/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-7
cd pd-7
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-7/pd-7.mdp -c ../npt-7/npt.gro				-t ../npt-7/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 8'
cd molStructure/mobley_210639/FE
mkdir steep-8
cd steep-8
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-8/steep-8.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-8
cd lbfg-8
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-8/lbfgs-8.mdp -c ../steep-8/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-8
cd nvt-8
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-8/nvt-8.mdp -c ../lbfg-8/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-8
cd npt-8
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-8/npt-8.mdp -c ../nvt-8/nvt.gro				-t ../nvt-8/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-8
cd pd-8
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-8/pd-8.mdp -c ../npt-8/npt.gro				-t ../npt-8/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 9'
cd molStructure/mobley_210639/FE
mkdir steep-9
cd steep-9
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-9/steep-9.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-9
cd lbfg-9
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-9/lbfgs-9.mdp -c ../steep-9/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-9
cd nvt-9
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-9/nvt-9.mdp -c ../lbfg-9/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-9
cd npt-9
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-9/npt-9.mdp -c ../nvt-9/nvt.gro				-t ../nvt-9/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-9
cd pd-9
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-9/pd-9.mdp -c ../npt-9/npt.gro				-t ../npt-9/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 10'
cd molStructure/mobley_210639/FE
mkdir steep-10
cd steep-10
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-10/steep-10.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-10
cd lbfg-10
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-10/lbfgs-10.mdp -c ../steep-10/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-10
cd nvt-10
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-10/nvt-10.mdp -c ../lbfg-10/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-10
cd npt-10
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-10/npt-10.mdp -c ../nvt-10/nvt.gro				-t ../nvt-10/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-10
cd pd-10
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-10/pd-10.mdp -c ../npt-10/npt.gro				-t ../npt-10/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 11'
cd molStructure/mobley_210639/FE
mkdir steep-11
cd steep-11
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-11/steep-11.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-11
cd lbfg-11
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-11/lbfgs-11.mdp -c ../steep-11/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-11
cd nvt-11
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-11/nvt-11.mdp -c ../lbfg-11/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-11
cd npt-11
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-11/npt-11.mdp -c ../nvt-11/nvt.gro				-t ../nvt-11/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-11
cd pd-11
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-11/pd-11.mdp -c ../npt-11/npt.gro				-t ../npt-11/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 12'
cd molStructure/mobley_210639/FE
mkdir steep-12
cd steep-12
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-12/steep-12.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-12
cd lbfg-12
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-12/lbfgs-12.mdp -c ../steep-12/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-12
cd nvt-12
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-12/nvt-12.mdp -c ../lbfg-12/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-12
cd npt-12
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-12/npt-12.mdp -c ../nvt-12/nvt.gro				-t ../nvt-12/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-12
cd pd-12
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-12/pd-12.mdp -c ../npt-12/npt.gro				-t ../npt-12/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 13'
cd molStructure/mobley_210639/FE
mkdir steep-13
cd steep-13
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-13/steep-13.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-13
cd lbfg-13
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-13/lbfgs-13.mdp -c ../steep-13/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-13
cd nvt-13
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-13/nvt-13.mdp -c ../lbfg-13/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-13
cd npt-13
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-13/npt-13.mdp -c ../nvt-13/nvt.gro				-t ../nvt-13/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-13
cd pd-13
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-13/pd-13.mdp -c ../npt-13/npt.gro				-t ../npt-13/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 14'
cd molStructure/mobley_210639/FE
mkdir steep-14
cd steep-14
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-14/steep-14.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-14
cd lbfg-14
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-14/lbfgs-14.mdp -c ../steep-14/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-14
cd nvt-14
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-14/nvt-14.mdp -c ../lbfg-14/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-14
cd npt-14
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-14/npt-14.mdp -c ../nvt-14/nvt.gro				-t ../nvt-14/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-14
cd pd-14
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-14/pd-14.mdp -c ../npt-14/npt.gro				-t ../npt-14/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 15'
cd molStructure/mobley_210639/FE
mkdir steep-15
cd steep-15
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-15/steep-15.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-15
cd lbfg-15
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-15/lbfgs-15.mdp -c ../steep-15/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-15
cd nvt-15
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-15/nvt-15.mdp -c ../lbfg-15/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-15
cd npt-15
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-15/npt-15.mdp -c ../nvt-15/nvt.gro				-t ../nvt-15/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-15
cd pd-15
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-15/pd-15.mdp -c ../npt-15/npt.gro				-t ../npt-15/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 16'
cd molStructure/mobley_210639/FE
mkdir steep-16
cd steep-16
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-16/steep-16.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-16
cd lbfg-16
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-16/lbfgs-16.mdp -c ../steep-16/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-16
cd nvt-16
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-16/nvt-16.mdp -c ../lbfg-16/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-16
cd npt-16
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-16/npt-16.mdp -c ../nvt-16/nvt.gro				-t ../nvt-16/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-16
cd pd-16
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-16/pd-16.mdp -c ../npt-16/npt.gro				-t ../npt-16/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 17'
cd molStructure/mobley_210639/FE
mkdir steep-17
cd steep-17
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-17/steep-17.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-17
cd lbfg-17
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-17/lbfgs-17.mdp -c ../steep-17/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-17
cd nvt-17
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-17/nvt-17.mdp -c ../lbfg-17/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-17
cd npt-17
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-17/npt-17.mdp -c ../nvt-17/nvt.gro				-t ../nvt-17/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-17
cd pd-17
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-17/pd-17.mdp -c ../npt-17/npt.gro				-t ../npt-17/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 18'
cd molStructure/mobley_210639/FE
mkdir steep-18
cd steep-18
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-18/steep-18.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-18
cd lbfg-18
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-18/lbfgs-18.mdp -c ../steep-18/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-18
cd nvt-18
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-18/nvt-18.mdp -c ../lbfg-18/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-18
cd npt-18
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-18/npt-18.mdp -c ../nvt-18/nvt.gro				-t ../nvt-18/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-18
cd pd-18
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-18/pd-18.mdp -c ../npt-18/npt.gro				-t ../npt-18/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 19'
cd molStructure/mobley_210639/FE
mkdir steep-19
cd steep-19
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-19/steep-19.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-19
cd lbfg-19
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-19/lbfgs-19.mdp -c ../steep-19/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-19
cd nvt-19
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-19/nvt-19.mdp -c ../lbfg-19/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-19
cd npt-19
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-19/npt-19.mdp -c ../nvt-19/nvt.gro				-t ../nvt-19/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-19
cd pd-19
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-19/pd-19.mdp -c ../npt-19/npt.gro				-t ../npt-19/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

echo 'Staring state 20'
cd molStructure/mobley_210639/FE
mkdir steep-20
cd steep-20
echo 'staring steep minimization'

gmx grompp -f ../../../../Freesol/lambda-20/steep-20.mdp -c ../../mobley_210639_pack.pdb				-p ../../TST_GMXsol.top -o steep.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm steep
wait
cd ../

mkdir lbfg-20
cd lbfg-20
echo 'starting L-BFGS minimization'

gmx grompp -f ../../../../Freesol/lambda-20/lbfgs-20.mdp -c ../steep-20/steep.gro				-p ../../TST_GMXsol.top -o lbfgs.tpr -maxwarn 2
gmx mdrun -nt 1 -v -deffnm lbfgs
wait
cd ../

mkdir nvt-20
cd nvt-20
echo 'Starting constant volume equilibration'

gmx grompp -f ../../../../Freesol/lambda-20/nvt-20.mdp -c ../lbfg-20/lbfgs.gro				-p ../../TST_GMXsol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm nvt
wait
cd ../

mkdir npt-20
cd npt-20
echo 'Starting constant pressure equilibration'

gmx grompp -f ../../../../Freesol/lambda-20/npt-20.mdp -c ../nvt-20/nvt.gro				-t ../nvt-20/nvt.cpt -p ../../TST_GMXsol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm npt
wait
cd ../

mkdir pd-20
cd pd-20
echo 'Starting production MD simulation'

gmx grompp -f ../../../../Freesol/lambda-20/pd-20.mdp -c ../npt-20/npt.gro				-t ../npt-20/npt.cpt -p ../../TST_GMXsol.top -o pd.tpr -maxwarn 2
gmx mdrun -v -nt $NSLOTS -deffnm pd
wait
cd ../

cd ../../../

