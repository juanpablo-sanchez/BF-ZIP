#!/bin/bash
outputdir=/scratch/juan/h2_micro_ZIP
cd $outputdir

rm -rf test_linear
mkdir test_linear
cd test_linear

length_sample=10000
burnin_sample=1000
lap=10



cat <<EOF > run_h2_simul.slurm
#!/bin/bash
#SBATCH -a 1-10
#SBATCH -J h2_rat
#SBATCH -e $outputdir/test_linear/FB_simul_h2.err
#SBATCH -o $outputdir/test_linear/FB_simul_h2.log
#SBATCH -p std
#SBATCH -t 0-10:00
#SBATCH -N 1                       #En un solo nodo
#SBATCH -n 1                       #Un solo proceso
#SBATCH -c 1                       #5 Cpus para ese proceso

pos=\$SLURM_ARRAY_TASK_ID

module load apps/R/3.5.1
module load toolchains/intel_mkl_ompi
module load tools/intel/comp_xe_2015
ulimit -s unlimited

export OMP_NUM_THREADS=\${SLURM_CPUS_PER_TASK}

mkdir rep_\${pos}
cd rep_\${pos}

cp /home/juan/FB/baycom/test_ZIP_ratios_INV/simul_raw.R        ./

echo ' DATAFILE '         > IMparameters
echo ' dat '             >> IMparameters
echo 'NUMBER_OF_TRAITS'  >> IMparameters
echo '       1'          >> IMparameters
echo 'NUMBER_OF_EFFECTS' >> IMparameters
echo '           6'      >> IMparameters
echo 'OBSERVATION(S)'    >> IMparameters
echo '   8 '            >> IMparameters
echo 'WEIGHT(S)'         >> IMparameters
echo ' '                 >> IMparameters
echo 'EFFECTS: POSITIONS_IN_DATAFILE NUMBER_OF_LEVELS TYPE_OF_EFFECT[EFFECT NESTED] ' >> IMparameters
echo '    2   3             cross   #TRAT         '  >> IMparameters
echo '    3   3             cross   #BANDA        '  >> IMparameters
echo '    4   2             cross   #MEDIDA       '  >> IMparameters
echo '    5  50             cross   #camada       '  >> IMparameters
echo '    6 100             cross   #jaula        '  >> IMparameters
echo '    7 620             cross   #aditivo      '  >> IMparameters
echo 'RANDOM_RESIDUAL VALUES  ' >> IMparameters
echo ' 1.0  '          >> IMparameters
echo ' RANDOM_GROUP  ' >> IMparameters
echo '     4  '        >> IMparameters
echo ' RANDOM_TYPE  '  >> IMparameters
echo ' diagonal  '     >> IMparameters
echo ' FILE  '         >> IMparameters
echo '  '              >> IMparameters
echo '(CO)VARIANCES  ' >> IMparameters
echo ' 1.0  '          >> IMparameters
echo ' RANDOM_GROUP  ' >> IMparameters
echo '     5  '        >> IMparameters
echo ' RANDOM_TYPE  '  >> IMparameters
echo ' diagonal  '     >> IMparameters
echo ' FILE  '         >> IMparameters
echo '  '              >> IMparameters
echo '(CO)VARIANCES  ' >> IMparameters
echo '  1.0  '         >> IMparameters
echo ' RANDOM_GROUP  ' >> IMparameters
echo '     6  '        >> IMparameters
echo ' RANDOM_TYPE  '  >> IMparameters
echo ' add_animal  '   >> IMparameters
echo ' FILE  '         >> IMparameters
echo ' ped '         >> IMparameters
echo '(CO)VARIANCES  ' >> IMparameters
echo '   1.0  '        >> IMparameters
echo 'OPTION solution mean' >> IMparameters


R  --vanilla  < simul_raw.R > out_simul 2> err_simul

echo -e " IMparameters \n $length_sample \n $burnin_sample   \n $lap " |  /home/juan/BLUPF90/f90-05/bin/gibbsf90test  > out_gibbs      2> err_gibbs
echo -e " IMparameters \n 1 \n $lap      \n 0                        " |  /home/juan/BLUPF90/f90-05/bin/postgibbsf90  > out_post_gibbs 2> err_post_gibbs


EOF



sbatch run_h2_simul.slurm

