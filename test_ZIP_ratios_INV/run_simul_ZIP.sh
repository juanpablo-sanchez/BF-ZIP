#!/bin/bash
outputdir=/home/juan/GitHub/BF-ZIP/test_ZIP_ratios_INV/
cd $outputdir

rm -rf test_ZIP
mkdir test_ZIP
cd test_ZIP

length_sample=1000
burnin_sample=100
lap=1 # siempre 1


for i in `seq 1 3`
do

echo 'REPLICA : ' ${i}

pos=$i


mkdir rep_${pos}
cd rep_${pos}

cp $outputdir/simul_raw.R        ./

echo 'dat                 #fichero de datos (caracter, factores fijos, factores aleatorios, animal=pedigree) ' > parametros_baycom
echo 'ped                 #fichero de pedigree ' >> parametros_baycom
echo '600                 #numero de datos en fichero datos ' >> parametros_baycom
echo '620                 #numero de animales en fichero pedigree' >> parametros_baycom
echo '5                   #numero de factores fijos y aleatorios (excepto el genético). Tienen que venir primero los fijos y luego los aleatorios (como en el fichero de datos)' >> parametros_baycom
echo '3 3 2 50 100        #niveles de los fijos y aleatorios   ' >> parametros_baycom
echo '0.0 0.0 0.0 0.2 0.1 #ratios de varianza (si son fijos=0)' >> parametros_baycom
echo '10.0                #varianza fenotipica' >> parametros_baycom
echo '0.3                 #heredabilidad' >> parametros_baycom
echo '0.5                 #p=prob de venir de Possion' >> parametros_baycom
echo '2.0                 #DT normal para general proposal en el MH para log(lambda)' >> parametros_baycom
echo $length_sample '     #nround     ' >> parametros_baycom
echo $burnin_sample '     #burnin' >> parametros_baycom
echo '5                   #n_procesadores' >> parametros_baycom


R  --vanilla  < simul_raw.R > out_simul 2> err_simul   
cat parametros_baycom | /home/juan/GitHub/BF-ZIP/bin/ZIP_FB_h2_ratios_INV   > out_FB 2> err_FB &

cd ..

done

