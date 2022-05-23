# BF-ZIP



PARAMETER FILE FOR ZIP MODEL (ZIP_FB_h2_ratios_INV)


d_rec               # data file (trait, fixed factors, random factor, animal=pedigree,Poisson_average)
p_rec               # pedigree file
500                 # number of rows in data file
100                 # number of rows in pedigree file numero de animales en fichero pedigree
5                   # number of fixed and random (except the genetic efect) First fixed them random
5 6 2 196 189       # levels of fixed and random factors
0.0 0.0 0.0 0.1 0.1 # starting values for ratios of variance (when fixed = 0)
5.0                 # starting value for Phenotypic variance (log_ind_poisson_average)
0.10                # starting value for heritability
0.5                 # starting value of the probability of Zero in the ZIP model 
2.0                 # Standard deviation of the Metropolis Hasting steps for the update of hte log_lambda parameters in the ZIP model 
10000               # number of rounds
1000                # burning
1                   # number of processors to use in the matrix multiplications and inversions
0                   # set the heritabilityto zero = 1, estimate the heritability = 0


PARAMETER FILE FOR NORMAL MODEL (FB_h2_ratios_INV)


d_rec               # data file (trait, fixed factors, random factor, animal=pedigree)
p_rec               # pedigree file
500                 # number of rows in data file
100                 # number of rows in pedigree file numero de animales en fichero pedigree
5                   # number of fixed and random (except the genetic efect) First fixed them random
5 6 2 196 189       # levels of fixed and random factors
0.0 0.0 0.0 0.1 0.1 # starting values for ratios of variance (when fixed = 0)
5.0                 # starting value for Phenotypic variance 
0.10                # starting value for heritability
10000               # number of rounds
1000                # burning
1                   # number of processors to use in the matrix multiplications and inversions
0                   # set the heritabilityto zero = 1, estimate the heritability = 0







