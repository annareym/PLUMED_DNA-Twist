# PLUMED_DNA-Twist
DNA twisting code to be used as a new collective variable in [PLUMED library](http://www.plumed.org) for [MD simulations](http://en.wikipedia.org/wiki/Molecular_dynamics).

## Authors
**Anna Reymer**[1], **Krystyna Zakrzewska**[2], **Richard Lavery**[2].

[1] University of Gothenburg, Department of Chemistry and Molecular Biology, Box 462, 40530 Göteborg, Sweden

[2] Institut de Biologie et Chimie des Protéines, Univ. Lyon I/CNRS UMR 5086, 7 passage du Vercors, Lyon 69367, France

## Description
This code can be used to monitor or control the value of total twist between any two base pair levels in a DNA fragment, alone or in complex with other molecules. The code can used throuhg PLUMED with a variety of MD packages, e.g. Amber or GROMACS. To use the code, it has to be placed in, for example, `plumed-2.2.1/src/colvar` folder, PLUMED has to be recompiled and linked to the corresponding MD software. For PLUMED installation please check: http://plumed.github.io/doc-v2.1/user-doc/html/_installation.html

## How to use
DNA twisting code works analogously to any collective variable (colvar) implemented in PLUMED.
TWIST colvar monitors or controls the value of total twist between any chosen base pair levels *i* and *j*.
As in imput to TWIST colvar, provide 48 atom numbers from bases that are restrained: (3×4 bases ×3 atoms: C1',N1/N9 and C6/C8 depending whether purine or pyrimidine), plus 6×2 auxiliary atoms, force constant, desired value of total twist,and how many turns the restricted DNA fragment has. The energy penalty will be added to the potential energy functional: `E_tw=0.5*k*(tw0-tw)^2`.

The numbers that represent DNA bases (desriptors) are taken from bases *i-1*, *i*, *i+1* and *j-1*, *j*, *j+1* from the Watson DNA strand (5'->3'), and from the Crick DNA strand (3'->5'). The auxiliary atoms are used as a support for a correct implementation of the periodic boundary condition, and have to be selected to be evenly destributed along the restricted DNA fragment. 

### Example of plumed.dat input file:
```
tw: TWIST ATOMS=44,42,45,74,72,75,106,104,107,167,169,261,263,358,360,423,421,424,455,453,456,487,485,488,5
79,577,580,612,610,613,644,642,645,704,706,801,803,896,898,960,958,961,993,991,994,1025,1023,1026 N_TURNS=1.0
tw_r: RESTRAINT ARG=tw KAPPA=0.25 AT=358.8
PRINT STRIDE=250 ARG=tw,tw_r.bias FILE=Twist_358.8
```

First 9 atom numbers represent bases *i-1*, *i*, *i+1* from the Watson strand, followed by 6 auxiliary atoms, follwed by 9 atoms from bases *j-1*, *j*, *j+1* from the Watson strand, follwed by 9 atoms from bases *j-1*, *j*, *j+1* from the Crick strand, followed by 6 auxiliary atoms, and finally the last 9 atoms represent bases *i-1*, *i*, *i+1* from the Crick strand.
