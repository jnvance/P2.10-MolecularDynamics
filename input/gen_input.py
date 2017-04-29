#!/usr/bin/env python3

import numpy as np

base = "# starting structure: \n\
inputfile crystal.xyz \n\
# final structure: \n\
outputfile output_{1}.xyz \n\
# temperature (Langevin thermostat) \n\
temperature {0} \n\
# timestep \n\
tstep 0.005 \n\
# friction coefficient (Langevin thermostat) \n\
friction 1.0 \n\
# cutoff for forces \n\
forcecutoff 2.5 \n\
# cutoff for neighbor lists \n\
listcutoff  3.0 \n\
# total number of steps \n\
nstep 2000 \n\
# stride for writing trajectory \n\
nconfig 10 trajectory_{1}.xyz \n\
# stride for writing energy file \n\
nstat   10 energies_{1}.dat \n\
# seed \n\
idum {2}\n\
"

istr = lambda k,n=3: str(k).zfill(n)

nrep = 4
minTemp = 0.7
maxTemp = 0.9

temps = np.linspace(minTemp,maxTemp,nrep)

for idtemp in range(nrep):

    filename = "in_{0}.params".format(istr(idtemp))
    text = base.format(temps[idtemp],istr(idtemp),-(idtemp+nrep))
    
    in_file = open(filename,"w")
    in_file.write(text)
    in_file.close()

    print("Saved parameters for temperature {0:12.8} to: {1}".format(temps[idtemp],filename))

