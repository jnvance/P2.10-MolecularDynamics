#!/usr/bin/env python3

import numpy as np
import os

base = "# starting structure:\n\
inputfile crystal.xyz\n\
# final structure:\n\
outputfile output_{1}.xyz\n\
# temperature (Langevin thermostat)\n\
temperature {0}\n\
# timestep\n\
tstep 0.005\n\
# friction coefficient (Langevin thermostat)\n\
friction 1.0\n\
# cutoff for forces\n\
forcecutoff 2.5\n\
# cutoff for neighbor lists\n\
listcutoff  3.0\n\
# total number of steps\n\
nstep 2000`\n\
# stride for writing trajectory\n\
nconfig 10 trajectory_{1}.xyz\n\
# stride for writing energy file\n\
nstat   10 energies_{1}.dat\n\
# seed\n\
idum {2}\n\
# stride for performing parallel tempering\n\
exchangestride 30\n\
"

for file in os.listdir():
    if (".params" in file): os.remove(file)

istr = lambda k,n=3: str(k).zfill(n)

nrep = 16
minTemp = 1.0
maxTemp = 1.5

# geometric spacing
r = np.power(maxTemp/minTemp,1.0/(nrep-1))
temps = np.power(r,range(nrep))

# # linear spacing
# temps = np.linspace(minTemp,maxTemp,nrep)

for idtemp in range(nrep):

    filename = "in_{0}.params".format(istr(idtemp))
    text = base.format(temps[idtemp],istr(idtemp),-(idtemp+nrep))

    in_file = open(filename,"w")
    in_file.write(text)
    in_file.close()

    print("Saved parameters for temperature {0:12.8} to: {1}".format(temps[idtemp],filename))

