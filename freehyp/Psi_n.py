#!/mnt/home/alexander/bin/python
# Calculates the wave function psi_n from wellenfkt.py (Morse potential nuclear WF)
# Usage: Psi_n.py [n]

# Importing modules
import numpy as np
from sys import argv, path
path.append('/mnt/home/alexander/eldest')
import sciconv as sc
import wellenfkt as wf

# Initialization
n = 0

# Read in possible cmd line args
if (len(argv) > 1):
    n = int(argv[1])

print("n =", n)


# Set Morse potential curve - provide alpha in a.u., Req in Angstrom, De in eV
alpha = 15.3994
Req = 6.
De = 0.5
# Reduced mass
red_mass = wf.red_mass_au(20.1797,20.1797)

Req_au = sc.angstrom_to_bohr(Req)
De_au = sc.ev_to_hartree(De)

# Write WF into file
with open("psi_n_{}.txt".format(n), "w") as f:
        f.write("# R / Angstrom     psi_n(R_au,n={},alpha={},Req_au=sc.angstrom_to_bohr({}),red_mass={},De=sc.ev_to_hartree({}))\n".format(n, alpha, Req, red_mass, De))
         
        for R in np.arange(1.5,30.,0.0001):
            print("R / Angstrom = ", R)
            f.write(str(R))
            f.write(" ")
            f.write(str(wf.psi_n(sc.angstrom_to_bohr(R),n,alpha,Req_au,red_mass,De_au)))
            f.write("\n")

print("n =", n)
