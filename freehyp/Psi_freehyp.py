#!/mnt/home/alexander/bin/python
# Calculates the wave function psi_freehyp from wellenfkt.py (free particle WF with energy = Va / R_start + Vb for R >= R_start, 0 for R < R_start)
# Usage: Psi_freehyp.py [R_start [phase_over_pi]]       ; provide R_start in Angstrom

# Importing modules
import numpy as np
from sys import argv, path
path.append('/mnt/home/alexander/eldest')
import sciconv as sc
import wellenfkt as wf
from fractions import Fraction

# Initialization
R_start = 6.            # In Angstrom
phase_over_pi = 0.      # Additional phase shift (in units of pi)

# Read in possible cmd line args
if (len(argv) > 1):
    R_start = float(argv[1])
if (len(argv) > 2):
    phase_over_pi = float(Fraction(argv[2]))

print("R_start / Angstrom =", R_start)
print("phase / pi =", phase_over_pi)

R_start_au = sc.angstrom_to_bohr(R_start)
phase = np.pi * phase_over_pi


# Set potential curve from which free-particle energy is chosen: V = Va / R_start + Vb
Va = sc.ev_to_hartree(sc.angstrom_to_bohr(1.2))
Vb = sc.ev_to_hartree(40.)
# Reduced mass used as particle mass
red_mass = wf.red_mass_au(20.1797,20.1797)


# Write complex WF into first file, imag part (starts with 0 at R = R_start for phase = 0) into second file
with open("psi_freehyp/norm_psi_freehyp_{}_{}.txt".format(R_start,phase_over_pi), "w") as f:
    with open("psi_freehyp/imag_norm_psi_freehyp_{}_{}.txt".format(R_start,phase_over_pi), "w") as g:
        f.write("# R / Angstrom     psi_freehyp(R_au,Va={},Vb={},red_mass={},R_start_au=sc.angstrom_to_bohr({}),phase=pi*{})\n".format(Va, Vb, red_mass, R_start, phase_over_pi))
        g.write("# R / Angstrom     imag(psi_freehyp(R_au,Va={},Vb={},red_mass={},R_start_au=sc.angstrom_to_bohr({}),phase=pi*{}))\n".format(Va, Vb, red_mass, R_start, phase_over_pi))
         
        for R in np.arange(1.5,10.,0.0001):
            print("R / Angstrom = ", R)
            psi = wf.psi_freehyp(sc.angstrom_to_bohr(R),Va,Vb,red_mass,R_start_au,phase)
             
            f.write(str(R))
            f.write(" ")
            f.write(str(psi))
            f.write("\n")
             
            g.write(str(R))
            g.write(" ")
            g.write(str(np.imag(psi)))
            g.write("\n")
