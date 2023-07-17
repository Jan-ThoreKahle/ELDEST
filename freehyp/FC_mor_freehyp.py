#!/mnt/home/alexander/bin/python
# Calculates the overlap integral S = <psi_n|psi_freehyp> and the Franck-Condon factor |S|**2 from wellenfkt.py 
# Usage: Psi_freehyp.py [R_start [phase_over_pi [n [mode [outfile]]]]]       ; provide R_start in Angstrom; mode = "R", "p" or other; keep a default by passing "d"
# Alexander Riegel, July 2023.

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
n = 0

# Read in possible cmd line args. Argument "d" skips the respective parameter and keeps its default value (hence "d").
if (len(argv) > 1):
    if (argv[1] != "d"): R_start = float(argv[1])
if (len(argv) > 2):
    if (argv[2] != "d"): phase_over_pi = float(Fraction(argv[2]))
if (len(argv) > 3):
    if (argv[3] != "d"): n = int(argv[3])

# Set free-particle potential curve 
Va = sc.ev_to_hartree(sc.angstrom_to_bohr(1.2))
Vb = sc.ev_to_hartree(40.)
# Set Morse potential curve - provide alpha in a.u., Req in Angstrom, De in eV
alpha = 15.3994
Req = 6.
De = 0.5
# Reduced mass
red_mass = wf.red_mass_au(20.1797,20.1797)

# FC integration boundaries
R_min = sc.angstrom_to_bohr(1.5)
R_max = sc.angstrom_to_bohr(30.)


R_start_au = sc.angstrom_to_bohr(R_start)
phi = np.pi * phase_over_pi
Req_au = sc.angstrom_to_bohr(Req)
De_au = sc.ev_to_hartree(De)


####################
# For looping over R_start or phase_over_pi, choose "R" or "p" and the outfile.
# Default (i. e., neither "R" nor "p" is chosen) is single-point calculation at fixed R_start or phase, no outfile needed.
mode = "s"
outfile = "fc-mor-freehyp-scan.out"

if (len(argv) > 4):
    mode = str(argv[4])
if (len(argv) > 5):
    outfile = str(argv[5])
####################

if (mode == "R"):
    with open(outfile, "w") as f:
        f.write("#S = wf.FCmor_freehyp(n,alpha,Req_au,De_au,red_mass,Va,Vb,R_start_au,R_min,R_max,phase,limit)\n")
        f.write(f"#  = wf.FCmor_freehyp({n},{alpha},{Req_au},{De_au},{red_mass},{Va},{Vb},R_start_au,{R_min},{R_max},phase={phi},limit=5000)\n")
        f.write(f"#R_start {'S':<55} |S|**2\n")

        for R in np.arange(0.0010, 0.0300, 0.0001):
            S = wf.FCmor_freehyp(n,alpha,Req_au,De_au,red_mass,Va,Vb,sc.angstrom_to_bohr(R),R_min,R_max,phase=phi,limit=5000)
            f.write(f"{R:<7.7} {S:< 20.20E}  {np.abs(S)**2}\n")
        for R in np.arange(0.0300, 10.000, 0.0010):
            S = wf.FCmor_freehyp(n,alpha,Req_au,De_au,red_mass,Va,Vb,sc.angstrom_to_bohr(R),R_min,R_max,phase=phi,limit=5000)
            f.write(f"{R:<7.7} {S:< 20.20E} {np.abs(S)**2}\n")
elif (mode == "p"):
    with open(outfile, "w") as f:
        f.write("#S = wf.FCmor_freehyp(n,alpha,Req_au,De_au,red_mass,Va,Vb,R_start_au,R_min,R_max,phase,limit)\n")
        f.write(f"#  = wf.FCmor_freehyp({n},{alpha},{Req_au},{De_au},{red_mass},{Va},{Vb},{R_start_au},{R_min},{R_max},phase=phase_over_pi*pi,limit=5000)\n")
        f.write(f"#phase_over_pi {'S':<55} |S|**2\n")

        for p in np.arange(0., 2., 0.01):
            S = wf.FCmor_freehyp(n,alpha,Req_au,De_au,red_mass,Va,Vb,R_start_au,R_min,R_max,phase=p*np.pi,limit=5000)
            f.write(f"{p:<13.13} {S:< 20.20E}  {np.abs(S)**2}\n")
else:
    S = wf.FCmor_freehyp(n,alpha,Req_au,De_au,red_mass,Va,Vb,R_start_au,R_min,R_max,phase=phi,limit=5000)
    print(f"R_start phase_over_pi n  {'S':<54} |S|**2")
    print(f"{R_start:<7.7}", f"{phase_over_pi:<13.13}", n, f"{S:< 20.20E}", np.abs(S)**2)
