"""Set global physical variables used in the irs.py script.
"""

# Numerical precision
ORDER = 3
N_F = int(1e4)
N_Y = int(8e4)

# Gyro lines
f = 430e6  # 1/s ~ radar frequency
n_e = 2e10  # 1/m^3 ~ electron number density
B = 3.5e-5  # T ~ magnetic field strength (towards Earth)
M = 16  # amu ~ ion mass
T_e = 200  # K ~ electron temperature
T_i = 200  # K ~ ion temperature
aspect = 135  # degree ~ radar pointing direction to magnetic field line
