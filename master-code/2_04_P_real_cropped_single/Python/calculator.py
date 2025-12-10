import math

c = 299792458  # m/s
m_p = 1.67262192e-27
q_e =  1.602e-19
B = 3

energy_mev = 100


# Convert MeV to eV
energy_ev = energy_mev * 1e6  # 1 MeV = 1e6 eV

# Convert eV to Joules
energy_joules = energy_ev * q_e

# formula: r = mv / qB
# v needs to be calculated first !
# gamma = 1 + (energy_mev/931.494) OR calculated with 
# E = gamma * m * c**2 <-> gamma = E/m*(c**2)

gamma = 1 + (energy_joules/(m_p * (c**2)))
#gamma = 1 + (energy_mev/931.494) used in MATLAB code
beta = math.sqrt(1-(1/gamma**2))
v = beta * c

r = (m_p * v)/(q_e * B)
r_cm = r * 100

print(r_cm)
