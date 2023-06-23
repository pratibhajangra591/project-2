# Dynamics of the various shells of DM halo around a PBH of mass, M_PBH=100 
# solar mass at different initial radii


import numpy as np
from scipy.integrate import odeint
import math
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib import rcParams

plt.rcParams.update({
  "text.usetex": True,
  "font.family": "serif"
})
 

π = np.pi
Ω_cdm = 0.85
G = 6.67e-11                # in units of m^3⋅kg^−1⋅s^−2
M_solar = 1.989e30        # in units of kg
A = 1.495e11                 # Astronomical unit
pc = 3.085e16               # in unist of meter. 
yr = 3.154e7              # in units of seconds
Λ = 2.036e-35                # in units of s⁻2
c = 3e8                     # in units of ms⁻¹


a_eq = 2.9374e-4
t_eq = 1.5923e+12   #in units of s, as per Eq.(4) of https://arxiv.org/pdf/1107.2025.pdf with                       a_eq=2.9374e-4
ρ_eq = 2.1548e-16    # in units of kg m^-3
ρ_meq = ρ_eq/2
σ_eq = 0.005

#M_H_eq = 3e12 * (t/1e-23)      # in units of solar mass, Horizon mass at z_eq, 
                                 #  https://arxiv.org/pdf/1706.10288.pdf 


t_m = 13.78e9 * yr            #in units of yrs corresponding to t_0=13.78Gyr    
t_0 = t_m        # age of the Universe today
t_1 = t_eq       # time at the end of radiation domination which is the time of MRE.
t_2 = 0.5 * t_0  # time at the end of matter-domination



z_th = 173        # Thermalization redshift
z_fr = 3.4e4      # DM freeze out redshift
z_rec = 1089.90   # recombination redshift. 
z_eq = ((1/a_eq)-1)      # matter-radiation equality, z ≈ 3403



h = 0.67
ρ_c0 = 1.9e-26 * (h**2)     #in units of kgm⁻³
Ω_r0 = 9.4e-5
Ω_m0 = 0.32
ρ_r0 = Ω_r0 * ρ_c0
ρ_m0 = Ω_m0 * ρ_c0
H_0 = np.sqrt((8 * π * G * ρ_c0)/3)




def ρ(z): #Density of the Universe in units of kgm⁻³
    ρ_r = ρ_r0 * ((1 + z )**4)
    ρ_b = ρ_b0 * ((1 + z )**3)
    if  z < z_th:              #after thermalization
        return ρ_b0
    elif z_th < z < z_rec:     #post recombination era 
        return ρ_b
    elif z_rec < z < z_eq:     #pre recombination era
        return  (ρ_r +  ρ_b)    
    elif z_eq < z < z_fr:      #post-DM freeze - out accretion
        return (ρ_r +  ρ_b)
    else :                     #pre-DM freeze - out accretion
        return (ρ_r +  ρ_b)

    

def c_s(z): #Speed of sound in units of ms⁻¹.
    if  z < z_th :               #after thermalization
        return 15 * (1 + z)
    elif z_th < z < z_rec:       #post recombination era 
        return 1.9e2 * np.sqrt(1 + z)
    elif z_rec < z < z_eq:       #pre recombination era 
        ρ_r = ρ_r0 * ((1 + z)**4)
        ρ_b = ρ_b0 * ((1 + z)**3)
        return  (c/np.sqrt(3)) * np.sqrt((4 * ρ_r)/(4 * ρ_r + 3 * ρ_b))    
    elif z_eq < z < z_fr:        #post-DM freeze - out accretion
        return c/np.sqrt(3)
    else :                       #pre-DM freeze - out accretion
        return c/np.sqrt(3)



    
#Redshift at which formation of PBH takes place as matter domination (Eq. 20) in https://arxiv.org/pdf/1706.10288.pdf   
def z_pbh(m):
    t_i = ( G * m)/( (c**3))
    gamma = 1/3
    return (((t_i/t_0)**(-1/2)) - 1) * gamma



# if v << c_s, Bondi radius of PBH of mass m_pbh
def r_B(m, z):
    return (2 * G * m)/(c_s(z)**2)


