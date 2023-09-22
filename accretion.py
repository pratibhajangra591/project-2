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
G = 4.4959e-15            #in units of M☉^-1 pc^3 yr^-2
c = 0.3068                #in units of pc yr^-1
pc = 3.0857e16            # in meters
yr = 3.154e7              # in units of seconds
t_0 = 13.78e9             #in units of yrs as the age of the Universe today, t_0=13.78Gyr
M_solar = 1.989e30        # in units of kg.         



σ_eq = 0.005
σ_eq = 0.005
a_eq = 2.9374e-4
ρ_eq = 3.1808e3           #in units of M☉ pc^-3 with ρ_eq = 2.1548e-16 kgm^-3 
ρ_meq = ρ_eq/2
t_eq = 1.5923e+12/yr      #as per Eq.(4) of https://arxiv.org/pdf/1107.2025.pdf with                                       a_eq=2.9374e-4


t_1 = t_eq              # time at the end of radiation domination which is the time of MRE.
t_2 = 0.5 * t_0         # time at the end of matter-domination



z_th = 173               # Thermalization redshift
z_fr = 3.4e4             # DM freeze out redshift
z_rec = 1089.90          # recombination redshift. 
z_eq = ((1/a_eq)-1)      # matter-radiation equality, z ≈ 3403



h = 0.67
ρ_c0 = 2.80659e-7 * (h**2)    #in units of M☉ pc^-3 with ρ_c0= 1.9e-26 * (h**2) kgm^-3   
Ω_r0 = 9.4e-5
Ω_m0 = 0.32
ρ_r0 = Ω_r0 * ρ_c0
ρ_m0 = Ω_m0 * ρ_c0
H_0 = np.sqrt((8 * π * G * ρ_c0)/3)

η_acc = 0.1  #Accretion efficiency
 


def ρ(z): #Density of the Universe in units of kgm⁻³
    ρ_r = ρ_r0 * ((1 + z)**4)
    ρ_m = ρ_m0 * ((1 + z)**3)
    if  z < z_th:              #after thermalization
        return ρ_m0
    elif z_th < z < z_rec:     #post recombination era 
        return ρ_m
    elif z_rec < z < z_eq:     #pre recombination era
        return  (ρ_r +  ρ_m)    
    elif z_eq < z < z_fr:      #post-DM freeze - out accretion
        return (ρ_r +  ρ_m)
    else :                     #pre-DM freeze - out accretion
        return (ρ_r +  ρ_m)

    

def c_s(z): #Speed of sound in units of ms⁻¹.
    if  z < z_th :               #after thermalization
        return 15 * (1 + z) * (yr/pc)
    elif z_th < z < z_rec:       #post recombination era 
        return 1.9e2 * np.sqrt(1 + z) * (yr/pc)
    elif z_rec < z < z_eq:       #pre recombination era 
        ρ_r = ρ_r0 * ((1 + z)**4)
        ρ_m = ρ_m0 * ((1 + z)**3)
        return  (c/np.sqrt(3)) * np.sqrt((4 * ρ_r)/(4 * ρ_r + 3 * ρ_m))    
    elif z_eq < z < z_fr:        #post-DM freeze - out accretion
        return c/np.sqrt(3)
    else :                       #pre-DM freeze - out accretion
        return c/np.sqrt(3)



#Redshift at which formation of PBH takes place as matter domination (Eq. 20) in https://arxiv.org/pdf/1706.10288.pdf  
# Here, γ is the ratio between the PBH mass and the horizon mass.
# t = np.sqrt(3/(4 * np.pi * G * ρ_eq)) * ((s**2)/2)  with s = a/a_eq.
def z_pbh(m, γ):
    def t_i(γ):
        return (G * m)/(γ * (c**3))
    return (np.sqrt(1/(2 * t_i(γ))) * ((3/(4 * π * G * ρ_eq))**(1/4)) * (1 + z_eq)) - 1 




# if v << c_s, Bondi radius of PBH of mass m_pbh
def r_B(m, z):
    return (2 * G * m)/(c_s(z)**2)
