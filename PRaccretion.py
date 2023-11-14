# Dynamics of the various shells of DM halo around a PBH of mass, M_PBH=100 
# solar mass at different initial radii


import numpy as np
from scipy.integrate import odeint
import math
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
from scipy.interpolate import interp1d
from scipy import interpolate


 


    
    
π = np.pi
Ω_cdm = 0.85
G = 4.4959e-15            #in units of M☉^-1 pc^3 yr^-2
c = 0.3068                #in units of pc yr^-1
pc = 3.0857e16            # in meters
yr = 3.154e7              # in units of seconds
t_0 = 13.78e9             #in units of yrs as the age of the Universe today, 
M_solar = 1.989e30        # in units of kg.         



σ_eq = 0.005
a_eq = 2.9374e-4
ρ_eq = 3.1808e3           #in units of M☉ pc^-3 with ρ_eq = 2.1548e-16 kgm^-3 
ρ_meq = ρ_eq/2
t_eq = 1.5923e+12/yr      #as per Eq.(4) of https://arxiv.org/pdf/1107.2025.pdf with                                       a_eq=2.9374e-4




z_fr = 3.4e4             # DM freeze out redshift
z_rec = 1089.90          # recombination redshift. 
z_eq = ((1/a_eq)-1)      # matter-radiation equality, z ≈ 3403
z_th = 173               # Thermalization redshift
z_dec = 130              # Redshift of the decoupling of matter and radiations




#Redshift at which formation of PBH takes place as matter domination (Eq. 20) in Philippa S. Cole et al. https://arxiv.org/pdf/1706.10288.pdf  
# Here, γ is the ratio between the PBH mass and the horizon mass.
# t = np.sqrt(3/(4 * np.pi * G * ρ_eq)) * ((s**2)/2)  with s = a/a_eq.
def z_pbh_Cole(m, γ):
    def t_i(γ):
        return (G * m)/(γ * (c**3))
    return (np.sqrt(1/(2 * t_i(γ))) * ((3/(4 * π * G * ρ_eq))**(1/4)) * (1 + z_eq)) - 1 



h = 0.67
ρ_c0 = 2.80659e-7 * (h**2)    #in units of M☉ pc^-3 with ρ_c0= 1.9e-26 * (h**2)kgm^-3  
Ω_r0 = 9.4e-5
Ω_m0 = 0.32
ρ_r0 = Ω_r0 * ρ_c0
ρ_m0 = Ω_m0 * ρ_c0
H_0 = np.sqrt((8 * π * G * ρ_c0)/3)

η_acc = 0.1  #Accretion efficiency
 


def ρ_Jared(z): #Density of the Universe in units of kgm⁻³
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

    
# Speed of sound in units of pcyr⁻¹ as per Fig.(B.5) of 
# Jared R. Rice et al. https://arxiv.org/abs/1702.08069v2
def c_s_Jared(z): 
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

    
    
    





# Reference paper Ricotti et al. https://arxiv.org/pdf/0709.0524.pdf  

#Redshift at which formation of PBH takes place as per
# (Eq. 1) of Ricotti et al. https://arxiv.org/pdf/0709.0524.pdf  
# Here, f_Hor is the ratio between the PBH mass and the horizon mass.
def z_pbh_Ricotti(m, f_Hor):
    M_H = 3.1e16 #Horizon mass at z_eq in units of solar mass
    return (1 + z_eq) * (((1/f_Hor) * (m/M_H))**(-1/2)) - 1 


# Mass of the DM halo accreted around the PBH
def M_halo(z, m):
    return 3 * m * (((1+z)/1000)**(-1))
 

# Radius of DM halo around the PBH which is 1/3rd of the turnaroud radius of the halo.
def r_halo(z, m):  # in units pf pc
    return 0.019 * M_halo(z, m) * (((1+z)/1000)**(-1))



#in units of M_solar yr⁻¹, from Eq.(12) of https://arxiv.org/pdf/2003.12589.pdf
def M_dot_Edd(z, m):  
    return 2e-3 * 1e-6 * m





#Digitized and extrapolated values of c_s, v_rel, v_effA, v_effB
z_cs_arr, cs_arr = np.loadtxt("c_s.txt", unpack=True)
z_rel_arr, v_rel_arr = np.loadtxt("vrel.txt", unpack=True)
z_arrA, v_effA_arr = np.loadtxt("veffA.txt", unpack=True)
z_arrB, v_effB_arr = np.loadtxt("veffB.txt", unpack=True)



x = z_cs_arr
y = cs_arr
f = interpolate.interp1d(x, y , fill_value = "extrapolate")
def c_s_extrapolation(variable): #in units of kms⁻¹
    if variable < 10:
        m = 0.9823
        b = -4.1193/np.log(10)
        return   (variable**m) * (10**b)
    elif  10 <= variable <= 1000:
        return f(variable)
    else:
        m =  0.5485
        b = -2.0570/np.log(10)
        return  (variable**m) * (10**b)
        




x_rel = z_rel_arr
y_rel = v_rel_arr
f_rel = interpolate.interp1d(x_rel, y_rel, fill_value = "extrapolate")
def vrel_extrapolation(variable):  #in units of kms⁻¹
    if variable < 10:
        m = 0.4712
        b = -0.3375/np.log(10)
        return  (variable**m) * (10**b)
    elif 10 <= variable <= 1000:
        return f_rel(variable)
    else:
        m = -0.4348
        b = 4.3650/np.log(10)
        return  (variable**m) * (10**b)
 

    

    
x_A = z_arrA
y_A = v_effA_arr
f_A = interpolate.interp1d(x_A, y_A, fill_value = "extrapolate")
def veffA_extrapolation(variable):  #in units of kms⁻¹
    if variable < 10:
        m = 0.7313
        b = -1.9190/np.log(10)
        return (variable**m) * (10**b)
    elif 10 <= variable <= 1000:
        return f_A(variable)
    else:
        m = 0.3039
        b = -0.1793/np.log(10)
        return (variable**m) * (10**b)
      

    

    
x_B = z_arrB
y_B = v_effB_arr
f_B = interpolate.interp1d(x_B, y_B, fill_value = "extrapolate")
def veffB_extrapolation(variable):
    if  variable < 10:
        m = 0.5308
        b = -0.6707/np.log(10)
        return (variable**m) * (10**b)
    elif 10 <= variable <= 1000:
        return f_B(variable)
    else:
        m = 0.3185
        b = -0.2794/np.log(10)
        return (variable**m) * (10**b)

  
    
def c_s_Ricotti(z):
    β = 1.72
    z_dec = 130
    return  5.74e3 * np.sqrt((1+z)/1000) * (((((1+z_dec)/(1 + z))**β) + 1)**(-1/(2 *              β))) * (yr/pc)


def v_rel_digitized_Ricotti(z):
    return  vrel_extrapolation(z) * 1e3 *  (yr/pc)


def v_eff_Ricotti(z): #Here, v_eff = v_eff,A, https://arxiv.org/pdf/0709.0524.pdf
    def Mach_number(z):
        return v_rel_digitized_Ricotti(z)/c_s_Ricotti(z)
    if Mach_number(z) > 1:
            return  c_s_Ricotti(z) * (((16/np.sqrt(2 * np.pi)) *                                           (Mach_number(z)**3))**(1/6))  
    else:
            return  c_s_Ricotti(z) * np.sqrt(1 + (Mach_number(z)**2))





def ρ(z):  # Density of the surrounding gas
    m_H = 1.67e-27  # Hydrogen mass in units of kg
    def n_gas(z):   # in units of pc⁻³
        return 2e8 * (((1+z)/1000)**3) * (pc**3)
    return (m_H/M_solar) * n_gas(z)

    
    
#Reference research paper Vivian Poulin, Pasquale D. Serpico et al. https://arxiv.org/pdf/1707.04206.pdf  
def v_L_Serpico(z):
    return np.minimum(1, ((1+z)/1000)) * 30e3 * (yr/pc)

def c_s_Serpico(z):
    return 6e3 * (yr/pc) * np.sqrt((1+z)/1000)


def v_eff_Serpico(z):
    return np.sqrt(c_s_Serpico(z)**2 + v_L_Serpico(z)**2)


#Conversion factor of dt/dz for dm/dz = dm/dt. dt/dz
def dt_dz(z):
    x = ((Ω_r0 * ((1 + z)**6)) + (Ω_m0 * ((1 + z)**5)))
    return  - np.sqrt(3/(8 * π * G * ρ_c0)) * (x**(-1/2))



def dt3_dz(z):
    first_term = - ((1 + z_eq)/((1+z)**2))
    second_term = - (1/2) * (1/np.sqrt(((1+z_eq)/(1+z)) + 1)) * ((1 + z_eq)/((1+z)**2))
    return   np.sqrt(3/(4 * π * G * ρ_eq)) * ((2/3) * first_term * np.sqrt(((1+z_eq)/(1+z)) + 1) + (2/3) * (((1+z_eq)/(1+z)) - 2) * second_term)