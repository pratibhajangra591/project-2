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





#Redshift at which formation of PBH takes place as matter domination (Eq. 20) in https://arxiv.org/pdf/1706.10288.pdf  
# Here, γ is the ratio between the PBH mass and the horizon mass.
# t = np.sqrt(3/(4 * np.pi * G * ρ_eq)) * ((s**2)/2)  with s = a/a_eq.
def z_pbh(m, γ):
    def t_i(γ):
        return (G * m)/(γ * (c**3))
    return (np.sqrt(1/(2 * t_i(γ))) * ((3/(4 * π * G * ρ_eq))**(1/4)) * (1 + z_eq)) - 1 
    


    
    
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
f_rel = interpolate.interp1d(x_rel, y_rel , fill_value = "extrapolate")
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
f_A = interpolate.interp1d(x_A, y_A , fill_value = "extrapolate")
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
        

    
def v_eff(z): #Here, v_eff = v_eff,A
    def c_s(z):
        β = 1.72
        z_dec = 130
        return  5.74e3 * np.sqrt((1+z)/1000) * (((((1+z_dec)/(1 + z))**β) + 1)**(-1/(2 * β)))            * (yr/pc)
    def v_rel_digitized(z):
        return  vrel_extrapolation(z) * 1e3 *  (yr/pc)
    def Mach_number(z):
        return v_rel_digitized(z)/c_s(z)
    if Mach_number(z) > 1:
            return  c_s(z) * (((16/np.sqrt(2 * np.pi)) * (Mach_number(z)**3))**(1/6))  
    else:
            return  c_s(z) * np.sqrt(1 + (Mach_number(z)**2))


        
        
def M_halo(z, m):
    return 3 * m * (((1+z)/1000)**(-1))
        
    