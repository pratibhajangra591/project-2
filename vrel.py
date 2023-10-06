import numpy as np
from scipy.interpolate import interp1d
from scipy import interpolate


z_arr, c_s_arr = np.loadtxt("c_s.txt", unpack=True)
z_rel_arr, v_rel_arr = np.loadtxt("vrel.txt", unpack=True)
z_arrA, v_effA_arr = np.loadtxt("veffA.txt", unpack=True)
z_arrB, v_effB_arr = np.loadtxt("veffB.txt", unpack=True)



x = z_arr
y = c_s_arr
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
        
