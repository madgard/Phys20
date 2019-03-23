from math import log, log10, sqrt
from scipy.constants import c
import numpy as np
import matplotlib.pyplot as mpl
from scipy.optimize import curve_fit
from scipy.integrate import quad
def distance_modulus_to_luminosity_distance(distance, uncertainty):
    d_parsecs = 10**(distance/5 + 1)
    dd_parsecs = 0.461*np.multiply(d_parsecs,uncertainty)
    d_Mpc = d_parsecs/10**6
    dd_Mpc = dd_parsecs/10**6
    return d_Mpc, dd_Mpc
# def luminosity_distance_to_distance_modulus(distance):
#     d_parsecs = distance*10**6
#     distance_modulus = [5*(log10(d)-1) for d in d_parsecs]
#     return distance_modulus
def Hubbles_Law(H0, redshift):
    luminosity_distance = c*redshift/H0
    return luminosity_distance
def non_Linear_Hubbles_Law(redshift, H0, q):
    luminosity_distance = c*redshift/H0*(1 + ((1 - q)/2)*redshift)
    return luminosity_distance
# def H02(luminosity_distance, redshift, H0):
#     calc_lum_dist = Hubbles_Law(redshift, H0)
#     differences = luminosity_distance-calc_lum_dist
#     the_fraction = differences**2/calc_lum_dist
#     X2 = np.sum(the_fraction)
#     return X2
# def X2(luminosity_distance, redshift, H0, q):
#     calc_lum_dist = non_Linear_Hubbles_Law(redshift, H0, q)
#     differences = luminosity_distance-calc_lum_dist
#     the_fraction = differences**2/calc_lum_dist
#     X2 = np.sum(the_fraction)
#     return X2
def integrate(z, Ohm_r, Ohm_m, Ohm_k, Ohm_A):
    integrate_this = 1/(Ohm_r*(1+z)**4+Ohm_m*(1+z)**3+Ohm_k*(1+z)**2+Ohm_A)
    return integrate_this
def FLRW(redshift,H0,Ohm_m,Ohm_k,Ohm_r):
    d_H = c/H0
    Ohm_A = 1 - Ohm_m - Ohm_k
    integral = np.array([quad(integrate, 0, z, args=(Ohm_r, Ohm_m, Ohm_k, Ohm_A,))[0] for z in redshift])
    d_C = d_H*integral
    d_M = d_C
    d_L = (1+redshift)*d_M
    return d_L
def Ohm_k_0(redshift,H0,Ohm_m,Ohm_r):
    return FLRW(redshift,H0,Ohm_m,0,Ohm_r)
#def best_q(H0, redshift, luminosity_distance):
#    for
f = open("SCPUnion2.1_mu_vs_z.txt")
lines=f.readlines()
col1 = np.empty(0)
col2 = np.empty(0)
col3 = np.empty(0)
col4 = np.empty(0)
col5 = np.empty(0)
for i in range(len(lines)):
    if i >= 5:
        values = lines[i].split("\t")
        col1 = np.append(col1,values[0])
        col2 = np.append(col2,float(values[1]))
        col3 = np.append(col3,float(values[2]))
        col4 = np.append(col4,float(values[3]))
        col5 = np.append(col5,float(values[4]))
#1: parse the file into relevant arrays, convert to luminosity distance
dists = col3
reds = col2
uncs = col4
lum_dists, lum_uncs = distance_modulus_to_luminosity_distance(dists,uncs)

#2: plot luminosity distance against redshift, distance modulus against redshift
mpl.subplot(121)
mpl.scatter(reds, lum_dists, s = 1, color = 'red')
mpl.xlabel('redshift')
mpl.ylabel('luminosity distance')
mpl.subplot(122)
mpl.scatter(reds, dists, s = 1, color = 'red')
mpl.ylabel('distance modulus')
mpl.xlabel('redshift')
mpl.show()
#mean, std:
#too big:
#62017.92595378953
#8208.293963027625
#0.3:
#68386.45936420794
#6239.004435703351
#0.5:
#68347.2301808295
#5754.287266619709

#3: Find the linear H0 for the first redshifts
H0s = np.empty(0)
for i in range(len(reds)):
    if reds[i]<0.05:
        H0s=np.append(H0s,c*reds[i]/lum_dists[i])
H0_lin = np.mean(H0s)
predicted_hubble = Hubbles_Law(H0_lin, reds)
difference = lum_dists-predicted_hubble
# X2 = H02(lum_dists, reds, H0_lin)
print("H0:", H0_lin)

#4: Plot with our linear regression with our H0
mpl.subplot(121)
mpl.scatter(reds, lum_dists, s = 1, color = 'red')
mpl.plot(reds, predicted_hubble, linewidth = 1, color = 'black')
mpl.xlabel('redshift')
mpl.ylabel('luminosity distance')
# mpl.show()
mpl.subplot(122)
mpl.scatter(predicted_hubble, difference/lum_uncs, s = 1, color = 'purple')
mpl.xlabel('predicted')
mpl.ylabel('normalized residuals')
mpl.show()

#5. Fit to our data with the non-linear Hubble's law with a deceleration parameter,
#   via the curve-fit sci-py method. Plot this fit and the normalized residuals
parameters, param_Covariance = curve_fit(non_Linear_Hubbles_Law, reds, lum_dists, sigma = lum_uncs)
H0_deceleration_scipy, q = parameters
#[6.01911023e+04 2.86107177e-01]
print("H0, q:",parameters)
points = np.arange(0.001,np.max(reds),.001)
calculated_with_deceleration = non_Linear_Hubbles_Law(points, H0_deceleration_scipy, q)

mpl.subplot(121)
mpl.scatter(reds, lum_dists, s = 1, color = 'red')
mpl.scatter(points, calculated_with_deceleration, s = 1, color = 'black')
mpl.xlabel('redshift')
mpl.ylabel('luminosity distance')
# mpl.show()
mpl.subplot(122)
# mpl.scatter(reds, dists, s = 1, color = 'red')
# mpl.scatter(points, luminosity_distance_to_distance_modulus(calculated_with_deceleration), s = 1, color = 'black')
# mpl.xlabel('redshift')
# mpl.ylabel('distance modulus')
# mpl.show()
#residuals
predicted_non_linear = non_Linear_Hubbles_Law(reds, H0_deceleration_scipy, q)
difference1 = lum_dists - predicted_non_linear
norm_resid = difference1/lum_uncs
mpl.scatter(predicted_non_linear, norm_resid, s = 1, color = 'purple')
mpl.xlabel('predicted')
mpl.ylabel('normalized residuals')
mpl.show()

#6. Fit to our data with the FLRW via the curve-fit sci-py method.
#   Plot this fit and the normalized residuals.
Ohm_k = 0
parameters2, param_Covariance2 = curve_fit(Ohm_k_0, reds, lum_dists, sigma = lum_uncs)
H0_FLRW_scipy, Ohm_m, Ohm_r = parameters2
calculated_with_FLRW = FLRW(points, H0_FLRW_scipy, Ohm_m, Ohm_k, Ohm_r)
mpl.subplot(121)
mpl.scatter(reds, lum_dists, s = 1, color = 'red')
mpl.scatter(points, calculated_with_FLRW, s = 1, color = 'black')
mpl.xlabel('redshift')
mpl.ylabel('luminosity distance')
# mpl.show()
mpl.subplot(122)
# mpl.scatter(reds, dists, s = 1, color = 'red')
# mpl.scatter(points, luminosity_distance_to_distance_modulus(calculated_with_FLRW), s = 1, color = 'black')
# mpl.xlabel('redshift')
# mpl.ylabel('distance modulus')
# mpl.show()
###residuals
predicted_FLRW = FLRW(reds, H0_FLRW_scipy, Ohm_m, Ohm_k, Ohm_r)
difference2 = lum_dists - predicted_FLRW
norm_resid2 = difference2/lum_uncs
mpl.scatter(predicted_FLRW, norm_resid2, s = 1, color = 'purple')#/(lum_dists*)param_Covariance2[0,1]
mpl.xlabel('predicted')
mpl.ylabel('normalized residuals')
mpl.show()
print ("H0, Ohm_m, Ohm_r:", parameters2)
#The statistical significance of Ohm_m and therefore Ohm_A.
print('Ohm_A: ', 1-Ohm_m, 'Covariance Ohm_m: ', sqrt(param_Covariance2[1,1]))
# Without using radiation: ('Ohm_A: ', 0.8854955179779508, 'Covariance Ohm_m: ', 0.0001755520748644085)
# With radiation: ('Ohm_A: ', 0.7113206558831211, 'Covariance Ohm_m: ', 0.12957568504213926)
# radiation changes the significance, so it is significant

#Again, but the universe is curved.
del Ohm_k
parameters21, param_Covariance21 = curve_fit(FLRW, reds, lum_dists, sigma = lum_uncs)
H0_FLRW_scipy, Ohm_m, Ohm_k, Ohm_r = parameters21
calculated_with_FLRW = FLRW(points, H0_FLRW_scipy, Ohm_m, Ohm_k, Ohm_r)
mpl.subplot(121)
mpl.scatter(reds, lum_dists, s = 1, color = 'red')
mpl.scatter(points, calculated_with_FLRW, s = 1, color = 'black')
mpl.xlabel('redshift')
mpl.ylabel('luminosity distance')
mpl.subplot(122)
predicted_FLRW = FLRW(reds, H0_FLRW_scipy, Ohm_m, Ohm_k, Ohm_r)
difference21 = lum_dists - predicted_FLRW
norm_resid21 = difference21/lum_uncs
mpl.scatter(predicted_FLRW, norm_resid21, s = 1, color = 'purple')
mpl.xlabel('predicted')
mpl.ylabel('normalized residuals')
mpl.show()
print ("H0, Ohm_m, Ohm_r:", parameters21)
#The statistical significance of Ohm_m and therefore Ohm_A.
print('Ohm_A: ', 1-Ohm_m-Ohm_k, 'Covariance Ohm_m: ', sqrt(param_Covariance21[1,1])+sqrt(param_Covariance21[2,2]))
f.close()
