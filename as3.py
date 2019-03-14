import numpy as np
import matplotlib.pyplot as mpl
from math import cos, sin
import sys
def exEuler(h, ts, xi, vi):
    xs = np.array([xi])
    vs = np.array([vi])
    for n in range(len(ts)-1):
        xi = xs[-1]
        vi = vs[-1]
        xi1 = xi + h*vi
        xs = np.append(xs, xi1)
        vi1 = vi - h*xi
        vs = np.append(vs, vi1)
    return xs, vs
def real(ts):
    realXs = np.cos(ts)
    realVs = -np.sin(ts)
    return realXs, realVs
def error(xs, vs, realXs, realVs):
    XError = abs(realXs-xs)
    VError = abs(realVs-vs)
    return XError, VError
def errorChange(hh, ns, xi, vi):
    k = np.arange(10)
    tt = ns*hh/2**k[-1]
    ts = ns*hh
    hs = hh/2**k
    realXs, realVs = real(ts)
    XErrorMaxs = np.empty(0)
    VErrorMaxs = np.empty(0)
    for h in hs:
        ts = np.arange(0,tt[-1]+h, h)
        realXs, realVs = real(ts)
        xs, vs = exEuler(h, ts, xi, vi)
        XError, VError = error(xs, vs, realXs, realVs)
        XErrorMaxs = np.append(XErrorMaxs, max(XError))
        VErrorMaxs = np.append(VErrorMaxs, max(VError))
    return hs, XErrorMaxs, VErrorMaxs
def imEuler(h, ts, xi, vi):
    xs = np.array([xi])
    vs = np.array([vi])
    for n in range(len(ts)-1):
        xi = xs[-1]
        vi = vs[-1]
        xi1 = (xi + h*vi)/(1+h**2)
        xs = np.append(xs, xi1)
        vi1 = (vi - h*xi)/(1+h**2)
        vs = np.append(vs, vi1)
    return xs, vs
def symEuler(h, ts, xi, vi):
    xs = np.array([xi])
    vs = np.array([vi])
    for n in range(len(ts)-1):
        xi = xs[-1]
        vi = vs[-1]
        xi1 = xi + h*vi
        xs = np.append(xs, xi1)
        vi1 = vi - h*xi1
        vs = np.append(vs, vi1)
    return xs, vs
def Energies(xs, vs):
    Es = xs**2 + vs**2
    return Es

h = 0.01
n = 100000
xi = 1
vi = 0

ns = np.arange(n+1)
ts = ns*h

# exs, evs = exEuler(h, ts, xi, vi)
# realXs, realVs = real(ts)
# ixs, ivs = imEuler(h, ts, xi, vi)
# sxs, svs = symEuler(h, ts, xi, vi)
#
# eXError, eVError = error(exs, evs, realXs, realVs)
# iXError, iVError = error(ixs, ivs, realXs, realVs)
#
# eEs = Energies(exs, evs)
# iEs = Energies(ixs, ivs)
# sEs = Energies(sxs, svs)

plot = int(sys.argv[1])
#positions and velocities for explicit
if plot == 1:
    exs, evs = exEuler(h, ts, xi, vi)
    mpl.subplot(121)
    mpl.plot(ts, exs)
    mpl.ylabel('xs')
    mpl.xlabel('t')
    mpl.subplot(122)
    mpl.plot(ts, evs)
    mpl.ylabel('vs')
    mpl.xlabel('t')
    mpl.savefig('plot1.pdf')
#error for Xs and Vs in explicit method
if plot == 2:
    exs, evs = exEuler(h, ts, xi, vi)
    realXs, realVs = real(ts)
    eXError, eVError = error(exs, evs, realXs, realVs)
    mpl.subplot(121)
    mpl.plot(ts, eXError)
    mpl.ylabel('XErrors')
    mpl.xlabel('t')
    mpl.subplot(122)
    mpl.plot(ts, eVError)
    mpl.ylabel('VErrors')
    mpl.xlabel('t')
    mpl.savefig('plot2.pdf')
#maximum error and h proportional for explicit
if plot == 3:
    hs, XErrorMaxs, VErrorMaxs = errorChange(h, ns, xi, vi)
    mpl.subplot(121)
    mpl.loglog(hs, XErrorMaxs)
    mpl.ylabel('XErrorsMaximums')
    mpl.xlabel('hs')
    mpl.subplot(122)
    mpl.loglog(hs, VErrorMaxs)
    mpl.ylabel('VErrorMaximums')
    mpl.xlabel('hs')
    mpl.savefig('plot3.pdf')
#energy evolution for explicit
if plot == 4:
    exs, evs = exEuler(h, ts, xi, vi)
    eEs = Energies(exs, evs)
    mpl.plot(ts, eEs)
    mpl.ylabel('Energies')
    mpl.xlabel('t')
    mpl.savefig('plot4.pdf')
#positions and velocities for implicit
if plot == 5:
    ixs, ivs = imEuler(h, ts, xi, vi)
    mpl.subplot(121)
    mpl.plot(ts, ixs)
    mpl.ylabel('xs')
    mpl.xlabel('t')
    mpl.subplot(122)
    mpl.plot(ts, ivs)
    mpl.ylabel('vs')
    mpl.xlabel('t')
    mpl.savefig('plot5.pdf')
#Energies, errors for implicit
if plot == 6:
    exs, evs = exEuler(h, ts, xi, vi)
    ixs, ivs = imEuler(h, ts, xi, vi)
    realXs, realVs = real(ts)
    eXError, eVError = error(exs, evs, realXs, realVs)
    iXError, iVError = error(ixs, ivs, realXs, realVs)
    eEs = Energies(exs, evs)
    iEs = Energies(ixs, ivs)
    mpl.subplot(321)
    mpl.plot(ts, iXError)
    mpl.ylabel('Implicit XErrors')
    mpl.xlabel('t')
    mpl.subplot(322)
    mpl.plot(ts, eXError)
    mpl.ylabel('Explicit XErrors')
    mpl.xlabel('t')
    mpl.subplot(323)
    mpl.plot(ts, iVError)
    mpl.ylabel('Implicit VErrors')
    mpl.xlabel('t')
    mpl.subplot(324)
    mpl.plot(ts, eVError)
    mpl.ylabel('Explicit VErrors')
    mpl.xlabel('t')
    mpl.subplot(325)
    mpl.plot(ts, iEs)
    mpl.ylabel('Implicit Energies')
    mpl.xlabel('t')
    mpl.subplot(326)
    mpl.plot(ts, eEs)
    mpl.ylabel('Explicit Energies')
    mpl.xlabel('t')
    mpl.savefig('plot6.pdf')
#phase space geometries for explicit and implicit
if plot == 7:
    exs, evs = exEuler(h, ts, xi, vi)
    ixs, ivs = imEuler(h, ts, xi, vi)
    realXs, realVs = real(ts)
    mpl.subplot(311)
    mpl.plot(exs, evs)
    mpl.ylabel('Explicit Vs')
    mpl.xlabel('Explicit Xs')
    mpl.subplot(312)
    mpl.plot(ixs, ivs)
    mpl.ylabel('Implicit Vs')
    mpl.xlabel('Implicit Xs')
    mpl.subplot(313)
    mpl.plot(realXs, realVs)
    mpl.ylabel('Analytic Vs')
    mpl.xlabel('Analytic Xs')
    mpl.savefig('plot7.pdf')
#positions and velocities with the symplectic method
if plot == 8:
    sxs, svs = symEuler(h, ts, xi, vi)
    mpl.subplot(121)
    mpl.plot(ts, sxs)
    mpl.ylabel('xs')
    mpl.xlabel('t')
    mpl.subplot(122)
    mpl.plot(ts, svs)
    mpl.ylabel('vs')
    mpl.xlabel('t')
    mpl.savefig('plot8.pdf')
#phase space geometries with symplectic
if plot == 9:
    exs, evs = exEuler(h, ts, xi, vi)
    ixs, ivs = imEuler(h, ts, xi, vi)
    sxs, svs = symEuler(h, ts, xi, vi)
    realXs, realVs = real(ts)
    mpl.subplot(221)
    mpl.plot(sxs, svs)
    mpl.ylabel('Symplectic Vs')
    mpl.xlabel('Symplectic Xs')
    mpl.subplot(223)
    mpl.plot(exs, evs)
    mpl.ylabel('Explicit Vs')
    mpl.xlabel('Explicit Xs')
    mpl.subplot(224)
    mpl.plot(ixs, ivs)
    mpl.ylabel('Implicit Vs')
    mpl.xlabel('Implicit Xs')
    mpl.subplot(222)
    mpl.plot(realXs, realVs)
    mpl.ylabel('Analytic Vs')
    mpl.xlabel('Analytic Xs')
    mpl.savefig('plot9.pdf')
#energy evolution for the symplectic method
if plot == 10:
    sxs, svs = symEuler(h, ts, xi, vi)
    sEs = Energies(sxs, svs)
    mpl.plot(ts, sEs)
    mpl.ylabel('Energies')
    mpl.xlabel('t')
    mpl.savefig('plot10.pdf')
