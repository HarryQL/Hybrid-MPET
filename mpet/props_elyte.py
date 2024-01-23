r"""
This module provides functions defining properties of the ion-conducting
phase -- the electrolyte Manage functions for the parameters involved in
Stefan-Maxwell based concentrated electrolyte transport theory for
binary electrolytes.

Each electrolyte set must output functions for the following as a
function of c (electrolyte concentration, M)
 - Dchem [m^2/s] = the prefactor for grad(c) in species conservation
 - sigma [S/m] = the conductivity
 - (1 + dln(f_\pm)/dln(c)) = the "thermodynamic factor"
 - t_+^0 = the transference number of the cations
T in these equations is nondimensionalized wrt 298K
"""

import numpy as np
from mpet.config import constants


def LiClO4_PC():
    """ Set of parameters from Fuller, Doyle, Newman 1994, with
    conductivity directly from dualfoil5.2.f
    """
    def tp0(c, T):
        return 0.2

    def D(c, T):
        return 2.58e-10  # m^2/s

    def therm_fac(c, T):
        # therm_fac adjusted to account for missing a factor of 2 in Eq A-2
        return 0.5

    def sigma(cin, T):
        c = cin * 1000  # mol/m^3
        p_max = 0.542
        p_u = 0.6616
        a = 0.855
        b = -0.08
        rho = 1.2041e3
        out = 0.0001 + c**a * (
            p_max*(1./(rho*p_u))**a
            * np.exp(b*(c/rho - p_u)**2
                     - (a/p_u)*(c/rho - p_u)))  # S/m
        return out
    Dref = D(constants.c_ref/1000, 1)

    def D_ndim(c, T):
        return D(c, T) / Dref

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*(constants.c_ref)))
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref


def Doyle96_EC_DMC_2_1():
    """ Set of parameters from Doyle, Newman, et al. 1996.
    """
    def tp0(c, T):
        return 0.363

    def D(c, T):
        return 7.5e-11  # m^2/s

    def therm_fac(c, T):
        return 1.

    def sigma(c, T):
        r1 = 4.1253e-4
        r2 = 5.007e-3
        r3 = 4.7212e-3
        r4 = 1.5094e-3
        r5 = 1.6018e-4
        k0 = r1 + r2*c - r3*c**2 + r4*c**3 - r5*c**4  # S/cm
        return(100*k0)

    Dref = D(constants.c_ref, 1)

    def D_ndim(c, T):
        return D(c, T) / Dref

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*constants.c_ref))
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref


def Doyle96_EC_DMC_1_2():
    """ Set of parameters from Doyle, Newman, et al. 1996.
    """
    def tp0(c, T):
        return 0.363

    def D(c, T):
        return 7.5e-11  # m^2/s

    def therm_fac(c, T):
        return 1.

    def sigma(c, T):
        r1 = 1.0793e-4
        r2 = 6.7461e-3
        r3 = 5.2245e-3
        r4 = 1.3605e-3
        r5 = 1.1724e-4
        k0 = r1 + r2*c - r3*c**2 + r4*c**3 - r5*c**4  # S/cm
        return(100*k0)

    Dref = D(constants.c_ref, 1)

    def D_ndim(c, T):
        return D(c, T) / Dref

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*constants.c_ref))
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref


def valoen_reimers():
    """ Set of parameters from Valoen and Reimers 2005 """
    def tp0(c, T):
        return 0.38

    def D(c, T):
        return 2*(
            10**(-4) * 10**(-4.43 - 54/((T*constants.T_ref) - (229 + 5*c)) - 0.22*c))  # m^2/s

    def therm_fac(c, T):
        tmp = 0.601 - 0.24*c**(0.5) + 0.982*(1 - 0.0052*((T*constants.T_ref) - 294))*c**(1.5)
        return tmp/(1-tp0(c, T))

    def sigma(c, T):
        (k00, k01, k02,
         k10, k11, k12,
         k20, k21) = (
             -10.5, 0.0740, -6.96e-5,
             0.668, -0.0178, 2.80e-5,
             0.494, -8.86e-4)
        out = c * (k00 + k01*(T*constants.T_ref) + k02*(T*constants.T_ref)**2
                   + k10*c + k11*c*(T*constants.T_ref) + k12*c*(T*constants.T_ref)**2
                   + k20*c**2 + k21*c**2*(T*constants.T_ref))**2  # mS/cm
        out *= 0.1  # S/m
        return out
    Dref = D(constants.c_ref/1000, 1)

    def D_ndim(c, T):
        return D(c, T) / Dref

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*(constants.c_ref)))
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref


def valoen_bernardi():
    """ Set of parameters from Bernardi and Go 2011, indirectly from
    Valoen and Reimers 2005. The only change from Valoen and Reimers
    is the conductivity.
    """
    D_ndim, Ign, therm_fac, tp0, Dref = valoen_reimers()

    def sigma(c, T):
        (k00, k01, k02,
         k10, k11, k12,
         k20, k21) = (
            -8.2488, 0.053248, -0.000029871,
            0.26235, -0.0093063, 0.000008069,
            0.22002, -0.0001765)
        out = c * (k00 + k01*(T*constants.T_ref) + k02*(T*constants.T_ref)**2
                   + k10*c + k11*c*(T*constants.T_ref) + k12*c*(T*constants.T_ref)**2
                   + k20*c**2 + k21*c**2*(T*constants.T_ref))**2  # mS/cm
        out *= 0.1   # S/m
        return out

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*(constants.c_ref)))
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref

def AEM():
    """ Set of parameters from Medtronic """
    def tp0(c, T):
        a = 0.423036316174248
        b = 0.0330735864437778
        d = -0.00720085690182939
        e = 0.00135507513887536
        f = -0.00648048504478617
        r = 0.238118744821976
        return a + b*c + d*c**2 + e*c**3 + f*(T*constants.T_ref-constants.T_ref)**r

    def D(c, T):
        a = 5.78146591951108
        b = 832.5127084508
        d = 62.6262620610715
        e = -4.44756287439382
        ee = -0.491427045954456
        f = 5.63853120634445
        ff =0.107101830129255
        return 3*10**(-a
                    -b/(T*constants.T_ref - (d +e*c +f*c**2))
                    +ee*c +ff*c**2)  # m^2/s

    def therm_fac(c, T):
        a = 0.696070566996749
        b = -5.80198978215293
        d = 24.820954868413
        e = -42.6835219535502
        f = 32.4889860113285
        g = -9.07900445005268
        h = 0.180490199371271
        return a + b*c + d*c**2 + e*c**3 + f*c**4 + g*c**5 + h*(1-0.0052*(T*constants.T_ref-constants.T_ref))*c**8

    def sigma(c, T):
        a = -5.61645264355501
        b = 1.95619966573733
        d = 0.0301557801774189
        e = -0.00922766463833395
        f = -0.111890908921336
        h = 0.000735394190843922
        k = -0.350822834815855
        s = 0.0000818632446102486
        return 1*(k-c+ c**0.634*(a + b*c + d*(T*constants.T_ref)  + e *c*(T*constants.T_ref)  + f*c**2 + h*c**2*(T*constants.T_ref) ) - s *5**(c))

    Dref = D(constants.c_ref/1000, 1)

    def D_ndim(c, T):
        return D(c, T) / Dref

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*(constants.c_ref)))
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref

def LIONSIMBA_nonisothermal():
    """ Set of parameters from LIONSIMBA validation. Torchio et al, 2016.
    """

    def tp0(c, T):
        return 0.364

    def sigma(c, T):
        c_dim = c*1000  # dimensionalized c
        T_dim = T*constants.T_ref
        r1 = -10.5
        r2 = 0.668e-3
        r3 = 0.494e-6
        r4 = 0.074
        r5 = -1.78e-5
        r6 = -8.86e-10
        r7 = -6.96e-5
        r8 = 2.8e-8
        sig_out = 1e-4 * c_dim * (r1 + r2*c_dim + r3*c_dim**2 + T_dim
                                  * (r4 + r5*c_dim + r6*c_dim**2)
                                  + T_dim**2 * (r7 + r8*c_dim))**2
        return sig_out  # m^2/s

    def D(c, T):
        c_dim = c*1000
        T_dim = T*constants.T_ref
        r1 = 4.43
        r2 = 54
        r3 = 229
        r4 = 5e-3
        r5 = 0.22e-3
        D_out = 1e-4 * 10**(-r1-r2/(T_dim-r3-r4*c_dim)-r5*c_dim)
        return D_out

    def therm_fac(c, T):
        return 1.

    Dref = D(constants.c_ref/1000, 1)

    def D_ndim(c, T):
        return D(c, T) / Dref

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*(constants.c_ref)))
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref


def LIONSIMBA_isothermal():
    """ Set of parameters from LIONSIMBA validation. Torchio et al, 2016.
    """
    def tp0(c, T):
        return 0.364

    def sigma(c, T):
        ce = c*1000  # dimensionalized c
        return (4.1253e-2 + 5.007e-4*ce - 4.7212e-7*ce**2
                + 1.5094e-10*ce**3 - 1.6018*1e-14*ce**4)  # S/m

    def D(c, T):
        return 7.5e-10  # m^2/s
    # isothermal at 298 K

    def therm_fac(c, T):
        return 1.

    Dref = D(constants.c_ref/1000, 1)

    def D_ndim(c, T):
        return D(c, T) / Dref

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*(constants.c_ref)))
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref
