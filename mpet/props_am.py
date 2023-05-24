"""This module handles properties associated with the active materials."""
import numpy as np

import mpet.geometry as geo
from mpet.config import constants


class Dfuncs():
    """This class returns the filling-fraction dependent variation of
    the transport coefficient, D(y), such that
    Flux = -D_ref*D(y)*grad(y) for solid solution transport or
    Flux = -D_ref*D(y)*grad(mu) for thermo-based transport
    where y here is the filling fraction, D_ref has dimensions of
    length^2/time, D(y) is dimensionless, and mu, the chemical
    potential, has been scaled to the thermal energy, k*Tref. For more
    details on the two forms, see muRfuncs which defines both solid
    solution (_ss) and materials based on simpler thermodynamic models.
    """

    def __init__(self, Dfunc):
        Dopts = {}
        Dopts['lattice'] = self.lattice
        Dopts['constant'] = self.constant
        self.Dfunc = Dopts[Dfunc]

    def constant(self, y):
        return 1.

    def lattice(self, y):
        return y*(1-y)


class muRfuncs():
    """ This class defines functions which describe the chemical
    potential of active materials.
    Each function describes a particular material. The
    chemical potential is directly related to the open circuit voltage
    (OCV) for solid solution materials.
    In each function, muR_ref is the offset (non-dimensional) chemical
    potential by which the returned value will be shifted (for
    numerical convenience).
    Each material function returns:
        muR -- chemical potential
        actR -- activity (if applicable, else None)
    """
    def __init__(self, config, trode, ind=None):
        """config is the full dictionary of
        parameters for the electrode particles, as made for the
        simulations. trode is the selected electrode.
        ind is optinally the selected particle, provided as (vInd, pInd)
        """
        self.config = config
        self.trode = trode
        self.ind = ind
        self.T = config['T']  # nondimensional
        # eokT and kToe are the reference values for scalings
        self.eokT = constants.e / (constants.k * constants.T_ref)
        self.kToe = 1. / self.eokT

        # Convert "muRfunc" to a callable function
        self.muRfunc = getattr(self, self.get_trode_param("muRfunc"))
        self.muRfunc2 = getattr(self, self.get_trode_param("muRfunc2"))

    def get_trode_param(self, item):
        """
        Shorthand to retrieve electrode-specific value
        """
        value = self.config[self.trode, item]
        # check if it is a particle-specific parameter
        if self.ind is not None and item in self.config.params_per_particle:
            value = value[self.ind]
        return value

    def get_muR_from_OCV(self, OCV, muR_ref):
        return -self.eokT*OCV + muR_ref

    ######
    # Solid solution functions
    # These are all obtained directly from fitting an OCV of the
    # material with a standard counter electrode.
    # They can all only return values at 298 K
    ######

    def LiMn2O4_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """ Doyle, Newman, 1996 """
        # OCV in V vs Li/Li+
        OCV = (4.19829 + 0.0565661*np.tanh(-14.5546*y + 8.60942)
               - 0.0275479*(1/((0.998432 - y)**(0.492465)) - 1.90111)
               - 0.157123*np.exp(-0.04738*y**8)
               + 0.810239*np.exp(-40*(y - 0.133875)))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def LiMn2O4_ss2(self, y, ybar, muR_ref, ISfuncs=None):
        """ Fuller, Doyle, Newman, 1994 """
        # OCV in V vs Li/Li+
        OCV = (4.06279 + 0.0677504*np.tanh(-21.8502*y + 12.8268)
               - 0.105734*(1/((1.00167 - y)**(0.379571)) - 1.575994)
               - 0.045*np.exp(-71.69*y**8)
               + 0.01*np.exp(-200*(y - 0.19)))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def LiC6_coke_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """ Doyle, Newman, 1996 """
        OCV = (-0.16 + 1.32*np.exp(-3.0*y) + 10.*np.exp(-2000.*y))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def LiC6_coke_ss2(self, y, ybar, muR_ref, ISfuncs=None):
        """ Fuller, Doyle, Newman, 1994 """
        c1 = -0.132056
        c2 = 1.40854
        c3 = -3.52312
        OCV = c1 + c2*np.exp(c3*y)
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def LiC6_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """ Safari, Delacourt 2011 """
        OCV = (0.6379 + 0.5416*np.exp(-305.5309*y)
               + 0.044*np.tanh(-(y - 0.1958)/0.1088)
               - 0.1978*np.tanh((y - 1.0571)/0.0854)
               - 0.6875*np.tanh((y + 0.0117)/0.0529)
               - 0.0175*np.tanh((y - 0.5692)/0.0875))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def LiC6_ss2(self, y, ybar, muR_ref, ISfuncs=None):
        """ Bernardi and Go 2011 """
        p1, p2, p3, p4 = (0.085, 0.120, 0.210, 3.5)
        sfac = 0.3
        OCV = (p1*step_down(y, 1., sfac*0.02)
               + (p2 - p1)*step_down(y, 0.5, 0.005)
               + (p3 - p2)*step_down(y, 0.1944, sfac*0.03571)
               + (p4 - p3)*step_down(y, 0., sfac*0.08333))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def LiC6_2step_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """
        Fit function to the OCV predicted by the phase separating
        2-variable graphite model (LiC6 function in this class).
        """
        Vstd = 0.1196
        Vstep = 0.0351733976
        edgeLen = 0.024
        lEdge = edgeLen
        rEdge = 1 - edgeLen
        width = 1e-4
        vshift = 1e-2
        lSide = -((np.log(y/(1-y)) - np.log(lEdge/(1-lEdge)) - vshift)*step_down(y, lEdge, width))
        rSide = -((np.log(y/(1-y)) - np.log(rEdge/(1-rEdge)) + vshift)*step_up(y, rEdge, width))
        OCV = (
            Vstd
            + Vstep*(step_down(y, 0.5, 0.013) - 1)
            + self.kToe*(lSide + rSide)
            )
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def Li_ss(self, y, ybar, muR_ref, ISfuncs=None):
        muR = 0.*y + muR_ref
        actR = 0.*y + 1.
        return muR, actR

    def LiSVO1_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """
        This function was obtained from Qiaohao Liang's fit of Gomadam et al. and Crespi et al.
        discharge data (not 298K but at 37 celcius) of silver vanadium oxide (SVO) as cathode.
        """
        a = 0.170853327061503
        c = 0.10575732541213
        cc = 0.0325999695981445
        d = -0.78083686412957
        dd = -0.240750600055772
        e = 2.45379781220411
        ee = 0.772378929022569
        f = -3.97315419342551
        ff = -1.42183809509377
        g = 2.18426393758425
        gg = 1.64989436658455
        h = 2.86023862411551
        hh = -1.56778412967418
        l = -2.38548221102615
        ll = 2.31136343789118
        m = -3.60162723777234
        mm = -1.79871093593143
        n = 1.31355252735059
        nn = -2.31932170090989
        o = 4.44535068501561
        oo = 4.33115673161338
        p = -2.61329695121938
        pp = -1.74507099094851

        OCV = +a*np.exp(-300*y) \
        +(c+  d*y+  e*(y**2)+ f*(y**3)+  g*(y**4)+  h*(y**5)+  l*(y**6)+  m*(y**7)+  n*(y**8)+  o*(y**9)+  p*(y**10)) \
        /(cc+dd*y+ee*(y**2)+ ff*(y**3)+ gg*(y**4)+ hh*(y**5)+ ll*(y**6)+ mm*(y**7)+ nn*(y**8)+ oo*(y**9)+ pp*(y**10))

        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        # EoKT = constants.e / (constants.k * 310.15)
        # muR = -EoKT*OCV + muR_ref
        # actR = None
        return muR, actR

    def LiVO1_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """
        Vanadium reduction OCV
        """

        a = 0.777202437736476
        b = 0.474744924663168
        d = 0.689210000197934
        dd = 0.231386236456565
        e = 0.900086230881169
        ee = 1.38659662497403
        f = 2.19792050603498
        ff = -0.570214017959082
        g = -0.247566668286967
        gg = 0.610299873935577
        h = -0.828002605941331
        hh = 0.882567947451821
        k = 0.424473584337631
        kk = 0.0423284925890506
        l = 0.828247009245531
        ll = -0.547376068752914
        m = 0.553745857338048
        mm = -0.94353674601898
        n = 0.916422586764535
        nn = 0.118145979574833
        o = -0.440939258125685
        oo = 1.01429721903779
        p = -0.468561937977976
        pp = 1.43586896994767

        OCV = +a +b*np.exp(-300*y)  \
              +(d + e*y+ f*y**2+g*y**3+h*y**4+k*y**5+l*y**6+m*y**7+n*y**8+o*y**9+p*y**10) \
              /(dd + ee*y+ ff*y**2+gg*y**3+hh*y**4+kk*y**5+ll*y**6+mm*y**7+nn*y**8+oo*y**9+pp*y**10)


        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        # EoKT = constants.e / (constants.k * 310.15)
        # muR = -EoKT*OCV + muR_ref
        # actR = None
        return muR, actR

    def LiVO2_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """
        Vanadium reduction OCV from Gomadam
        """

        OCV = +0.823078467*np.exp(-20*(4*y)) \
            +(3.176921533+5.802419895*((4*y)**2)+0.191983977*((4*y)**4)-0.16084978*((4*y)**6)+0.00900142*((4*y)**8)) \
            /(1+2.462745769*((4*y)**2)-0.024605961*((4*y)**4)-0.04188333*((4*y)**6)+0.001617672*((4*y)**8)+0.000062746*((4*y)**10))



        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        # EoKT = constants.e / (constants.k * 310.15)
        # muR = -EoKT*OCV + muR_ref
        # actR = None
        return muR, actR

    def LiCF_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """
        CFx (x=1) discharge OCV from Gomadam
        """
        b = -0.02
        c = 6.64598657006137e-06
        cc =2.12258603268465e-06
        d = -0.00553346316001277
        dd =-0.00182420160990677
        e = 2.54146296708079
        ee =0.8598721909619
        f = 5.81715187206634
        ff =1.63460834212519
        g =-3.50302948659574
        gg =0.000616092918954162
        h = -2.03863850063247
        hh =-3.81158221764385
        l = -4.9616107272193
        ll =3.66430362872806
        m = -2.15485220570863
        mm =-3.93194442432665
        n = -0.708878434759009
        nn =-2.46912720449879
        o = -0.214783881560837
        oo =3.3241731587368
        p = 5.17109932227109
        pp = 0.80576977271189

        OCV = b*np.log(y/(1-y)) \
                +(c+  d*y+  e*(y**2)+ f*(y**3)+  g*(y**4)+  h*(y**5)+  l*(y**6)+  m*(y**7)+  n*(y**8)+  o*(y**9)+  p*(y**10)) \
                /(cc+dd*y+ee*(y**2)+ ff*(y**3)+ gg*(y**4)+ hh*(y**5)+ ll*(y**6)+ mm*(y**7)+ nn*(y**8)+ oo*(y**9)+ pp*(y**10))



        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        # EoKT = constants.e / (constants.k * 310.15)
        # muR = -EoKT*OCV + muR_ref
        # actR = None
        return muR, actR

    def LiCF2_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """
        CFx (x=1) discharge OCV from Gomadam
        """
        b = -0.02
        c = -4.14639585277983e-07
        cc = -1.68721890317438e-06
        d = 0.0600539082629088
        dd = 0.0207303415063955
        e = 2.18701626123853
        ee = 0.717550319820223
        f = -1.37732821051454
        ff = -0.488765757256375
        g = -0.478977238048069
        gg = -0.164137424184735
        h = -1.31167352341344
        hh = -0.026312377902526
        l = -0.606919331393176
        ll = -0.692372665120114
        m = -0.0519265218115834
        mm = -0.315412340050245
        n = -0.102424920148515
        nn = 0.421979431767184
        o = 0.561841248734846
        oo = 0.51602268062699
        p = 1.27473684931705
        pp = 0.0524058323859283

        OCV = b*np.log(y/(1-y)) \
                +(c+  d*y+  e*(y**2)+ f*(y**3)+  g*(y**4)+  h*(y**5)+  l*(y**6)+  m*(y**7)+  n*(y**8)+  o*(y**9)+  p*(y**10)) \
                /(cc+dd*y+ee*(y**2)+ ff*(y**3)+ gg*(y**4)+ hh*(y**5)+ ll*(y**6)+ mm*(y**7)+ nn*(y**8)+ oo*(y**9)+ pp*(y**10))

        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        # EoKT = constants.e / (constants.k * 310.15)
        # muR = -EoKT*OCV + muR_ref
        # actR = None
        return muR, actR

    def LiCF3_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """
        CFx (x=1) discharge OCV from Gomadam
        """
        bra = -0.5508820669 + (5.566158288-1476.163146*y + 240126.5532*(y**2)-110008.0909*(y**3)-168834.5645*(y**4)) \
                          /(1-404.0847731*y + 75655.38173*(y**2)-32080.42547*(y**3)-55461.3021*(y**4))

        OCV = 0.02+ 0.08*(0.3-y)+bra - (1/0.57)*((310.15*1.38e-23)/(1.602e-19))*np.log((3*1.908/(10e-6*2.75e6)) * 5e-4 * y * (1-y)**(1/3))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        # EoKT = constants.e / (constants.k * 310.15)
        # muR = -EoKT*OCV + muR_ref
        # actR = None
        return muR, actR


    def NCA_ss1(self, y, ybar, muR_ref, ISfuncs=None):
        """
        This function was obtained from Dan Cogswell's fit of Samsung
        data.
        """
        OCV = (3.86 + 1.67*y - 9.52*y**2 + 15.04*y**3 - 7.95*y**4
               - 0.06*np.log(y/(1-y)))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def NCA_ss2(self, y, ybar, muR_ref, ISfuncs=None):
        """
        Li_q Ni(0.8)Co(0.15)Al(0.05)O2
        as a function of y. Here, y actually represents a practical
        utilization of 70% of the material, so the material is "empty"
        (y=0) when q=0.3 and full (y=1) when q=1.
        This function was obtained from a fit by Raymond B. Smith
        of Samsung data of a LiC6-NCA cell discharged at C/100.
        """
        OCV = (-self.kToe*np.log(y/(1-y))
               + 4.12178 - 0.2338*y - 1.24566*y**2 + 1.16769*y**3
               - 0.20745*y**4)
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def testIS_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """Ideal solution material for testing."""
        OCV = -self.kToe*np.log(y/(1-y))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def testRS_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """
        Regular solution material which phase separates at binodal points,
        for modeling as a solid solution. For testing.
        """
        # Based Omg = 3*k*T_ref
        yL = 0.07072018
        yR = 0.92927982
        OCV_rs = -self.kToe*self.reg_sln(y, self.get_trode_param("Omega_a"), ISfuncs)
        width = 0.005
        OCV = OCV_rs*step_down(y, yL, width) + OCV_rs*step_up(y, yR, width) + 2
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    ######
    # Functions based on thermodynamic models
    ######

    def ideal_sln(self, y, ISfuncs=None):
        """ Helper function: Should not be called directly from
        simulation. Call a specific material instead. """
        T = self.T
        if ISfuncs is not None:
            # Input must be a vector when using ISfuncs
            muR = T*np.array([ISfuncs[i]() for i in range(len(y))])
        else:
            muR = T*np.log(y/(1-y))
        return muR


    def reg_sln(self, y, Omga, ISfuncs=None):
        """ Helper function """
        muR_IS = self.ideal_sln(y, ISfuncs=ISfuncs)
        enthalpyTerm = Omga*(1-2*y)
        muR = muR_IS + enthalpyTerm
        return muR

    def graphite_2param_homog(self, y, Omga, Omgb, Omgc, EvdW, ISfuncs=None):
        """ Helper function """
        y1, y2 = y
        if ISfuncs is None:
            ISfuncs1, ISfuncs2 = None, None
        else:
            ISfuncs1, ISfuncs2 = ISfuncs
        muR1 = self.reg_sln(y1, Omga, ISfuncs1)
        muR2 = self.reg_sln(y2, Omga, ISfuncs2)
        muR1 += Omgb*y2 + Omgc*y2*(1-y2)*(1-2*y1)
        muR2 += Omgb*y1 + Omgc*y1*(1-y1)*(1-2*y2)
        muR1 += EvdW * (30 * y1**2 * (1-y1)**2)
        muR2 += EvdW * (30 * y2**2 * (1-y2)**2)
        return (muR1, muR2)

    def graphite_1param_homog(self, y, Omga, Omgb, ISfuncs=None):
        """ Helper function """
        width = 5e-2
        tailScl = 5e-2
        muLtail = -tailScl*1./(y**(0.85))
        muRtail = tailScl*1./((1-y)**(0.85))
        slpScl = 0.45
        muLlin = slpScl*Omga*4*(0.26-y)*step_down(y, 0.5, width)
        muRlin = (slpScl*Omga*4*(0.74-y) + Omgb)*step_up(y, 0.5, width)
        muR = muLtail + muRtail + muLlin + muRlin
        return muR

    def graphite_1param_homog_2(self, y, Omga, Omgb, ISfuncs=None):
        """ Helper function """
        width = 5e-2
        tailScl = 5e-2
        slpScl = 0.45
        muLtail = -tailScl*1./(y**(0.85))
        muRtail = tailScl*1./((1-y)**(0.85))
        muLlin = (slpScl*Omga*12*(0.40-y)
                  * step_down(y, 0.49, 0.9*width)*step_up(y, 0.35, width))
        muRlin = (slpScl*Omga*4*(0.74-y) + Omgb)*step_up(y, 0.5, width)
        muLMod = (0.
                  + 40*(-np.exp(-y/0.015))
                  + 0.75*(np.tanh((y-0.17)/0.02) - 1)
                  + 1.0*(np.tanh((y-0.22)/0.040) - 1)
                  )*step_down(y, 0.35, width)
        muR = muLMod + muLtail + muRtail + muLlin + muRlin
        return muR

    def graphite_1param_homog_3(self, y, Omga, Omgb, ISfuncs=None):
        """ Helper function with low hysteresis and soft tail """
        width = 5e-2
        tailScl = 5e-2
        muLtail = -tailScl*1./(y**(0.85))
        muRtail = tailScl*1./((1-y)**(0.85))
        muRtail = 1.0e1*step_up(y, 1.0, 0.045)
        muLlin = (0.15*Omga*12*(0.40-y**0.98)
                  * step_down(y, 0.49, 0.9*width)*step_up(y, 0.35, width))
        muRlin = (0.1*Omga*4*(0.74-y) + 0.90*Omgb)*step_up(y, 0.5, 0.4*width)
        muLMod = (0.
                  + 40*(-np.exp(-y/0.015))
                  + 0.75*(np.tanh((y-0.17)/0.02) - 1)
                  + 1.0*(np.tanh((y-0.22)/0.040) - 1)
                  )*step_down(y, 0.35, width)
        muR = 0.18 + muLMod + muLtail + muRtail + muLlin + muRlin
        return muR

    def non_homog_rect_fixed_csurf(self, y, ybar, B, kappa, ywet):
        """ Helper function """
        N = len(y)
        ytmp = np.empty(N+2, dtype=object)
        ytmp[1:-1] = y
        ytmp[0] = ywet
        ytmp[-1] = ywet
        dxs = 1./N
        curv = np.diff(ytmp, 2)/(dxs**2)
        muR_nh = -kappa*curv + B*(y - ybar)
        return muR_nh

    def non_homog_round_wetting(self, y, ybar, B, kappa, beta_s, shape, r_vec):
        """ Helper function """
        dr = r_vec[1] - r_vec[0]
        Rs = 1.
        curv = geo.calc_curv(y, dr, r_vec, Rs, beta_s, shape)
        muR_nh = B*(y - ybar) - kappa*curv
        return muR_nh

    def general_non_homog(self, y, ybar):
        """ Helper function """
        ptype = self.get_trode_param("type")
        mod1var, mod2var = False, False
        if isinstance(y, np.ndarray):
            mod1var = True
            N = len(y)
        elif (isinstance(y, tuple) and len(y) == 2
                and isinstance(y[0], np.ndarray)):
            mod2var = True
            N = len(y[0])
        else:
            raise Exception("Unknown input type")
        if ("homog" not in ptype) and (N > 1):
            shape = self.get_trode_param("shape")
            kappa = self.get_trode_param("kappa")
            B = self.get_trode_param("B")
            if shape == "C3":
                if mod1var:
                    cwet = self.get_trode_param("cwet")
                    muR_nh = self.non_homog_rect_fixed_csurf(
                        y, ybar, B, kappa, cwet)
                elif mod2var:
                    raise NotImplementedError("no 2param C3 model known")
            elif shape in ["cylinder", "sphere"]:
                beta_s = self.get_trode_param("beta_s")
                r_vec = geo.get_unit_solid_discr(shape, N)[0]
                if mod1var:
                    muR_nh = self.non_homog_round_wetting(
                        y, ybar, B, kappa, beta_s, shape, r_vec)
                elif mod2var:
                    muR1_nh = self.non_homog_round_wetting(
                        y[0], ybar[0], B, kappa, beta_s, shape, r_vec)
                    muR2_nh = self.non_homog_round_wetting(
                        y[1], ybar[1], B, kappa, beta_s, shape, r_vec)
                    muR_nh = (muR1_nh, muR2_nh)
        else:  # homogeneous particle
            if mod1var:
                muR_nh = 0*y
            elif mod2var:
                muR_nh = (0*y[0], 0*y[1])
        return muR_nh

    def LiAgO(self, y, ybar, muR_ref, ISfuncs=None):
        """ Harry QL """
        muRtheta = -self.eokT*3.24
        muRhomog = self.reg_sln(y, self.get_trode_param("Omega_a"), ISfuncs)
        muR = muRhomog
        actR = np.exp(muR/self.T)
        muR += muRtheta + muR_ref
        return muR, actR


    def LiAgO_5(self, y, ybar, muR_ref, ISfuncs=None):
        """ Harry QL """
        muRtheta = -self.eokT*3.24
        muR_IS = self.ideal_sln(y, ISfuncs=ISfuncs)

        Omga = self.get_trode_param("Omega_a")
        Omgb = self.get_trode_param("Omega_b")
        Omgc = self.get_trode_param("Omega_c")

        muR_enthalpy = Omga*(1-1.4*y) + Omgb*(y**2) + Omgc*(y**3)

        muR = muR_IS + muR_enthalpy
        actR = np.exp(muR/self.T)
        muR += muRtheta + muR_ref
        return muR, actR

    def LiAgO_6(self, y, ybar, muR_ref, ISfuncs=None):
        """ Harry QL """
        muRtheta = -self.eokT*3.24
        muR_IS = self.ideal_sln(y, ISfuncs=ISfuncs)

        Omga = self.get_trode_param("Omega_a")
        Omgb = self.get_trode_param("Omega_b")
        Omgc = self.get_trode_param("Omega_c")
        Omgd = self.get_trode_param("Omega_d")


        muR_enthalpy = Omga + Omgb*y + Omgc*(y**1.5) + Omgd*(y**2)

        muR = muR_IS + muR_enthalpy
        actR = np.exp(muR/self.T)
        muR += muRtheta + muR_ref
        return muR, actR

    def LiAgO_7(self, y, ybar, muR_ref, ISfuncs=None):
        """ Harry QL """
        muRtheta = -self.eokT*3.24
        muR_IS = self.ideal_sln(y, ISfuncs=ISfuncs)

        Omga = self.get_trode_param("Omega_a")
        Omgb = self.get_trode_param("Omega_b")
        Omgc = self.get_trode_param("Omega_c")
        Omgd = self.get_trode_param("Omega_d")

        muR_enthalpy = Omga + Omgb*y + Omgc*(y**2) + Omgd*(y**3)

        muR = muR_IS + muR_enthalpy
        actR = np.exp(muR/self.T)
        muR += muRtheta + muR_ref
        return muR, actR



    def LiFePO4(self, y, ybar, muR_ref, ISfuncs=None):
        """ Bai, Cogswell, Bazant 2011 """
        muRtheta = -self.eokT*3.422
        muRhomog = self.reg_sln(y, self.get_trode_param("Omega_a"), ISfuncs)
        muRnonHomog = self.general_non_homog(y, ybar)
        muR = muRhomog + muRnonHomog
        actR = np.exp(muR/self.T)
        muR += muRtheta + muR_ref
        return muR, actR

    def LiC6(self, y, ybar, muR_ref, ISfuncs=(None, None)):
        """ Ferguson and Bazant 2014 """
        muRtheta = -self.eokT*0.12
        muR1homog, muR2homog = self.graphite_2param_homog(
            y, self.get_trode_param("Omega_a"), self.get_trode_param("Omega_b"),
            self.get_trode_param("Omega_c"), self.get_trode_param("EvdW"), ISfuncs)
        muR1nonHomog, muR2nonHomog = self.general_non_homog(y, ybar)
        muR1 = muR1homog + muR1nonHomog
        muR2 = muR2homog + muR2nonHomog
        actR1 = np.exp(muR1/self.T)
        actR2 = np.exp(muR2/self.T)
        muR1 += muRtheta + muR_ref
        muR2 += muRtheta + muR_ref
        return (muR1, muR2), (actR1, actR2)

    def SVO_hybrid(self, y, ybar, muR_ref, ISfuncs=(None, None)):
        """Harry 2022 Medtronic"""

        y1, y2 = y
        # muRtheta2 = -self.eokT*2.59

        muR1homog, actR1 = self.LiAgO(y1, ybar, muR_ref, ISfuncs=None)

        muR2, actR2 = self.LiVO2_ss(y2, ybar, muR_ref, ISfuncs=None)


        return (muR1homog, muR2), (actR1, actR2)





    def LiC6_1param(self, y, ybar, muR_ref, ISfuncs=None):
        muRtheta = -self.eokT*0.12
        muRhomog = self.graphite_1param_homog_3(
            y, self.get_trode_param("Omega_a"), self.get_trode_param("Omega_b"), ISfuncs)
        muRnonHomog = self.general_non_homog(y, ybar)
        muR = muRhomog + muRnonHomog
        actR = np.exp(muR/self.T)
        muR += muRtheta + muR_ref
        return muR, actR

    def testRS(self, y, ybar, muR_ref, ISfuncs=None):
        muRtheta = 0.
        muR = self.reg_sln(y, self.get_trode_param("Omega_a"), ISfuncs)
        actR = np.exp(muR/self.T)
        muR += muRtheta + muR_ref
        return muR, actR

    def testRS_ps(self, y, ybar, muR_ref, ISfuncs=None):
        muRtheta = -self.eokT*2.
        muRhomog = self.reg_sln(y, self.get_trode_param("Omega_a"), ISfuncs)
        muRnonHomog = self.general_non_homog(y, ybar)
        muR = muRhomog + muRnonHomog
        actR = np.exp(muR/self.T)
        muR += muRtheta + muR_ref
        return muR, actR

    def LiCoO2_LIONSIMBA(self, y, ybar, muR_ref, ISfuncs=None):
        """ Torchio et al, 2016. """
        T = self.T
        r1 = 4.656
        r2 = 88.669
        r3 = 401.119
        r4 = 342.909
        r5 = 462.471
        r6 = 433.434
        r7 = 1
        r8 = 18.933
        r9 = 79.532
        r10 = 37.311
        r11 = 73.083
        r12 = 95.96
        OCV_ref = (-r1 + r2*y**2 - r3*y**4 + r4*y**6 - r5*y**8 + r6 * y**10) / \
            (-r7 + r8*y**2 - r9*y**4 + r10*y**6 - r11*y**8 + r12*y**10)
        k1 = -0.001
        k2 = 0.199521039
        k3 = -0.928373822
        k4 = 1.364550689000003
        k5 = -0.6115448939999998
        k6 = 1
        k7 = -5.661479886999997
        k8 = 11.47636191
        k9 = -9.82431213599998
        k10 = 3.048755063
        dUdT = k1*(k2+k3*y+k4*y**2+k5*y**3)/(k6+k7*y+k8*y**2+k9*y**3+k10*y**4)
        OCV = OCV_ref + dUdT*(T-1)*constants.T_ref
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def LiC6_LIONSIMBA(self, y, ybar, muR_ref, ISfuncs=None):
        """ Torchio et al, 2016. """
        T = self.T
        r1 = 0.7222
        r2 = 0.1387
        r3 = 0.029
        r4 = 0.0172
        r5 = 0.0019
        r6 = 0.2808
        r7 = 0.7984
        OCV_ref = r1 + r2*y + r3*y**0.5 - r4 * \
            y**(-1) + r5*y**(-1.5) + r6*np.exp(0.9-15*y) - r7*np.exp(0.4465*y-0.4108)
        k1 = 0.001
        k2 = 0.005269056
        k3 = 3.299265709
        k4 = -91.79325798
        k5 = 1004.911008
        k6 = -5812.278127
        k7 = 19329.7549
        k8 = -37147.8947
        k9 = 38379.18127
        k10 = -16515.05308
        k11 = 1
        k12 = -48.09287227
        k13 = 1017.234804
        k14 = -10481.80419
        k15 = 59431.3
        k16 = -195881.6488
        k17 = 374577.3152
        k18 = -385821.1607
        k19 = 165705.8597
        dUdT = k1*(k2+k3*y+k4*y**2+k5*y**3+k6*y**4+k7*y**5+k8*y**6+k9*y**7+k10*y**8) / \
            (k11+k12*y+k13*y**2+k14*y**3+k15*y**4+k16*y**5+k17*y**6+k18*y**7+k19*y**8)
        OCV = OCV_ref + dUdT*(T-1)*constants.T_ref
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR


def step_down(x, xc, delta):
    return 0.5*(-np.tanh((x - xc)/delta) + 1)


def step_up(x, xc, delta):
    return 0.5*(np.tanh((x - xc)/delta) + 1)
