"""This module handles properties associated with the active materials."""
import numpy as np

import mpet.geometry as geo
from mpet.config import constants

class sigmafuncs():

    def __init__(self, sigmafunc):
        sigmaopts = {}
        sigmaopts['constant'] = self.constant
        sigmaopts['custom1'] = self.custom1
        sigmaopts['custom1_mid'] = self.custom1_mid
        sigmaopts['custom1_low'] = self.custom1_low
        sigmaopts['custom1_mix'] = self.custom1_mix
        sigmaopts['custom_exp'] = self.custom_exp
        sigmaopts['custom_exp2'] = self.custom_exp2
        sigmaopts['custom_exp3'] = self.custom_exp3
        sigmaopts['custom_design_low'] = self.custom_design_low
        self.sigmafunc = sigmaopts[sigmafunc]

    def constant(self, y):
        return 1.

    def custom1(self, y):
        a1 =-1.11239062854913
        b1 =0.315723658193075
        c1 =0.0310390627398842
        a2 =11973.9288219987
        b2 =2.73645936763567
        c2 =0.544339695594774
        a3 =0
        b3 =-9.81089780262787
        c3 =0.00262913388145669
        a4 =7.08755376180519
        b4 =0.264458459887164
        c4 =0.126447829811828
        a5 =1.04903014525017
        b5 =0.597969984519493
        c5 =0.132082914736046
        return a1*np.exp(-((y-b1)/c1)**2) + a2*np.exp(-((y-b2)/c2)**2) +a3*np.exp(-((y-b3)/c3)**2) + a4*np.exp(-((y-b4)/c4)**2) +a5*np.exp(-((y-b5)/c5)**2)

    def custom1_mid(self, y):
        a1=0.631134107695653
        b1=0.272191928046006
        c1=0.0260856071547557
        a2=6.61092284887808
        b2=0.246113978563317
        c2=0.112403788444996
        a3=0.0
        b3=-22.3885818663665
        c3=3.06526598740341
        a4=1.07443206507807
        b4=0.576128821541168
        c4=0.225323862100736
        a5=1.13019158587323
        b5=0.369778392700667
        c5=0.0336035139775979

        return a1*np.exp(-((y-b1)/c1)**2) + a2*np.exp(-((y-b2)/c2)**2) +a3*np.exp(-((y-b3)/c3)**2) + a4*np.exp(-((y-b4)/c4)**2) +a5*np.exp(-((y-b5)/c5)**2)


    def custom1_low(self, y):
        a1 =0.614103757698879
        b1 =0.803718464626244
        c1 =0.528093418833514
        a2 =2.97161938719491
        b2 =0.319276230290737
        c2 =0.107096305077622
        a3 =0.226731861002967
        b3 =0.525996116218349
        c3 =0.0668192003512049
        a4 =5.08939503071113
        b4 =0.226451822227103
        c4 =0.0928002488435098
    #     a5 =4.33486260623102
    #     b5 =0.271934666442693
    #     c5 =0.132291340213731

        return a1*np.exp(-((y-b1)/c1)**2) + a2*np.exp(-((y-b2)/c2)**2) +a3*np.exp(-((y-b3)/c3)**2) + a4*np.exp(-((y-b4)/c4)**2)

    def custom1_mix(self, y):
        return (self.custom1(y) + self.custom1_mid(y) + self.custom1_low(y))/ 3

    def custom_exp(self, y):
        a1 =5.1161252172975
        b1 =0.216270537722703
        c1 =0.238193259582719
        a2 =8.66799996734116
        b2 =0.627303203266681
        c2 =0.400666436484809
        return a1*np.exp(-((y-b1)/c1)**2) + a2*np.exp(-((y-b2)/c2)**2)

    def custom_exp2(self, y):
        a1 =9.06100497361239
        b1 =0.439089187035436
        c1 =0.381559908648431
        a2 =1.71871983516742
        b2 =0.175441510245089
        c2 =0.165650987480684
        return a1*np.exp(-((y-b1)/c1)**2) + a2*np.exp(-((y-b2)/c2)**2)

    def custom_exp3(self, y):
        a1 =10.535144971817857
        b1 =0.814830305069119
        c1 =0.534797612892019
        a2 =4.877164282385152
        b2 =0.204530993417729
        c2 =0.209613899005338
        a3 =1.461815254718133
        b3 =0.395469849194017
        c3 =0.117833355192806
        return a1*np.exp(-((y-b1)/c1)**2) + a2*np.exp(-((y-b2)/c2)**2) + a3*np.exp(-((y-b3)/c3)**2)

    def custom_design_low(self, y):
        a0=16.2286581023396
        a1=0.365737915446761
        b1=-20.0052442755857
        a2=-15.2018714985282
        b2=-5.99054149464293
        a3=-5.04723196601268
        b3=6.84458148562386
        a4=2.72493456982476
        b4=2.73311205030866
        a5=1.01544097000561
        b5=-0.324579365923269
        w=6.55562654116442

        return a0 + a1*np.cos(y*w)+ b1*np.sin(y*w)+ a2*np.cos(2*y*w)+ b2*np.sin(2*y*w)+ a3*np.cos(3*y*w)+ b3*np.sin(3*y*w)+ a4*np.cos(4*y*w)+ b4*np.sin(4*y*w)+ a5*np.cos(5*y*w)+ b5*np.sin(5*y*w)



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
        Dopts['custom3'] = self.custom3
        Dopts['custom1_mid'] = self.custom1_mid
        Dopts['custom2_mid'] = self.custom2_mid
        Dopts['custom3_mid'] = self.custom3_mid
        Dopts['custom4_mid'] = self.custom4_mid
        Dopts['custom5_mid'] = self.custom5_mid
        Dopts['custom6_mid'] = self.custom6_mid
        Dopts['custom7_mid'] = self.custom7_mid
        Dopts['custom1_low'] = self.custom1_low
        Dopts['custom1_mix'] = self.custom1_mix
        self.Dfunc = Dopts[Dfunc]

    def constant(self, y):
        return 1.

    def lattice(self, y):
        return y*(1-y)

    def custom3(self, y):
        a0 = 1.47265551030696*10**(-15)
        a1 = 7.84011415965682*10**(-16)
        b1 = -1.08634742846511*10**(-15)
        a2 = -1.26849184549397*10**(-16)
        b2 = -3.50368625478569*10**(-16)
        w = 4.07999045920752
        return (a0 + a1*np.cos(y*w)+ b1*np.sin(y*w)+ a2*np.cos(2*y*w)+ b2*np.sin(2*y*w))/(5*10**(-16))
    def custom1_mid(self, y):
        a0 = 5.24046359930457*10**(-16)
        a1 = 6.71641548863033*10**(-17)
        b1 = -4.70437257502231*10**(-17)
        a2 = -3.25290637416018*10**(-18)
        b2 = -3.50707382157863*10**(-17)
        w = 4.00713348672168
        return (a0 + a1*np.cos(y*w)+ b1*np.sin(y*w)+ a2*np.cos(2*y*w)+ b2*np.sin(2*y*w))/(5*10**(-16))


    def custom2_mid(self, y):
        a0 = 1.02759348991012*10**(-15)
        a1 = -2.96157327252685*10**(-16)
        b1 = -8.35605061031513*10**(-16)
        a2 = -3.21777650934669*10**(-16)
        b2 = 2.92120601568872*10**(-16)
        a3 = 9.24655388670179*10**(-17)
        b3 = 3.97685114030883*10**(-17)
        w = 2.74494770955858
        return (a0 + a1*np.cos(y*w)+ b1*np.sin(y*w)+ a2*np.cos(2*y*w)+ b2*np.sin(2*y*w) + a3*np.cos(3*y*w)+ b3*np.sin(3*y*w))/(5*10**(-16))

    def custom3_mid(self, y):

        p1 =34.2286369254162
        p2 =-88.8027328121971
        p3 =93.3124628565447
        p4 =-54.7228246554853
        p5 =21.6595614609822
        p6 =-6.41408107148133
        p7 =1.35196667363861
        p8 =-0.19128969116873
        p9 =-15.306847193439
        return (10**(p1*y**8 + p2*y**7 + p3*y**6 + p4*y**5 + p5*y**4 + p6*y**3 + p7*y**2 + p8*y + p9)) /(5*10**(-16))

    def custom4_mid(self, y):
        p1 =-15.5340908523484
        p2 =33.160056744686
        p3 =-13.6443943716574
        p4 =-15.6995792032625
        p5 =15.1907969904133
        p6 =-3.61169759197661
        q1 =-2.13096823891331
        q2 =0.867824275098675
        q3 =1.0218735602065
        q4 =-0.983204077681832
        q5 =0.233439498299931
        return 10**((p1*y**5 + p2*y**4 + p3*y**3 + p4*y**2+ p5*y+ p6)/(y**5 + q1*y**4 + q2*y**3 + q3*y**2+ q4*y+ q5)) /(5*10**(-16))

    def custom5_mid(self, y):
        p1 =-15.3839910572732
        p2 =-22.3498341092738
        p3 =2.25918213486115
        p4 =39.9344135710357
        p5 =17.6206074995782
        p6 =-27.1242374271543
        q1 =1.45502681505354
        q2 =-0.131778063592997
        q3 =-2.61491426340857
        q4 =-1.15019355512244
        q5 =1.77066594929798
        return 10**((p1*y**5 + p2*y**4 + p3*y**3 + p4*y**2+ p5*y+ p6)/(y**5 + q1*y**4 + q2*y**3 + q3*y**2+ q4*y+ q5)) /(5*10**(-16))

    def custom6_mid(self, y):
        p1 =-15.534284568355782
        p2 =-11.554870941970011
        p3 =4.918792282832964
        p4 =18.324321469185563
        p5 =18.599864249560195
        p6 =-20.998099023131541
        q1 =0.791373664318053
        q2 =-0.339326895555531
        q3 =-1.198292436203694
        q4 =-1.230170377304700
        q5 =1.383417091213740
        return 10**((p1*y**5 + p2*y**4 + p3*y**3 + p4*y**2+ p5*y+ p6)/(y**5 + q1*y**4 + q2*y**3 + q3*y**2+ q4*y+ q5)) /(5*10**(-16))

    def custom7_mid(self, y):
        p1 =-15.389158274806022
        p2 =31.140004480433017
        p3 =-8.302705179538862
        p4 =-24.442830384447589
        p5 =22.624127792999833
        p6 =-5.930495796561357
        q1 =-2.015402956660859
        q2 =0.510831465873819
        q3 =1.627178795342191
        q4 =-1.493814996264665
        q5 =0.390812139966707

        return 10**((p1*y**5 + p2*y**4 + p3*y**3 + p4*y**2+ p5*y+ p6)/(y**5 + q1*y**4 + q2*y**3 + q3*y**2+ q4*y+ q5))

    def custom1_low(self, y):
        a0 = 1.14907525660876*10**(-15)
        a1 = 4.48545291585507*10**(-16)
        b1 = -7.05449089826032*10**(-16)
        a2 = -9.74118198082059*10**(-17)
        b2 = -1.4268753889296*10**(-16)
        w = 3.95168887243999
        return (a0 + a1*np.cos(y*w)+ b1*np.sin(y*w)+ a2*np.cos(2*y*w)+ b2*np.sin(2*y*w))/(5*10**(-16))

    def custom1_mix(self, y):
        return (self.custom3(y) + self.custom1_mid(y) + self.custom1_low(y))/ 3

class kLifuncs():

    def __init__(self, kLifunc):
        kLiopts = {}
        kLiopts['constant'] = self.constant
        kLiopts['custom_mid'] = self.custom_mid
        kLiopts['custom2_mid'] = self.custom2_mid
        kLiopts['custom5_mid'] = self.custom5_mid
        kLiopts['custom6_mid'] = self.custom6_mid
        kLiopts['custom7_mid'] = self.custom7_mid
        self.kLifunc = kLiopts[kLifunc]

    def constant(self, y):
        return 1

    def custom_mid(self, y):
        a = 0.990754082767668
        b = 9.60529528427668
        c =-16.7458287647686
        m = -28.4923425695893
        n = 0.501544182876847
        r = -1.77212306347445

        return 10**(2.6/(a+np.exp(-b*y+c)) -2.2/(n+np.exp(-m*y+r)))

    def custom2_mid(self, y):
        a = 1.02942500588641
        b = 8.99863845919693
        c =-16.2509939788027
        m =-25.3929903996158
        n = 0.348721104422145
        r = -1.41423445971999

        return 10**(3/(a+np.exp(-b*y+c)) -2.2/(n+np.exp(-m*y+r)))

    def custom5_mid(self, y):
        a = 0.898307777851485
        b = 36.6532681948738
        c =2.119116695964
        m =0.0820834543614888
        n =4.35293222785387
        r =-0.338036018293395

        return 10**(3/(a+np.exp(-b*y+c)) -2.2/(n+np.exp(-m*y+r)))

    def custom6_mid(self, y):
        a = 0.998991455415107
        b = 8.670386813179991
        c = -25.942946335557011
        m =-31.448461485490228
        n = 0.311387105405065
        r = -2.623498506020950

        return 10**(2.6/(a+np.exp(-b*y+c)) -1.2/(n+np.exp(-m*y+r)))

    def custom7_mid(self, y):
        a = 0.785249980029520
        b = 6.555323684927352
        c = -16.564952401956099
        m = -34.630433392544049
        n = 0.281738006742055
        r = -3.045705181603137

        return 10**(2.6/(a+np.exp(-b*y+c)) -1.2/(n+np.exp(-m*y+r)))



class Rfilmfuncs():

    def __init__(self, Rfilmfunc):
        Rfilmopts = {}
        Rfilmopts['constant'] = self.constant
        Rfilmopts['custom_mid'] = self.custom_mid
        Rfilmopts['custom2_mid'] = self.custom2_mid
        Rfilmopts['custom3_mid'] = self.custom3_mid
        Rfilmopts['custom4_mid'] = self.custom4_mid
        self.Rfilmfunc = Rfilmopts[Rfilmfunc]

    def constant(self, y):
        return 1

    def custom_mid(self, y):
        a0 = 0.196398223205289
        a1 = -0.219593351032325
        b1 = -0.235110246185945
        a2 = -0.011549382462294
        b2 = 0.172539926621049
        a3 = 0.0430690080295571
        b3 = -0.0356883957437408
        a4 = -0.00829710911565918
        b4 = -0.000887770550027742
        w = 3.32654876491931
        return a0 + a1*np.cos(y*w)+ b1*np.sin(y*w)+ a2*np.cos(2*y*w)+ b2*np.sin(2*y*w) + a3*np.cos(3*y*w)+ b3*np.sin(3*y*w) + a4*np.cos(4*y*w)+ b4*np.sin(4*y*w)

    def custom2_mid(self, y):
        a0 = -4.45707464610095
        a1 = -0.314587327364548
        b1 = -0.702461771788783
        a2 = -0.153754085249938
        b2 = -0.0567509023480456
        a3 = -0.156339572248037
        b3 = 0.046764255783031
        a4 = -0.0327049839944974
        b4 = -0.00690390154452286
        a5 = -0.00867054907535169
        b5 = 0.0581858856748892
        w = 8.79044332072256
        return 10**(a0 + a1*np.cos(y*w)+ b1*np.sin(y*w)+ a2*np.cos(2*y*w)+ b2*np.sin(2*y*w) + a3*np.cos(3*y*w)+ b3*np.sin(3*y*w) + a4*np.cos(4*y*w)+ b4*np.sin(4*y*w)+ a5*np.cos(5*y*w)+ b5*np.sin(5*y*w))

    def custom3_mid(self, y):
        p1 =-2.71649169521244
        p2 =0.538469290139425
        p3 =2.64314601040722
        p4 =-2.16378616733113
        p5 =0.660411476458034
        p6 =-0.0734683342209862
        q1 =-0.692741838529293
        q2 =-0.221751923738117
        q3 =0.346845357296604
        q4 =-0.118575048968982
        q5 =0.0137143584680668


        return 10**((p1*y**5 + p2*y**4 + p3*y**3 + p4*y**2+ p5*y+ p6)/(y**5 + q1*y**4 + q2*y**3 + q3*y**2+ q4*y+ q5))

    def custom4_mid(self, y):
        p1 = -2.359053359214419
        p2 = 0.753900208305855
        p3 = 1.597847465607573
        p4 = -1.261476708681014
        p5 = 0.344672107553635
        p6 = -0.033676396760198
        q1 = -0.864844823595020
        q2 = 0.059040778287494
        q3 = 0.154800174506188
        q4 = -0.056988627978934
        q5 = 0.006209820678795

        return 10**((p1*y**5 + p2*y**4 + p3*y**3 + p4*y**2+ p5*y+ p6)/(y**5 + q1*y**4 + q2*y**3 + q3*y**2+ q4*y+ q5))




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

    def LiCF4_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """
        CFx (x=1) discharge OCV from Gomadam Comsol
        """
        goma_coms = (3.426177769+1317.958327*y+9596.394751*y**2-12208.74244*y**3)/(1+527.3952345*y+3570.78942*y**2-4231.792944*y**3-361.6223807*y**4)

        # OCV = goma_coms - (1/0.57)*((310.15*1.38e-23)/(1.602e-19))*np.log((3*1.908/(3.8277e-9*2.75e6)) * 1e-7 * y * (1-y)**(1/3))
        OCV = 0+goma_coms - (1/0.57)*((310.15*1.38e-23)/(1.602e-19))*np.log((3*1.908/(1e-6*2.75e6)) * 1e-5 * y * (1-y)**(1/3))
        # OCV = goma_coms - (1/0.57)*((310.15*1.38e-23)/(1.602e-19))*np.log((3*1.908/(self.get_trode_param("mean2")*2.75e6)) * self.get_trode_param("k2") * y * (1-y)**(1/3))


        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        # EoKT = constants.e / (constants.k * 310.15)
        # muR = -EoKT*OCV + muR_ref
        # actR = None
        return muR, actR

    def LiCF5_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """
        CFx (x=1) discharge OCV from Gomadam Comsol
        """
        goma_coms = (3.426177769+1317.958327*y+9596.394751*y**2-12208.74244*y**3)/(1+527.3952345*y+3570.78942*y**2-4231.792944*y**3-361.6223807*y**4)

        # OCV = goma_coms - (1/0.57)*((310.15*1.38e-23)/(1.602e-19))*np.log((3*1.908/(3.8277e-9*2.75e6)) * 1e-7 * y * (1-y)**(1/3))
        OCV = + 0.01 - 0.08*(0.18-y)+goma_coms - (1/0.57)*((310.15*1.38e-23)/(1.602e-19))*np.log((3*1.908/(1e-6*2.75e6)) * 1e-5 * y * (1-y)**(1/3))
        # OCV = goma_coms - (1/0.57)*((310.15*1.38e-23)/(1.602e-19))*np.log((3*1.908/(self.get_trode_param("mean2")*2.75e6)) * self.get_trode_param("k2") * y * (1-y)**(1/3))


        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        # EoKT = constants.e / (constants.k * 310.15)
        # muR = -EoKT*OCV + muR_ref
        # actR = None
        return muR, actR

    def LiSi_DeLi_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """
        Si Delithiation 2023
        """
        a =-0.00588731077882632
        b =0.948060714071423
        bb =0.362020695953494
        c =-1.0929719146036
        cc =0.230243434621118
        d =2.88578320867775
        dd =-2.02658195341187
        e = -1.66957510133694
        ee =1.56813539818296
        f = -2.13252949510444
        ff =1.18127524811439
        g =0.528989739814381
        gg =1.04637994155614
        h =1.8952245260347
        hh =-2.24873010707051
        k =-0.509263799951346
        OCV =  a*np.log(y/(1-y)) + b + (c*y + d*y**2+ e*y**3+ f*y**4+ g*y**5+ h*y**6 +
                         k*y**7)/(bb + cc*y + dd*y**2+ ee*y**3+ ff*y**4+ gg*y**5+hh*y**6)

        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        # EoKT = constants.e / (constants.k * 310.15)
        # muR = -EoKT*OCV + muR_ref
        # actR = None
        return muR, actR

    def LiSi_Li_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """
        Si Lithiation 2023
        """
        a =-0.0836650628158838
        b =0.284240463531336
        bb =0.00734884061487265
        c =0.0216516264965879
        cc =0.13124492283084
        d =-0.710622455588877
        dd =1.1580840418446
        e = 2.67311331694277
        ee =-1.12044805260863
        f = -3.7620463488108
        ff =-0.289863721558452
        g =0.24585013210871
        gg =0.789899264851894
        h =3.58781848387559
        hh =-0.656888102352765
        k =-2.04969882841769
        OCV =  a*np.log(y/(1-y)) + b + (c*y + d*y**2+ e*y**3+ f*y**4+ g*y**5+ h*y**6 +
                         k*y**7)/(bb + cc*y + dd*y**2+ ee*y**3+ ff*y**4+ gg*y**5+hh*y**6)
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

    def graphite_1param_homog_Liang(self, y, Omga, Omgb, ISfuncs=None):
        """ 2023 Si/C"""
        width = 5e-2
        tailScl = 5e-2
        muLtail = -tailScl*1./(y**(0.55))
        muRtail = 0.7e1*step_up(y, 1.0, 0.015)
        muLlin = (0.15*Omga*12*(0.17-y**0.98)
                  * step_down(y, 0.55, 0.9*width)*step_up(y, 0.38, width))
        muRlin = (0.1*Omga*4*(0.74-y) + 0.55*Omgb - 2*(1-y))*step_up(y, 0.6, 0.8*width)
        muLMod = (0.
                  + (30*(-np.exp(-y/0.025))- 2*(1-y))*step_down(y, 0.035, width)
                  + 0.7*(np.tanh((y-0.37)/0.075) - 1)
                  + 0.8*(np.tanh((y-0.2)/0.06) - 1)
                  + 0.38*(np.tanh((y-0.14)/0.015) - 1)
                  )*step_down(y, 0.42, width)
        muR = -0.02 + muLMod + muLtail + muRtail + muLlin + muRlin
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
        muRnonHomog = self.general_non_homog(y, ybar)
        muR = muRhomog + muRnonHomog
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
        y1_bar, y2_bar = ybar
        # muRtheta2 = -self.eokT*2.59

        muR1, actR1 = self.LiAgO(y1, y1_bar, muR_ref, ISfuncs=None)

        muR2, actR2 = self.LiVO2_ss(y2, y2_bar, muR_ref, ISfuncs=None)


        return (muR1, muR2), (actR1, actR2)





    def LiC6_1param(self, y, ybar, muR_ref, ISfuncs=None):
        muRtheta = -self.eokT*0.12
        muRhomog = self.graphite_1param_homog_3(
            y, self.get_trode_param("Omega_a"), self.get_trode_param("Omega_b"), ISfuncs)
        muRnonHomog = self.general_non_homog(y, ybar)
        muR = muRhomog + muRnonHomog
        actR = np.exp(muR/self.T)
        muR += muRtheta + muR_ref
        return muR, actR

    def LiC6_Liang_1param(self, y, ybar, muR_ref, ISfuncs=None):
        muRtheta = -self.eokT*0.12
        muRhomog = self.graphite_1param_homog_Liang(
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
