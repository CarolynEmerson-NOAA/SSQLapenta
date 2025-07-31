import numpy as np
import numpy.ma as ma

p_0 = 1000
Rd = 287 # J/K Kg
Cp = 1004 # J/K

e_0 = 6.1173 # std atm vapor pressure (hPa)
Lv = 2501000 # J/Kg
Rv = 461.5 # J/K Kg
t_0 = 273.16 # K

def calc_t(th, p):
    """
    Calculates temperature for a np.ndarray

    Adapted from similar wrftools routines (https:https://github.com/keltonhalbert/wrftools/)
    Copyright (c) 2015, Kelton Halbert

    Input:

        th - np.ndarray of potential temperature (K)
        p - np.ndarray of pressure (hPa)

    Returns:

        th - np.ndarray of temperature (K)
    """
    t = th * (p / p_0)**(Rd / Cp)

    nans = np.isnan(t)
    t[nans] = 0.

    return t


def calc_thv(th, qv, qt):
    """
    Calculate virtual temperature (or virtual potential temperature)

    Adapted from simialr wrftools routine (https:https://github.com/keltonhalbert/wrftools/)
    Copyright (c) 2015, Kelton Halbert

    Args:
    -----------------------------
        th, numpy array, Temperature or Potential Temperature (in K)
        qt, numpy array, total water mixing ratio (in Kg/Kg)
        qv, numpy arrau, Water vapor mixing ratio (in Kg/Kg)

    Returns:
    ------------------------------
        thv, numpy array, virtual (potential) temperature (in K)

    """
    thv = th * (1. + 0.61 * qv - qt)
    thv[np.isnan(thv)] = 0.

    return thv

def relative_humidity(t, p, qv):
    """ Returns relative humidity computed from the saturation mixing ratio

    Args:
        t, np.array
            Temperature (in K)
        p, np.array
            Pressure (in hPa)
        qv, np.array
            Mixing ratio (Kg/Kg)
    """
    e_s = e_0 * np.exp((Lv / Rv) * ((1. / t_0) - (1. / t)))  #sat vapor pressure via Clasius Clapeyron equation (hPa)
    w_s = (0.622 * e_s) / (p - e_s) #sat mixing ratio (Kg/Kg)
    rh = (qv / w_s) * 100. #relative humidity

    return rh




def calc_visibility(p,t,u,v,qv,qc,qr,qi,qs,qg,czen):

### INFO
#   Started with Stoelinga-Warner algorithm for hydrometeors only.
#   Added coefficients for graupel.
#    Added algorithm for clear-air RH-based visibility.
#
#   This routine computes horizontal visibility (in km) at the
#   surface or lowest model layer, from qc, qr, qi, qs, and qg.
#   qv--water vapor mixing ratio (kg/kg)
#   qc--cloud water mixing ratio (kg/kg)
#   qr--rain water mixing ratio  (kg/kg)
#   qi--cloud ice mixing ratio   (kg/kg)
#   qs--snow mixing ratio        (kg/kg)
#   qg--graupel mixing ratio     (kg/kg)
#   u/v - u/v wind components    (m/s)
#   tb --            temperature (k)
#   pp--pressure                 (Pa)
#   rhb-- relative humidity      (0-100%)
#   aextc55--aerosol extinction coefficient (m**-1)
#
#
#   Independent of the above definitions, the scheme can use different
#   assumptions of the state of hydrometeors:
#        meth='r': Uses the four mixing ratios qrain, qsnow, qclw,
#           and qclice
#
#   The routine uses the following
#   expressions for extinction coefficient, beta (in km**-1),
#   with C being the mass concentration (in g/m**3):
#
#      cloud water:  beta = 144.7 * C ** (0.8800)
#      rain water:   beta =  2.24 * C ** (0.7500)
#      cloud ice:    beta = 327.8 * C ** (1.0000)
#      snow:         beta = 10.36 * C ** (0.7776)
#      graupel:      beta =  8.0  * C ** (0.7500)
#
#   These expressions were obtained from the following sources:
#
#      for cloud water: from Kunkel (1984)
#      for rainwater: from M-P dist'n, with No=8e6 m**-4 and
#         rho_w=1000 kg/m**3
#      for cloud ice: assume randomly oriented plates which follow
#         mass-diameter relationship from Rutledge and Hobbs (1983)
#      for snow: from Stallabrass (1985), assuming beta = -ln(.02)/vis
#      for graupel: guestimate by John Brown and Stan Benjamin,
#         similar to snow, but a smaller extinction coef seemed
#         reasonable.  27 Aug 99
#
#   calculated, and then all applicable betas are summed to yield
#   a single beta. Then the following relationship is used to
#   determine visibility (in km), where epsilon is the threshhold
#   of contrast, usually taken to be .02:
#
#      vis = -ln(epsilon)/beta      [found in Kunkel (1984)]
#
#   The 'aextc55' field is 3-D and is derived from the 'aod_3d' field
#   by dividing by dz, the vertical thickness of that model level in 'm'.
#   This can be handled as a 2-D field if needed to save resources.
               #Method used for clear-air visibility with extension for aerosols
  method = 1
#                RH-only method (1),
#                Aerosol method (2),
#                Smoke added to RH method for clear air  (3)
#                   3 - option to add reducted visibility from smoke-based aerosols.


  CELKEL     = 273.15
  TICE       = CELKEL-10.
  COEFLC     = 144.7
  COEFLP     =   2.24
  COEFFC     = 327.8
  COEFFP     =  10.36

# - modified number - Stan B. - Dec 2007
#     after quick talks with John Brown and Ismail Gultepe
  coeffp_dry =  10.0
  coeffp_wet =   6.0
#COEFFg     =   8.0
# - values from Roy Rasmussen - Dec 2003
#    Rasmussen et al. 2003,  J. App. Meteor.
#    Snow Nowcasting Using a Real-Time Correlation of Radar Reflectivity with Snow Gauge Accumulation
  coeffg     =   4.0
  EXPONLC    =   0.8800
  EXPONLP    =   0.7500
  EXPONFC    =   1.0000
#    EXPONFP    =   0.7776
# - new value from Roy Rasmussen - Dec 2003
  EXPONFP    =   1.0
  EXPONFg    =   0.75
#  CONST1=-np.log(0.02)
#if(MODELNAME == 'RAPR'):
  CONST1= 3.000

# Total visibility is minimum of vis-rh  (developed by Benjamin, Brown, Smirnova)
#   and vis-hydrometeors from Stoelinga/Warner

  RHOICE=917.
  RHOWAT=1000.

  vis_min = 1.e6
  visrh_min = 1.e6
  coef_snow = 0.0

  nx = p.shape[2]
  ny = p.shape[1]
  nz = p.shape[0]

  # CALCULATE STUFF
  qtotal = qr+qc+qg+qs+qi
  tv = calc_thv(t, qv, qtotal)
  tk = calc_t(t, p/100)
  rh = np.ma.masked_invalid(relative_humidity(t, p/100, qv))
#  aextc55 = aod3d[0:nz-1,:,:]/ph[1:nz,:,:]

  out_vis=np.zeros(shape=(3,ny,nx))
  #visibility=np.zeros(ny,nx)

  for j in range(0, ny):
   for i in range(0, nx):

      #  - take max RH of levels 1 and 2 near the sfc
      aaa=np.float64(rh[1,j,i])
      bbb=np.float64(rh[2,j,i])
      #maxrh= np.max(aaa,bbb)
      maxrh=max(aaa,bbb)
      qrh = max(0.0,min(0.8,(maxrh/100.-0.15)))
      visrh = 90. * np.exp(-2.5*qrh)
      #visrh = 60. * np.exp(-2.5*qrh)

      #CALCUALTE VISIBILTY FROM FIRST THREE LEVELS
      for k in range(0, 3):

#  -- add term to increase RH vis term for
#     low-level wind shear increasing from 4 to 6 ms-1
#     (using Evan Kuchera's paper as a guideline)

# -- calculate term for shear between levels 1 and 4
#   (about 25 hPa for HRRR/RAP)
        shear = np.sqrt( (u[4,j,i]-u[k,j,i])**2 + (v[4,j,i]-v[k,j,i])**2  )

        shear_fac = min(1.,max(0.,(shear-4.)/2.) )
        if (visrh < 10.): visrh = visrh + (10.-visrh)*shear_fac

        RHOAIR=p[k,j,i]/(Rd*tv[k,j,i])
        VOVERMD=(1.+qv[k,j,i])/RHOAIR+(qc[k,j,i]+qr[k,j,i])/RHOWAT+(qg[k,j,i]+qi[k,j,i]+qs[k,j,i])/RHOICE

        CONCLC=qc[k,j,i]/VOVERMD*1000.
        CONCLP=qr[k,j,i]/VOVERMD*1000.
        CONCFC=qi[k,j,i]/VOVERMD*1000.
        CONCFP=qs[k,j,i]/VOVERMD*1000.
        CONCFg=qg[k,j,i]/VOVERMD*1000.

        temp_fac = min(1.,max((tk[k,j,i]-271.15),0.) )
        coef_snow = coeffp_dry*(1.-temp_fac)+ coeffp_wet* temp_fac

# Key calculation of attenuation from each hydrometeor type (cloud, snow, graupel, rain, ice)
        BETAV = COEFFC*CONCFC**EXPONFC + coef_snow*CONCFP**EXPONFP + COEFLC*CONCLC**EXPONLC + COEFLP*CONCLP**EXPONLP + coeffg*CONCFg**EXPONFg +1.E-10

# Addition of attenuation from aerosols if option selected
        if(method == 2 or method == 3):   # aerosol method
            BETAV = BETAV + np.abs(aextc55[k,j,i])*1000.
            #print(i,j,k,aextc55[k,j,i],ph[k,j,i])

# Calculation of visibility based on hydrometeor and aerosols.  (RH effect not yet included.)
        #out_vis[k,j,i]=min(90.,CONST1/BETAV+extcof55[k,j,i])      # max of 90km
        out_vis[k,j,i]=min(90.,CONST1/BETAV)

        if (out_vis[k,j,i] < vis_min): vis_min = out_vis[k,j,i]
        if (visrh < visrh_min):        visrh_min = visrh

#-- Dec 2003 - Roy Rasmussen (NCAR) expression for night vs. day vis1.609 factor is number of km in mile.
        vis_night = 1.69 * ((out_vis[k,j,i]/1.609)**0.86) * 1.609
        zen_fac = min(0.1,max(czen[j,i],0.))/ 0.1
        out_vis[k,j,i] = zen_fac * out_vis[k,j,i] + (1.-zen_fac)*vis_night

        #print(i,j,k,out_vis[k,j,i])

        if(method == 1 or method == 3): # RH method (if lower vis)
         out_vis[k,j,i] = min(out_vis[k,j,i],visrh)

  visibility = out_vis[1,:,:]
  # CONVERT FROM KM to Miles
  visibility[visibility>24.135] = 24.135
  visibility = visibility*0.621371

  return visibility
