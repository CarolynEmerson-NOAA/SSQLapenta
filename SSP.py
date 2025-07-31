import numpy as np
import metpy.calc as mpcalc
from metpy.units import units

#Snow Squall Parameter in Banacos et. al 2014
def build_highres(hgt2km,p2km,t2km,td2km,rh2km,wind2km):
    #Get High Resolution Data to Exaclty 2km
    H_new=np.arange(0,2010,10) #interpolate from 0m (surface) to 2000m or 2km
    p_highres=np.interp(H_new,hgt2km,p2km)
    t_highres=np.interp(H_new,hgt2km,t2km)
    td_highres=np.interp(H_new,hgt2km,td2km)
    rh_highres=np.interp(H_new,hgt2km,rh2km)
    hgt_highres=H_new
    wind_highres=np.interp(H_new,hgt2km,wind2km)
    return p_highres,t_highres,td_highres,rh_highres,hgt_highres,wind_highres

def wetbulb_calc(t_sfc,td_sfc,p_sfc):
    L=2.5e6
    es=611*np.exp(L/461.5*(1/273.16-1/t_sfc))
    zlcl=125*(t_sfc-td_sfc)
    tlcl=t_sfc-9.8/1000*zlcl
    rs=0.622*es/(p_sfc-es)
    Gs=9.8065/1004.7*(1+rs*L/287.0/t_sfc)/(1+L**2*rs*0.622/1004.7/287/t_sfc**2)
    wetbulb_val = tlcl+Gs*zlcl
    return wetbulb_val

def snow_squall_param(p_3d,t_3d,td_3d,rh_3d,hgt_3d,wind_3d,
        msl,p_sfc,t_sfc,td_sfc,rh_sfc,wind_sfc):

    count = 0.
    snsq = np.zeros([np.shape(p_sfc)[0],np.shape(p_sfc)[1]]) #2D aray for surface
    snsq = np.full_like(snsq,np.nan) #fill with nans for processing

    for lvl in range(p_3d.shape[0]):#loop through the atmosphere to build the 2km layers
        count +=1
        #Find where the nans are in the surface 2D array. If nans exists, then it must be underground in some locations
        nan_values = np.where(np.isnan(snsq))

        if not nan_values[0].size == 0: #if Nans exists still... lets try to calculate numbers

            ##########################Begin Calculations##################################
            for i,j in zip(nan_values[0],nan_values[1]): #Only do calculations in grid locations of nans. (i.e. where calculations haven't been done yet)

                agl = hgt_3d[lvl,i,j] - msl[i,j] #check to see if you are above ground. Otherwise skip this grid point till later.

                if agl > 0: #only do the calculation if you are above ground (i.e. postive heights)
                    #Check Surface Wetbulb Temperature
                    wetbulb_thresh = wetbulb_calc(t_sfc[i,j],td_sfc[i,j],p_sfc[i,j])
                    if wetbulb_thresh <= 274.15: #If wetbulb is above 1°C then skip calculation

                        for k in range(np.shape(hgt_3d)[0]): #loop up to 2km
                            k +=lvl #Start at the surface level and go up. Be sure you have enough vertical data or this will crash

                            agl_2kmlayer = hgt_3d[k,i,j] - msl[i,j] #Find how high off the ground you are
                            if agl_2kmlayer > 2000: #Once we find 2km... run calculations over the layer
                                #Build 2km arrays
                                k_end = k+1 #mark the indice to stop at. (i.e. just past the 2km indice)
                                p2km = p_3d[lvl:k_end,i,j] #create a new 1D array with just 2km data
                                p2km = np.insert(p2km,0,p_sfc[i,j]) #put the surface data in the front of the array
                                t2km = t_3d[lvl:k_end,i,j] #create a new 1D array with just 2km data
                                t2km = np.insert(t2km,0,t_sfc[i,j]) #put the surface data in the front of the array
                                td2km = td_3d[lvl:k_end,i,j] #create a new 1D array with just 2km data
                                td2km = np.insert(td2km,0,td_sfc[i,j]) #put the surface data in the front of the array
                                rh2km = rh_3d[lvl:k_end,i,j] #create a new 1D array with just 2km data
                                rh2km = np.insert(rh2km,0,rh_sfc[i,j]) #put the surface data in the front of the array
                                hgt2km = hgt_3d[lvl:k_end,i,j] #create a new 1D array with just 2km data
                                hgt2km = np.insert(hgt2km,0,0) #put the surface data in the front of the array 0-meters = surface
                                wind2km = wind_3d[lvl:k_end,i,j] #create a new 1D array with just 2km data
                                wind2km = np.insert(wind2km,0,wind_sfc[i,j]) #put the surface data in the front of the array

                                #send to Numba function to speed up interpolation
                                p_highres,t_highres,td_highres,rh_highres,hgt_highres,wind_highres=build_highres(hgt2km,p2km,t2km,td2km,rh2km,wind2km)

                                #Use MetPy to calculate theta_e at surface. The magnitude function strips the units off the number
                                sfc_theta_e = (mpcalc.equivalent_potential_temperature(((p_sfc[i,j]/100)*units.hPa),(t_sfc[i,j]*units.degK),(td_sfc[i,j]*units.degK))).magnitude
                                theta_e_2km = (mpcalc.equivalent_potential_temperature(((p_highres[-1]/100)*units.hPa),(t_highres[-1]*units.degK),(td_highres[-1]*units.degK))).magnitude

                                #Calculate SNSQ Terms
                                term1 = (np.mean(rh_highres)-60)/15 #Relative Humidity Term
                                term2 = (4 -(theta_e_2km-sfc_theta_e))/4 #Equivalent Potential Temperature Term
                                term3 = np.mean(wind_highres)/9 #Wind term

                                if term1 < 0 or term2 < 0: #If terms go negative set it all to 0
                                    snsq[i,j] = 0
                                else:
                                    snsq[i,j] = term1*term2*term3 #Final calculation

                                break #once we find 2km we are done with this grid point
                    else: #mark the grid points we want to mask out
                        snsq[i,j] = -99
        else: #If nans no longer exist. Don't calculate anything. We have found the surface at all grid points and finished calculations
            print("Too high in the Atmosphere")
            break

    snsq_grid = np.where(snsq < -1, np.nan, snsq) #put nans back in for -99
    return snsq_grid
