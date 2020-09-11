# DMCA-based methodology to estimate catchment response time:
#this routine finds the catchment response time given a rainfall 
#and a streamflow timeseries aligned in time and with the same
#length and time step.

#HOW TO CALL THE FUNCTION FROM ANOTHER SCRIPT
#from Pyfun_Tr_DMCA import Tr_DMCA

#INPUT
#rain_array: rainfall timeseries in a 0 x n numpy array 
#flow_array: streamflow timeseries in a 0 x n numpy array 
#max_window: maximum window tested. Set it sensibly according to the resolution of your data (e.g. hourly data, mac_window=300 means that Tc can be maximum 300/2=150 hours ~ 6days)


#function
def Tr_DMCA(rain_array, flow_array, max_window):
    import numpy as np
    import math

    #ROUTINE
    rain_int=np.nancumsum(rain_array, axis=0) #cumulating rainfall timeseries (Eq.1)
    flow_int=np.nancumsum(flow_array, axis=0) #cumulating streamflow timeseries (Eq.2)
    T=len(rain_array) #length of the timeseries
    n=0 #counter for the windows tested
    rho=np.empty(int(max_window/2))*np.nan #initializing rho coefficient as a NaN vector
    
    for window in range(3, max_window,2):
        rain_mean= np.convolve(rain_int, np.ones(window), 'valid') / window #moving average on the integrated rainfall timeseries (Eq.5)
        flow_mean= np.convolve(flow_int, np.ones(window), 'valid') / window #moving average on the integrated streamflow timeseries (Eq.6)
        flutt_rain=rain_int[int(float(window)/2+0.5-1):int(len(rain_int)-float(window)/2+0.5)]-rain_mean
        F_rain=(1/float(T-window+1))*np.nansum((np.square(flutt_rain))) #Squared rainfall fluctuations (Eq.3)
        flutt_flow=flow_int[int(float(window)/2+0.5-1):int(len(flow_int)-float(window)/2+0.5)]-flow_mean
        F_flow=(1/float(T-window+1))*np.nansum((np.square(flutt_flow))) #Squared streamflow fluctuations (Eq.4)
        F_rain_flow=(1/float(T-window+1))*np.nansum(np.multiply(flutt_rain,flutt_flow)) #Bivariate rainfall-streamflow fluctuations (Eq.7)
        if np.logical_or(F_rain==0 ,F_rain==0):
            rho[n]=np.nan #avoiding division by 0
        else:
            rho[n]=F_rain_flow/(math.sqrt(F_rain)*math.sqrt(F_flow)) #DMCA-based correlation coefficent (Eq.8)
        n=n+1
        
    #OUTPUT
    position_minimum=np.where(rho==np.nanmin(rho))
    catchment_response_time=float(position_minimum[0][0])+1
    return catchment_response_time;      