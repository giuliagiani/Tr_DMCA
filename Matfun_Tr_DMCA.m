function catchment_response_time=Matfun_Tr_DMCA(rain, flow, max_window)
%DMCA-based methodology to estimate catchment response time:
%this function finds the catchment response time given a rainfall 
%and a streamflow timeseries aligned in time and with the same
%length and time step.

% INPUT
%max_window: Maximum window tested. Set it sensibly according to the resolution of your data (e.g. hourly data, max_window=300 means that time of concentration can be maximum 300hours/2 = 150hours =~ 6days)
%rain: your rainfall timeserie as a row vector 1 x n
%flow: your streamflow timeseries as a row vector 1 x n


%ROUTINE
rain_int=cumsum(rain, 'omitnan'); %cumulating rainfall timeseries (Eq.1)
flow_int=cumsum(flow, 'omitnan'); %cumulating streamflow timeseries (Eq.2)
T=length(rain); %length of the timeseries

for window=3:2:max_window
    rain_mean((window-1)/2,:)=movmean(rain_int, window); %moving average on the integrated rainfall timeseries (Eq.5)
    flow_mean((window-1)/2,:)=movmean(flow_int, window); %moving average on the integrated streamflow timeseries (Eq.6)
    flutt_rain((window-1)/2,:)=rain_int-rain_mean((window-1)/2,:);
    F_rain((window-1)/2)=(1/(T-window+1))*nansum((flutt_rain((window-1)/2,window-0.5*(window-1):T-0.5*(window-1))).^2); %Squared rainfall fluctuations (Eq.3)
    flutt_flow((window-1)/2,:)=flow_int-flow_mean((window-1)/2,:);
    F_flow((window-1)/2)=(1/(T-window+1))*nansum((flutt_flow((window-1)/2,window-0.5*(window-1):T-0.5*(window-1))).^2); %Squared streamflow fluctuations (Eq.4)
    F_rain_flow((window-1)/2)=(1/(T-window+1))*nansum(flutt_rain((window-1)/2,window-0.5*(window-1):T-0.5*(window-1)).*flutt_flow((window-1)/2,window-0.5*(window-1):T-0.5*(window-1))); %Bivariate rainfall-streamflow fluctuations (Eq.7)
    rho((window-1)/2)=F_rain_flow((window-1)/2)/(sqrt(F_rain((window-1)/2))*sqrt(F_flow((window-1)/2))); %DMCA-based correlation coefficent (Eq.8)
end


%OUTPUT
position_minimum=find(rho==nanmin(rho));
catchment_response_time=position_minimum;

