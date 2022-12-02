function pos = getStationPos(t, stationNum)
% tracking station locations (need to pass in t and n (station number)
% xi = @(t,n)rE*cos(omegaE*t + (n-1)*pi/6);
% yi = @(t,n)rE*sin(omegaE*t + (n-1)*pi/6);
% xidot = -rE*omegaE*sin(omegaE*t + (n-1)*pi/6);
% yidot = rE*omegaE*cos(omegaE*t + (n-1)*pi/6);
% pos = [Xi; Xidot; Yi; Yidot]  
    % Constants
    rE = 6378;
    omegaE = 2*pi/86400;
    
    theta = (stationNum-1) * pi/6;
    pos = rE* [cos(omegaE*t+theta);
        -omegaE*sin(omegaE*t+theta);
        sin(omegaE*t+theta);
        omegaE*cos(omegaE*t+theta)];
end