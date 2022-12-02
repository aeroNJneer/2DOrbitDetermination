function pos = getStationPos(t, stationNum)
    % Constants
    rE = 6378;
    omegaE = 2*pi/86400;
    
    theta = (stationNum-1) * pi/6;
    pos = rE* [cos(omegaE*t+theta); sin(omegaE*t+theta)];
end