function Htilde = getHtilde(xNom, stationPos)
    Htilde = [(xNom(1)-stationPos(1))/((xNom(1)-stationPos(1))^2 + (xNom(3)-stationPos(3))^2)^0.5 0 (xNom(3)-stationPos(3))/((xNom(2)-stationPos(1))^2 + (xNom(3)-stationPos(3))^2)^0.5 0;...
        ]; % Not yet complete, only first row is here
end