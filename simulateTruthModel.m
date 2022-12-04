function [truthX, truthY] = simulateTruthModel(t, x0, Q, R, stationNums)
    % This function performs a complete nonlinear simulation of satellite
    % dynamics and tracking station measurements with process and
    % measurement noise. There is no forcing input.
    %
    % Inputs:
    % t = Matrix of time steps that you want to simulate at
    % x0 = 4x1 initial state matrix
    % Q = process noise covariance matrix
    % R = measurement noise covariance matrix
    % stationNums = matrix containing the IDs of the stations that are
    % generating measurements
    %
    % Outputs:
    % truthX = time-tagged truth satellite state vector
    % truthY = cell array of containing time-tagged maps of truth nonlinear
    % measurements from tracking stations
    %
    % truthY is a cell array containing maps for each tracking station. For
    % example, truthY[1] is the tracking data for station 1, and truthY[10]
    % is the tracking data for station 10. Each element is a map whose key
    % is the timestep and value is the 3x1 measurement matrix.

    % Constant
    mu = 398600;
    
    %% Compute nonlinear dynamics
    % Define nonlinear model as anonymous function and set options
    xdot = @(x, u, w) [x(2);
        -mu*x(1)/sqrt((x(1)^2+x(3)^2)^3) + u(1) + w(1);
        x(4);
        -mu*x(3)/sqrt((x(1)^2+x(3)^2)^3) + u(2) + w(2)];
    options = odeset('RelTol', 1e-12);

    % Initialize truth state vector and Cholesky decomp of Q matrix
    Sv = chol(Q, 'lower');
    truthX = x0;
    
    % Integrate orbit trajectory stepwise with process noise
    for i = 1:length(t)-1
        % Draw random noise
        w = Sv * randn(size(Sv, 1), 1);
        
        % Define wrapper function for use in ODE45
        xdotwrap = @(t, x) xdot(x, [0; 0], w);

        % Compute next step using ODE45
        [~, xout] = ode45(xdotwrap, [t(i) t(i+1)], truthX(:, end), options);

        % Store truth state output
        truthX = [truthX xout(end, :)'];
    end

    %% Compute nonlinear measurements
    % Initialize tracking data maps
    truthY = {};
    for i = 1:length(stationNums)
        truthY = [truthY; {containers.Map('KeyType', 'double','ValueType', 'any')}];
    end

    % Perform Cholesky decomp on R matrix
    Sv = chol(R, 'lower');

    for i = 1:length(t)
        for j = 1:length(stationNums)
            % Identify the current station's state vector
            stationID = stationNums(j);
            stationState = getStationState(t(i), stationID);
    
            % Compute satellite elevation angle and station angle to check
            % when satellite becomes visible
            satElevation = atan2(truthX(3, i)-stationState(3), truthX(1, i)-stationState(1));
            stationTheta = atan2(stationState(3), stationState(1));
    
            % Compute smallest angle between satellite and tracking station
            angleDiff = satElevation - stationTheta;
            angleDiff = abs(mod((angleDiff + pi), 2*pi) - pi);
    
            % Check to see if satellite is visible to station
            if angleDiff <= pi/2
                % Compute measurement with process noise
                rho = sqrt((truthX(1,i)-stationState(1))^2 + (truthX(3,i)-stationState(3))^2);
                outY = [rho;...
                    ((truthX(1,i)-stationState(1))*(truthX(2,i)-stationState(2)) + (truthX(3,i)-stationState(3))*(truthX(4,i)-stationState(4)))/rho;...
                    atan2(truthX(3, i)-stationState(3), truthX(1, i)-stationState(1))]...
                    + Sv*randn(size(R, 1), 1);
                
                % Query station map, add tracking data, and return to map array
                stationData = truthY{stationNums==stationID};
                stationData(t(i)) = outY;
                truthY{stationNums==stationID} = stationData;
            end
        end
    end
end