close all
clear
mu = 398600; % gravitational constant
omegaE = 2*pi/86400;  % angular velocity of earth
rE = 6378; % radius of earth
dt = 10;
t = 0:dt:14000;
stationNums = 1:12;
% initial conditions
x0 = 6678;
y0 = 0;
r0 = sqrt(x0^2+y0^2);
x0dot = 0;
y0dot = r0 * sqrt(mu/r0^3);
X0 = [x0, x0dot, y0, y0dot]';
omega = sqrt(mu/r0^3);          % angular velocity of satellite
A0 = getAtilde(X0);
B = [0 0; 1 0; 0 0; 0 1]; 
G = dt * B; 
% end initial conditions

%% Nonlinear model simulation
pert = [0, 0.075, 0, -0.021]';  % perturbation 
X0 = X0 +  pert;                % Define initial total state
[tout, xnl] = nlmodel(t, X0, mu, [0 0]);
xnl = xnl';                     % Take transpose
% simulate measurements
Y = {};
for i = 1:length(stationNums)
    Y = [Y; {containers.Map('KeyType', 'double','ValueType', 'any')}];
end

Y = getMeasurements(Y, tout, xnl, stationNums, 0);
plotStates(tout, xnl, "Total States vs. Time, Nonlinear Dynanmics - ODE45");
plotMeasurements(tout, Y, "Nonlinear Model Data Simulation");
% end of NL simulation

%% Approximate linearized simulation
% Define empty perturbation state and measurement vectors as well as zero control input
x = zeros(4, length(t));
dx = zeros(4, length(t));
u = zeros(2, length(t));

% Compute nominal trajectory as a function of true anomaly
trueAnom = omega * t;
nomX0 = [r0*cos(trueAnom); -r0*omega*sin(trueAnom); r0*sin(trueAnom); r0*omega*cos(trueAnom)];

% Initialize perturbation and total state
dx(:, 1) = [0, 0.075, 0, -0.021]';
x(:, 1) = nomX0(:,1) + dx(:,1);

% Initialize array of tracking data maps and specify which stations I want to use
trackingData = {};
for i = 1:length(stationNums)
    trackingData = [trackingData; {containers.Map('KeyType', 'double','ValueType', 'any')}];
end

% Run simulation
for i = 2:length(t)
    % Compute perturbation states
    F = eye(4) + getAtilde(nomX0(:, i-1))*dt;
    dx(:, i) = F*dx(:, i-1) + G*u(:, i-1);

    % Compute total states
    x(:, i) = nomX0(:, i) + dx(:,i);

    % Compute measurement perturbations for each station
    for j = 1:length(stationNums)
        % Identify the current station's state vector
        stationID = stationNums(j);
        stationState = getStationState(t(i), stationID);

        % Compute satellite elevation angle and station angle using truth
        % (ODE45) data when satellite actually becomes visible
        satElevation = atan2(xnl(3, i)-stationState(3), xnl(1, i)-stationState(1));
        stationTheta = atan2(stationState(3), stationState(1));

        % Compute angle between satellite and tracking station
        angleDiff = satElevation - stationTheta;
        angleDiff = abs(mod((angleDiff + pi), 2*pi) - pi);

        % Check to see if satellite is visible to station
        if angleDiff <= pi/2
            % Compute nominal measurement
            nomMeas = getNLMeasurement(t(i), nomX0(:, i), stationID);

            % Compute a measurement for station j at time t(i)
            Htilde = getHtilde(nomX0(:, i), stationState);
            dy = Htilde*dx(:, i);
            
            % Query station map, add tracking data, and return to map array
            stationData = trackingData{find(stationNums==stationID)};
            % Data is measurement in first column and Htilde is
            % right-appended.
            stationData(t(i)) = [nomMeas+dy Htilde];
            trackingData{find(stationNums==stationID)} = stationData;
        end
    end
end

% Plot
% plotLinearizedData(t,x,dx,xnl,trackingData);
plotStates(t,x, "Total States vs. Time, Linearized DT SS Model");
plotStates(t,dx, "Approx Perturbations vs. Time, Linearized DT SS Model");
plotStates(t,x-xnl, "State Error vs. Time, Linearized DT SS Model vs. ODE45");
plotMeasurements(t,trackingData, "Approximate Linearized Model Data Simulation");
%% LKF Implementation
% Tuning inputs
numTrials = 1;
givenData = load('orbitdeterm_finalproj_KFdata.mat');
Qkf = givenData.Qtrue;
pPlus0 = diag([.1 .001 .1 .001]);    % Selected such that the satellite doesn't crash into Earth

% Setup
stationNums = 1:12;
dt = 10;
t = 0:dt:14000;
r0 = 6678;
mu = 398600;
nomX0 = [r0; 0; 0; sqrt(mu/r0)];
gamma = [0 0; 1 0; 0 0; 0 1]; 
omegaTilde = dt * gamma; 

% Initialize outputs. Each row is a different trial
filteredPerturbations = {};
filteredTotalStates = {};
truthTotalStates = {};
simulatedMeasurements = {};
estimationErrs = {};
estimationErrsCov = {};
measInnovations = {};
measInnovationsCov = {};

for i = 1:numTrials
    disp(['Running LKF trial ' int2str(i) '...']);
    % Initialize priors
    dxPlus = zeros(4,1);
    pPlus = pPlus0;    % Selected such that the satellite generally stays outside the Earth
    Sv = chol(pPlus, 'lower');

    % Draw intial state with perturbation and conduct truth model sim
    X0 = nomX0 + Sv * randn(size(Sv, 1), 1);
    [nomX, nomY] = simulateTruthModel(t, nomX0, zeros(2), zeros(3), stationNums);
    [truthX, measY] = simulateTruthModel(t, X0, givenData.Qtrue, givenData.Rtrue, stationNums);
    
    % Initialize this run's outputs
    trialPerturbations = {};
    trialTotalStates = {};
    trialEstimationErrs = {};
    trialEstimationErrsCov = {};
    trialMeasInnovations = {};
    trialMeasInnovationsCov = {};

    % Begin LKF filtering
    for j = 2:length(t)
        % Update step
        Ftilde = eye(4) + getAtilde(nomX(:, j-1))*dt;
        dxMinus = Ftilde*dxPlus;
        pMinus = Ftilde*pPlus*Ftilde' + omegaTilde*Qkf*omegaTilde';
        
        % Get stacked station measurements
        [y, stationIDs] = getStackedMeasAtTime(t(j), measY);
        yStar = [];
        HtildeStacked = [];
        RStacked = [];

        if ~isempty(y)
            % If there are measurements, then perform correction step by
            % stacking nominal measurements y* and getting corresponding
            % stacked Htilde and R matrices
            for k=1:length(stationIDs)
                yStar = [yStar; getNLMeasurement(t(j), nomX(:, j), stationIDs(k))];
                HtildeStacked = [HtildeStacked; getHtilde(nomX(:, j), getStationState(t(j), stationIDs(k)))];
                RSize = size(givenData.Rtrue, 1);   % R is square
                RStacked = [RStacked [zeros((k-1)*RSize, RSize); givenData.Rtrue; zeros((length(stationIDs)-k)*RSize, RSize)]];
            end
    
            % Compute Kalman gain
            K = pMinus * HtildeStacked' / (HtildeStacked*pMinus*HtildeStacked' + RStacked);
    
            % Update estimation error covariance and compute measurement
            % innovation as well as new perturbation state
            pPlus = (eye(size(K, 1))-K*HtildeStacked)*pMinus;
            dy = y-yStar;
            dxPlus = dxMinus + K * (y-yStar - HtildeStacked*dxMinus);

            % Compute measurement innovation covariance matrix
            S = HtildeStacked*pMinus*HtildeStacked' + RStacked;
        else
            % No measurements available, so carry forward pMinus and state
            pPlus = pMinus;
            dxPlus = dxMinus;

            % dy and S are null
            dy = [nan; nan; nan];
            S = [nan nan nan; nan nan nan; nan nan nan];
        end
        
        % Store outputs
        trialPerturbations = [trialPerturbations dxPlus];
        trialTotalStates = [trialTotalStates nomX(:, j)+dxPlus];
        trialEstimationErrs = [trialEstimationErrs truthX(:, j)-trialTotalStates{:, end}];
        trialEstimationErrsCov = [trialEstimationErrsCov pPlus];
        trialMeasInnovations = [trialMeasInnovations dy];
        trialMeasInnovationsCov = [trialMeasInnovationsCov S];
    end

    % Store trial output
    filteredPerturbations = [filteredPerturbations; trialPerturbations];
    filteredTotalStates = [filteredTotalStates; trialTotalStates];
    truthTotalStates = [truthTotalStates; truthX];
    simulatedMeasurements = [simulatedMeasurements; measY'];
    estimationErrs = [estimationErrs; trialEstimationErrs];
    estimationErrsCov = [estimationErrsCov; trialEstimationErrsCov];
    measInnovations = [measInnovations; trialMeasInnovations];
    measInnovationsCov = [measInnovationsCov; trialMeasInnovationsCov];
end

% Plot Monte Carlo simulated states and measurements
plotStates(t, nomX, "Truth States vs. Time, Monte Carlo Run");
plotStates(t, truthX, "Noisy Truth States vs. Time, Monte Carlo Run");
plotMeasurements(t, measY, "Noisy Simulated Measurements, Monte Carlo Run")
%plotSimulationData(t, nomX, truthXnom, nomY, measY, stationNums);

% Plot typical simulation from first sim run
plotStates(t, truthTotalStates{1, :}, "Noisy Truth States vs. Time, Typical LKF Run");
%plotStates(t, filteredTotalStates{1, :}, "Noisy Simulated States vs. Time, Typical LKF Run");

plotMeasurements(t, simulatedMeasurements(1, :), "Noisy Simulated Measurements, Typical LKF Run");
plotEstimationErrors(t, cell2mat(estimationErrs(1, :)), cell2mat(estimationErrsCov(1,:)) , "State Estimation Errors vs. Time, Typical LKF Run")
%plotLKFData(t, truthTotalStates{1, :}, cell2mat(estimationErrs(1, :)), cell2mat(estimationErrsCov(1,:)), simulatedMeasurements(1, :));


%% EKF Implementation
% Tuning inputs
% stationNums = 1:12;
% t = 0:dt:14000;
clear
close all
numTrials = 2;
givenData = load('orbitdeterm_finalproj_KFdata.mat');
Qkf = givenData.Qtrue;
pPlus0 = diag([.1 .001 .1 .001]);    % Selected such that the satellite doesn't crash into Earth

% Setup
rand(100);
stationNums = 1:12;
dt = 10;
t = 0:dt:14000;
r0 = 6678;
mu = 398600;
nomX0 = [r0; 0; 0; sqrt(mu/r0)];    % initial conditions w/no noise
gamma = [0 0; 1 0; 0 0; 0 1]; 
omegaTilde = dt * gamma; 

% Initialize outputs. Each row is a different trial
filteredPerturbations = {};
filteredTotalStates = {};
truthTotalStates = {};
simulatedMeasurements = {};
estimationErrs = {};
estimationErrsCov = {};
measInnovations = {};
measInnovationsCov = {};

for i = 1:numTrials
    disp(['Running EKF trial ' int2str(i) '...']);
    % step 1: Initialize priors
    xhatPlus = zeros(4,1);
    pPlus = pPlus0;    % Selected such that the satellite generally stays outside the Earth
    Sv = chol(pPlus, 'lower');
    
    % Draw intial state with perturbation and conduct truth model sim
    X0 = nomX0 + Sv * randn(size(Sv, 1), 1);

    % simulate with no noise
    [xnominal, ~] = simulateTruthModel(t, nomX0, zeros(2), zeros(3), stationNums);
    % Monte Carlo simulation for measurements
    [xnoisy, measY] = simulateTruthModel(t, X0, givenData.Qtrue, givenData.Rtrue, stationNums);
   
    % compute xhatPlus with no noise at nominal point
    [xhatPlus, ~] = simulateTruthModel([t(1) t(2)], nomX0, zeros(2), zeros(3), stationNums);
    
    xhatPlus = xhatPlus(:,end);
    % Initialize this run's outputs
    trialPerturbations = {};
    trialTotalStates = {};
    trialEstimationErrs = {};
    trialEstimationErrsCov = {};
    trialMeasInnovations = {};
    trialMeasInnovationsCov = {};

    % Begin EKF filtering
    for j = 2:length(t)
        % Update step
        Ftilde = eye(4) + getAtilde(xhatPlus)*dt;
        pMinus = Ftilde*pPlus*Ftilde' + omegaTilde*Qkf*omegaTilde';
        xhatMinus = Ftilde * xhatPlus;

        % Get stacked station measurements
        [y, stationIDs] = getStackedMeasAtTime(t(j), measY);
        dy = [];
        yhat = [];
        HtildeStacked = [];
        RStacked = [];
        if ~isempty(y)
            % If there are measurements, then perform correction step by
            % stacking nominal measurements y* per station and getting corresponding
            % stacked Htilde and R matrices
            for k=1:length(stationIDs)
                % Compute a measurement for station j at time t(i)
                yhat = [yhat; getNLMeasurement(t(j), xhatMinus, stationIDs(k))];
                HtildeStacked = [HtildeStacked; getHtilde(xhatMinus, getStationState(t(j), stationIDs(k)))];
                RSize = size(givenData.Rtrue, 1);   % R is square
                RStacked = [RStacked [zeros((k-1)*RSize, RSize); givenData.Rtrue; zeros((length(stationIDs)-k)*RSize, RSize)]];
            end

            % Compute Kalman gain
            K = pMinus * HtildeStacked' / (HtildeStacked*pMinus*HtildeStacked' + RStacked);
    
            % Update estimation error covariance and compute measurement innovation 
            pPlus = (eye(size(K, 1))-K*HtildeStacked)*pMinus;
            dy = y-yhat;
            xhatPlus = xhatMinus + (K * y);

            % Compute measurement innovation covariance matrix
            S = HtildeStacked*pMinus*HtildeStacked' + RStacked;
        else
            % No measurements available, so carry forward pMinus and state
            pPlus = pMinus;
            xhatPlus = xhatMinus;
            % dy and S are null
            dy = [nan; nan; nan];
            S = [nan nan nan; nan nan nan; nan nan nan];
        end
        
        % Store outputs
        trialTotalStates = [trialTotalStates xhatPlus];
        trialEstimationErrs = [trialEstimationErrs xnominal(:, j)-xhatPlus];
        trialEstimationErrsCov = [trialEstimationErrsCov pPlus];
        trialMeasInnovations = [trialMeasInnovations dy];
        trialMeasInnovationsCov = [trialMeasInnovationsCov S];
    end

    % Store trial output
    filteredPerturbations = [filteredPerturbations; trialPerturbations];
    filteredTotalStates = [filteredTotalStates; trialTotalStates];
    truthTotalStates = [truthTotalStates; xnominal];
    simulatedMeasurements = [simulatedMeasurements; yhat];
    estimationErrs = [estimationErrs; trialEstimationErrs];
    estimationErrsCov = [estimationErrsCov; trialEstimationErrsCov];
    measInnovations = [measInnovations; trialMeasInnovations];
    measInnovationsCov = [measInnovationsCov; trialMeasInnovationsCov];
end

% Plot typical simulation from first sim run
plotStates(t, xnoisy , "Noisy Truth States vs. Time, Typical EKF Run");
%plotStates(t, trialTotalStates{1, :} , "Noisy Truth States vs. Time, Typical EKF Run");
plotEstimationErrors(t, cell2mat(estimationErrs(1, :)), cell2mat(estimationErrsCov(1,:)), "State Estimation Errors vs. Time, Typical EKF Run" );


