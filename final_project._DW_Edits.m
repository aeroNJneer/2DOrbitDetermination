close all
clear
rng(100)
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
omega = sqrt(mu/r0^3); % angular velocity of satellite
A0 = getAtilde([x0, x0dot, y0, y0dot]');
B = [0 0; 1 0; 0 0; 0 1]; 
G = dt * B; 
% end initial conditions

%% Nonlinear model simulation
pert = [0, 0.075, 0, -0.021]'; % perturbation 
xdot = @(x, u) [x(2);
    -mu*x(1)/sqrt((x(1)^2+x(3)^2)^3) + u(1);
    x(4);
    -mu*x(3)/sqrt((x(1)^2+x(3)^2)^3) + u(2)];
% Define initial total state
X0 = [x0, x0dot, y0, y0dot]' +  pert; 

% Define wrapper function for use in ODE45
xdotwrap = @(t, x) xdot(x, [0; 0]);
% Compute total state using ODE45
options = odeset('RelTol', 1e-12); % Set ODE45 options
[tout, xnl] = ode45(xdotwrap, t, X0, options);
xnl = xnl'; % Take transpose

% simulate measurements
Y = {};
for i = 1:length(stationNums)
    Y = [Y; {containers.Map('KeyType', 'double','ValueType', 'any')}];
end

Y = getMeasurements(Y, tout, xnl, stationNums, 0);
% plotNLData(tout, xnl, zeros(2), zeros(3), Y);
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

%% LKF Implementation
% Tuning inputs
numTrials = 1;
givenData = load('orbitdeterm_finalproj_KFdata.mat');
Qkf = givenData.Qtrue;
% Qkf = [-.1 0;0 -.01];
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

NEES = {};
NIS = {};

for i = 1:numTrials
    disp(['Running LKF trial ' int2str(i) '...']);
    % Initialize priors
    dxPlus = zeros(4,1);
    pPlus = pPlus0;    % Selected such that the satellite generally stays outside the Earth
    Sv = chol(pPlus, 'lower');

    % Draw intial state with perturbation and conduct truth model sim
    X0 = nomX0 + Sv * randn(size(Sv, 1), 1);
    [nomX, ~] = simulateTruthModel(t, nomX0, zeros(2), zeros(3), stationNums);
    [truthX, measY] = simulateTruthModel(t, X0, givenData.Qtrue, givenData.Rtrue, stationNums);
    
    % Initialize this run's outputs
    trialPerturbations = {};
    trialTotalStates = {};
    trialEstimationErrs = {};
    trialEstimationErrsCov = {};
    trialMeasInnovations = {};
    trialMeasInnovationsCov = {};
    trialNEES = {};
    trialNIS = {};
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
        trialNEES = [trialNEES (truthX(:, j)-trialTotalStates{:, end})' * (pPlus)^(-1) * (truthX(:, j)-trialTotalStates{:, end})];
        trialNIS = [trialNIS dy'*(S)^(-1)*dy];
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
    NEES = [NEES;trialNEES];
    NIS = [NIS;trialNIS];
    
end



% NEES_bar = mean(cell2mat(NEES),1);
% NIS_bar = mean(cell2mat(NIS),1);
% alpha = .05;
% r1 = chi2inv(alpha/2,4*numTrials)/numTrials;
% r2 = chi2inv(1-alpha/2,4*numTrials)/numTrials;
% 
% figure
% subplot(2,1,1)
% scatter(t(2:end),NEES_bar)
% xlabel('Time (s)')
% ylabel('NEES statistic')
% 
% subplot(2,1,2)
% scatter(t(2:end),NIS_bar)
% xlabel('Time (s)')
% ylabel('NIS statistic')

% Plot typical simulation from first sim run
plotLKFData(t, truthTotalStates{1, :}, cell2mat(estimationErrs(1, :)), cell2mat(estimationErrsCov(1,:)), simulatedMeasurements(1, :));

% EKF Simulation