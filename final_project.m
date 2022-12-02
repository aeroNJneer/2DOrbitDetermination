close all
clear
mu = 398600; % gr62)avitational constant
omegaE = 2*pi/86400;  % angular velocity of earth
rE = 6378; % radius of earth
dt = 10;
t = 0:dt:14000;
% initial conditions
x0 = 6678;
y0 = 0;
r0 = sqrt(x0^2+y0^2);
x0dot = 0;
y0dot = r0 * sqrt(mu/r0^3);
%period = 2 * pi * sqrt(r0^3/mu); % orbital period
omega = sqrt(mu/r0^3); % angular velocity of satellite
% Atilde is time varying so need to create paramaterized variable.... 
% -- pass in r and t (r0 for nominal)
% Atilde = [0 1 0 0; 
%     mu*(2*(cos(omega*t))^2-(sin(omega*t))^2)/r^3 0 3*mu*sin(omega*t)*cos(omega*t)/r^3 0;
%     0 0 0 1;
%     3*mu*sin(omega*t)*cos(omega*t)/r^3 0 mu*(2*(sin(omega*t)^2 -(cos(omega*t))^2))/r^3 0];
A0 = [0 1 0 0 ; 2*mu/r0^3 0 0 0 ; 0 0 0 1; 0 0 -mu/r0^3 0];
B = [0 0; 1 0; 0 0; 0 1]; 
G = dt * B; 

% Define nonlinear model as anonymous function
pert = [0, 0.075, 0, -0.021]'; % perturbation 
xdot = @(x, u) [x(2);
    -mu*x(1)/sqrt((x(1)^2+x(3)^2)^3) + u(1);
    x(4);
    -mu*x(3)/sqrt((x(1)^2+x(3)^2)^3) + u(2)];
% Define initial total state
X0 = [x0, x0dot, y0, y0dot]' +  pert; 
% 
% Define wrapper function for use in ODE45
xdotwrap = @(t, x) xdot(x, [0; 0]);
% Compute total state using ODE45
options = odeset('RelTol', 1e-12); % Set ODE45 options
[tout, xnl] = ode45(xdotwrap, t, X0, options);
xnl = xnl'; % Take transpose so that dimensions match with xnom

% Plot
figure;
subplot(4, 1, 1);
plot(t, xnl(1, :));
title('Total States vs. Time, Nonlinear Dynanmics - ODE45');
ylabel('X (km)');
grid on;
subplot(4, 1, 2);
plot(t, xnl(2, :));
ylabel('Xdot (km/s)');
grid on;
subplot(4, 1, 3);
plot(t, xnl(3, :));
ylabel('Y (km)');
grid on;
subplot(4, 1, 4);
plot(t, xnl(4, :));
ylabel('Ydot (km/s)');
xlabel('Time t (s)');
grid on;
% end of NL simulation

% % linearized nominal simulation
% Define empty perturbation state and measurement vectors as well as zero control input
x = zeros(4, length(t));
dx = zeros(4, length(t));
u = zeros(2, length(t));

% Compute nominal state as a function of true anomaly
theta = omega * t;
nomState = [r0*cos(theta); -r0*omega*sin(theta); r0*sin(theta); r0*omega*cos(theta)];

% Initialize perturbation and total state
dx(:, 1) = [0, 0.075, 0, -0.021]';
x(:, 1) = nomState(:,1) + dx(:,1);

for i = 2:length(t)
    % Compute perturbation states
    F = eye(4) + getAtilde(nomState(:, i))*dt;
    dx(:, i) = F*dx(:, i-1) + G*u(:, i-1);

    % Compute total states
    x(:, i) = nomState(:, i) + dx(:,i);

    % Compute measurement perturbations for all stations
    curStepMeas = {};
    for j = 1:12
        stationState = getStationState(t, j);

        % Check to see if satellite is visible to station
        satElevation = atan((x(3, i)-stationState(3))/(x(1, i)-stationState(1)));
        stationTheta = atan(stationState(3)/stationState(1));
        if satElevation>=-pi/2+stationTheta && satElevation<=pi/2+stationTheta
            % Compute a measurement for station j at time t
            
        end
    end
end

% Plot
figure;
subplot(4, 1, 1);
plot(t, x(1, :));
title('Total States vs. Time, Linearized DT SS Model');
ylabel('X (km)'); ylim([-10000, 10000])
grid on;
subplot(4, 1, 2);
plot(t, x(2, :));
ylabel('Xdot (km/s)'); ylim([-10, 10])
grid on;
subplot(4, 1, 3);
plot(t, x(3, :));
ylabel('Y (km)'); ylim([-10000, 10000])
grid on;
subplot(4, 1, 4);
plot(t, x(4, :));
ylabel('Ydot (km/s)'); ylim([-10, 10])
xlabel('Time t (s)');
grid on;

figure;
subplot(4,1,1);
plot(t,dx(1,:));
title('Approx Perturbations vs. Time, Linearized DT SS Model');
ylabel('\deltaX (km)'); ylim([-500, 500])
grid on
subplot(4,1,2);
plot(t,dx(2,:));
ylabel('\deltaXdot (km/s)'); ylim([-1, 1])
grid on
subplot(4,1,3);
plot(t,dx(3,:));
ylabel('\deltaY (km)'); ylim([-1000, 1000])
grid on
subplot(4,1,4);
plot(t,dx(4,:));
ylabel('\deltaYdot (km/s)'); ylim([-1, 1])
xlabel('Time t (s)');
grid on


figure;
subplot(4, 1, 1);
plot(t, x(1, :) - xnl(1, 1:length(t)));
title('State Error vs. Time, Linearized DT SS Model vs. ODE45');
ylabel('x (km)');
grid on;
subplot(4, 1, 2);
plot(t, x(2, :) - xnl(2, 1:length(t)));
ylabel('dx (km/s)');
grid on;
subplot(4, 1, 3);
plot(t, x(3, :) - xnl(3, 1:length(t)));
ylabel('y (km)');
grid on;
subplot(4, 1, 4);
plot(t, x(4, :) - xnl(4, 1:length(t)));
ylabel('dy (km/s)');
xlabel('Time t (s)');
grid on


