function plotLKFData(t, x, estErr, estErrCov, trackingData) 
    %% Plot nonlinear truth dynamics model simulation

    % Total states
    figure;
    subplot(4, 1, 1);
    plot(t, x(1, :));
    title('Noisy Truth States vs. Time, Typical LKF Run');
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
    
    %% Plot nonlinear model data simulation
    figure;
    for j = 1:length(trackingData)
        % Instantiate arrays for scatter plot for time and each measurement
        time = [];
        range = [];
        rangeRate = [];
        elevation = [];

        % Get tracking data for current station
        data = trackingData{j};
        for i=1:length(t)
            % Iterate through map to get all the times where a measurement
            % is available
            if isKey(data, t(i))
                trackPoint = data(t(i));
                time = [time; t(i)];
                range = [range; trackPoint(1,1)];
                rangeRate = [rangeRate; trackPoint(2,1)];
                elevation = [elevation; trackPoint(3,1)];
            end
        end

        % Plot with with single color for this tracking station
        subplot(3,1,1);
        hold on
        scatter(time, range);
        subplot(3,1,2);
        hold on
        scatter(time, rangeRate);
        subplot(3,1,3);
        hold on
        scatter(time, elevation);
    end

    % Add labels
    subplot(3,1,1);
    title('Noisy Simulated Measurements, Typical LKF Run');
    ylabel('Range (km)');
    subplot(3,1,2);
    ylabel('Range Rate (km/s)');
    subplot(3,1,3);
    ylabel('Elevation (rad)');
    xlabel('Time t (s)');
    hold off;

    %% Plot LKF estimation errors
    figure;
    threeSigma1 = [];
    sigma1High = [];
    sigma1Low = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma1 = [threeSigma1 3*sqrt(estErrCov(1,i*size(estErrCov, 1)-3))];
        sigma1High = [sigma1High estErr(1,i) + threeSigma1(i)];
        sigma1Low = [sigma1Low estErr(1,i) - threeSigma1(i)];
    end
    subplot(4,1,1);
    hold on;
    plot(t(2:end), estErr(1,:), 'b');
    plot(t(2:end), sigma1High, '--b');
    plot(t(2:end), sigma1Low, '--b');
    title('State Estimation Errors vs. Time, Typical LKF Run');
    ylabel('X Est. Error [km]');
    legend('Est. Error', '+/- 3\sigma bounds')
    hold off;

    threeSigma2 = [];
    sigma2High = [];
    sigma2Low = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma2 = [threeSigma2 3*sqrt(estErrCov(2,i*size(estErrCov, 1)-2))];
        sigma2High = [sigma2High estErr(2,i) + threeSigma2(i)];
        sigma2Low = [sigma2Low estErr(2,i) - threeSigma2(i)];
    end
    subplot(4,1,2);
    hold on;
    plot(t(2:end), estErr(2,:), 'b');
    plot(t(2:end), sigma2High, '--b');
    plot(t(2:end), sigma2Low, '--b');
    ylabel('Xdot Est. Error [km/s]');
    hold off;

    threeSigma3 = [];
    sigma3High = [];
    sigma3Low = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma3 = [threeSigma3 3*sqrt(estErrCov(3,i*size(estErrCov, 1)-1))];
        sigma3High = [sigma3High estErr(3,i) + threeSigma3(i)];
        sigma3Low = [sigma3Low estErr(3,i) - threeSigma3(i)];
    end
    subplot(4,1,3);
    hold on;
    plot(t(2:end), estErr(3,:), 'b');
    plot(t(2:end), sigma3High, '--b');
    plot(t(2:end), sigma3Low, '--b');
    ylabel('Y Est. Error [km]');
    hold off;

    threeSigma4 = [];
    sigma4High = [];
    sigma4Low = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma4 = [threeSigma4 3*sqrt(estErrCov(4,i*size(estErrCov, 1)))];
        sigma4High = [sigma4High estErr(4,i) + threeSigma4(i)];
        sigma4Low = [sigma4Low estErr(4,i) - threeSigma4(i)];
    end
    subplot(4,1,4);
    hold on;
    plot(t(2:end), estErr(4,:), 'b');
    plot(t(2:end), sigma4High, '--b');
    plot(t(2:end), sigma4Low, '--b');
    ylabel('Ydot Est. Error [km/s]');
    xlabel('Time t (s)');
    hold off;
end