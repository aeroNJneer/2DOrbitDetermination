function plotNLData(t,xnl, Q, R, trackingData)
% Plot states
    figure;
    sgtitle('Total States vs. Time, Nonlinear Dynanmics - ODE45');
    subplot(4, 1, 1);
        plot(t, xnl(1, :));
        ylabel('X (km)');
        grid on;
    subplot(4, 1, 2);
        plot(t, xnl(2, :));
        ylabel('\deltax (km/s)');
        grid on;
    subplot(4, 1, 3);
        plot(t, xnl(3, :));
        ylabel('Y (km)');
        grid on;
    subplot(4, 1, 4);
        plot(t, xnl(4, :));
        ylabel('\deltay (km/s)');
        xlabel('Time t (s)');
    grid on;
% Plot measurements
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
    sgtitle('Nonlinear Model Data Simulation');
    subplot(3,1,1);
        ylabel('Range (km)');
    subplot(3,1,2);
        ylabel('Range Rate (km/s)');
    subplot(3,1,3);
        ylabel('Elevation (rad)');
        xlabel('Time t (s)');
    hold off;
end
        