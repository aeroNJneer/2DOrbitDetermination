function plotMeasurements(t, trackingData, title)
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
    sgtitle(title);
    subplot(3,1,1);
        ylabel('\rho Range (km)');
    subplot(3,1,2);
        ylabel('\rhodot Range Rate (km/s)');
    subplot(3,1,3);
        ylabel('\phi Elevation (rad)');
        xlabel('Time t (s)');
    hold off;
end
        