function plotLinearizedData(t, x, dx, xnl, trackingData)    
    figure(1);
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
    
    figure(2);
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
    
    
    figure(3);
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
    
    figure(4);
    for j = 1:length(trackingData)
        time = [];
        range = [];
        rangeRate = [];
        elevation = [];
        data = trackingData{j};
        for i=1:length(t)
            if isKey(data, t(i))
                trackPoint = data(t(i));
                time = [time; t(i)];
                range = [range; trackPoint(1,1)];
                rangeRate = [rangeRate; trackPoint(2,1)];
                elevation = [elevation; trackPoint(3,1)];
            end
        end
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
    subplot(3,1,1);
    title('Approximate Linearized Model Data Simulation');
    ylabel('Range (km)');
    subplot(3,1,2);
    ylabel('Range Rate (km/s)');
    subplot(3,1,3);
    ylabel('Elevation (rad)');
    xlabel('Time t (s)');
    hold off;
end