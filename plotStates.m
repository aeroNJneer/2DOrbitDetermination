function plotStates(t, x, title)
    figure;
    sgtitle(title);
    subplot(4, 1, 1);
        plot(t, x(1, :));
        ylabel('X (km)');
        grid on;
    subplot(4, 1, 2);
        plot(t, x(2, :));
        ylabel('\deltaX (km/s)');
        grid on;
    subplot(4, 1, 3);
        plot(t, x(3, :));
        ylabel('Y (km)');
        grid on;
    subplot(4, 1, 4);
        plot(t, x(4, :));
        ylabel('\deltaY (km/s)');
        xlabel('Time t (s)');
    grid on;
end
