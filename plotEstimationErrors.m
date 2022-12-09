function plotEstimationErrors(t, estErr, estErrCov, title)
   %% Plot LKF estimation errors
    figure;
    threeSigma1 = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma1 = [threeSigma1 3*sqrt(estErrCov(1,i*size(estErrCov, 1)-3))];
    end
    subplot(4,1,1);
    hold on;
    plot(t(2:end), estErr(1,:), 'b');
    plot(t(2:end), threeSigma1, '--b');
    plot(t(2:end), -threeSigma1, '--b');
    sgtitle(title);
    ylabel('X Est. Error [km]');
    legend('Est. Error', '+/- 3\sigma bounds')
    hold off;

    threeSigma2 = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma2 = [threeSigma2 3*sqrt(estErrCov(2,i*size(estErrCov, 1)-2))];
    end
    subplot(4,1,2);
    hold on;
    plot(t(2:end), estErr(2,:), 'b');
    plot(t(2:end), threeSigma2, '--b');
    plot(t(2:end), -threeSigma2, '--b');
    ylabel('Xdot Est. Error [km/s]');
    hold off;

    threeSigma3 = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma3 = [threeSigma3 3*sqrt(estErrCov(3,i*size(estErrCov, 1)-1))];
    end
    subplot(4,1,3);
    hold on;
    plot(t(2:end), estErr(3,:), 'b');
    plot(t(2:end), threeSigma3, '--b');
    plot(t(2:end), -threeSigma3, '--b');
    ylabel('Y Est. Error [km]');
    hold off;

    threeSigma4 = [];
    for i = 1:size(estErrCov, 2)/size(estErrCov, 1)
        threeSigma4 = [threeSigma4 3*sqrt(estErrCov(4,i*size(estErrCov, 1)))];
    end
    subplot(4,1,4);
    hold on;
    plot(t(2:end), estErr(4,:), 'b');
    plot(t(2:end), threeSigma4, '--b');
    plot(t(2:end), -threeSigma4, '--b');
    ylabel('Ydot Est. Error [km/s]');
    xlabel('Time t (s)');
    hold off;

end