function t = myplot2(tOut,x)
    t = tiledlayout(4,1);
    title(t, {"$\bf Response$"; "From: u"}, 'Interpreter', 'Latex')
    ylabel(t, 'Amplitude')
    nexttile;
    stairs(tOut,x(:,1))
    ylabel('to: x')
    nexttile;
    stairs(tOut,x(:,6)) % xhat = x - xtilde
    ylabel('to: $x_{err}$', 'Interpreter','latex')
    nexttile;
    stairs(tOut,x(:,3)*180/pi)
    ylabel('to: \theta [deg]')
    nexttile;
    stairs(tOut,x(:,8)*180/pi)
    ylabel('to: $\theta_{err}$ [deg]', 'Interpreter','latex')

end