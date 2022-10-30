% [1/s] friction damping coefficient
k = [0.1, 0.5, 0.75, 1, 5, 10, 50, 100]; % logspace(-2,3,10);
l2_iter = k*l1; % [m]
vars_iter = zeros(size(vars,1), size(vars,2),length(k));
A0_iter = zeros(6,6,length(k));
B0_iter = zeros(6,1,length(k));
gains = zeros(1,6,length(k));
Abar = zeros(size(A0_iter));
t = 0:0.01:4;
states = {'x' 'theta_1' 'theta_1' 'x_dot' 'theta_dot_1' 'theta_dot_2'};
inputs = {'u'};
outputs = {'x'; 'theta_1'; 'theta_2'};
clear clsys_iter
for i = 1:length(k)
    vars_iter(:,:,i) = [d0_real, d1_real, d2_real, g_real, l1_real, l2_iter(i), m1_real, m2_real, m_q_real];
    A0_iter(:,:,i) = double(subs(A0,vars,vars_iter(:,:,i)));
    B0_iter(:,:,i) = double(subs(B0,vars,vars_iter(:,:,i)));
    % placing poles at the same locations... will see if that has an effect
    gains(:,:,i) = acker(A0_iter(:,:,i),B0_iter(:,:,i),p);
    % u = -Kx + r => xdot = (A - B*K)*x + B*r
    Abar(:,:,i) = A0_iter(:,:,i)-B0_iter(:,:,i)*gains(:,:,i);
    % Closed Loop system
    clsys_iter(:,:,i) = ss(Abar(:,:,i),B0_iter(:,:,i),C,D,'statename',states,'inputname',inputs,'outputname',outputs);
    [y(:,:,i), ~, x_iter(:,:,i)] = initial(clsys_iter(:,:,i),[0.05; -10*pi/180; -2*pi/180; 0; 0; 0],t);
end

%%

figure(20);
tiledlayout(3,1);
ax1 = nexttile; hold on; grid on;
ax1.ColorOrder = [1 0 0; 0 0 1; 0 0.7 0; 0 0 0; 0.929 0.694 0.125];
ax1.LineStyleOrder = {'-','--'};
ax2 = nexttile; hold on; grid on;
ax2.ColorOrder = [1 0 0; 0 0 1; 0 0.7 0; 0 0 0; 0.929 0.694 0.125];
ax2.LineStyleOrder = {'-','--'};
ax3 = nexttile; hold on; grid on;
ax3.ColorOrder = [1 0 0; 0 0 1; 0 0.7 0; 0 0 0; 0.929 0.694 0.125];
ax3.LineStyleOrder = {'-','--'};
title(ax1,'Response to Initial Conditions');
ylabel(ax1,'To: x [m]');
ylabel(ax2,'To: \theta_1 [deg]');
ylabel(ax3,'To: \theta_2 [deg]'); xlabel('Time (seconds)');
for i = 1:length(k)
    plot(ax1,t, y(:,1,i)); 
    
    plot(ax2,t,y(:,2,i)*180/pi); 

    plot(ax3,t,y(:,3,i)*180/pi); 

end
maxmin = 5;
degmax = 30;
ylim(ax1,[-1,1]*maxmin)
ylim(ax2,[-1 1]*degmax)
ylim(ax3,[-1 1]*degmax)
axes(ax1)
legend(num2str(k'))