%% Final Project Scratchwork
% Summary of example objective

%% Modeling Pendulum, no friction

m = 1; % kg
l = .1; % m
a = 9.8; % m/s^2
% EOM: m*l^2*theta_ddot = l*m*a*sin(theta)
tspan = linspace(0,10,10^4);
theta0 = [pi/2-pi/6 0];
[t,theta] = ode45(@(t,theta) pend_nofriction(t,theta,a,m,l),tspan,theta0);
figure;
plot(t,theta)

% note: fanimator could work, but have to figure out how to get the time
% noted
% x_pos = sin(thetaSolPlot);
% y_pos = -cos(thetaSolPlot);
% fanimator(@fplot,x_pos,y_pos,'ko','MarkerFaceColor','k','AnimationRange',[0 5*T]);
% hold on;
% fanimator(@(t) plot([0 x_pos(t)],[0 y_pos(t)],'k-'),'AnimationRange',[0 5*T]);
% fanimator(@(t) text(-0.3,0.3,"Timer: "+num2str(t,2)+" s"),'AnimationRange',[0 5*T]);


%% Animation, No Friction

th = linspace(0,2*pi) ; % the points for generating the plot of the circle mass
r = 0.01; % radius of circle
cx = l*cos(theta(:,1)); % x-axis location of the center of the circle
cy = l*sin(theta(:,1)); % y-axis location of the center of the circle
x = cx(1) + r*cos(th) ; 
y = cy(1) + r*sin(th) ;
bar = [0 0; cx(1) cy(1)];
figure;
h = plot(x,y,'r') ;
hold on
k = plot(bar(:,1),bar(:,2),'b');
axis equal 
%%animation
for i = 2:length(cx)
    bar = [0 0; cx(i) cy(i)];
    x = cx(i)+r*cos(th) ;
    y = cy(i)+r*sin(th) ;
    set(h,'XData',x,'YData',y) ;
    set(k,'XData',bar(:,1),'YData',bar(:,2));
    axis([-l l -l l]*2)    
    drawnow
    pause(0.0002)
end

%% Modeling Pendulum, friction

m = 1; % kg
l = .1; % m
a = 9.8; % m/s^2
mu = 0.05; % friction factor (can also capture all losses in drive train)
% EOM: m*l^2*theta_ddot = l*m*a*sin(theta)
tspan = linspace(0,10,10^3);
theta0 = [pi/2-pi/6 0];
[t,theta] = ode45(@(t,theta) pend_friction(t,theta,a,m,l,mu),tspan,theta0);
figure;
plot(t,rad2deg(theta(:,1)),t,theta(:,2))

%% Animation, Friction

th = linspace(0,2*pi) ; % the points for generating the plot of the circle mass
r = 0.01; % radius of circle
cx = l*cos(theta(:,1)); % x-axis location of the center of the circle
cy = l*sin(theta(:,1)); % y-axis location of the center of the circle
x = cx(1) + r*cos(th) ; 
y = cy(1) + r*sin(th) ;
bar = [0 0; cx(1) cy(1)];
figure;
h = plot(x,y,'r') ;
hold on
k = plot(bar(:,1),bar(:,2),'b');
axis equal 
%%animation
for i = 2:length(cx)
    bar = [0 0; cx(i) cy(i)];
    x = cx(i)+r*cos(th) ;
    y = cy(i)+r*sin(th) ;
    set(h,'XData',x,'YData',y) ;
    set(k,'XData',bar(:,1),'YData',bar(:,2));
    axis([-l l -l l]*2)    
    drawnow
    pause(0.0005)
end

%% Linearization and State Space Feedback
% now adding in input to the system

% Linearized model:
m = 1; % kg
l = .1; % m
a = 9.8; % m/s^2
mu = 0.05; % friction factor (can also capture all losses in drive train)
A = [0 1; a/l 0];
B = [0; 1/(m*l^2)];
C = [1 0];

% State Feedback, placing poles at -1 and -2
poles = [-1; -2];
K = acker(A,B,poles) % gains of controller u = -K*xtilde + r
Abar = A-B*K; % closed-loop A matrix
sysgain = -C*inv(A-B*K)*B % system gain
figure;
olsys = ss(A,B,C,0);
step(olsys)
figure;
clsys = ss(Abar,B,C,0);
step(clsys)

yd = pi/2; % desired output of theta is at 90deg, or pi/2
r = -inv(C*inv(Abar)*B)*yd;
figure;
% step response with input so that output is at desired location
step(clsys*r) 

% simulating linearized CL system with initial conditions theta0
t = linspace(0,10,10^3)';
u = repmat(r,[length(t),1]);
theta0 = [pi/2-pi/4; pi];
[ytilde,t,theta] = lsim(clsys,u,t,theta0);
figure;
plot(t,theta)

%% Animation, Linearized and State Feedback

th = linspace(0,2*pi) ; % the points for generating the plot of the circle mass
r = 0.01; % radius of circle
cx = l*cos(theta(:,1)); % x-axis location of the center of the circle
cy = l*sin(theta(:,1)); % y-axis location of the center of the circle
x = cx(1) + r*cos(th) ; 
y = cy(1) + r*sin(th) ;
bar = [0 0; cx(1) cy(1)];
figure;
h = plot(x,y,'r') ;
hold on
k = plot(bar(:,1),bar(:,2),'b');
axis equal 
%%animation
for i = 2:length(cx)
    bar = [0 0; cx(i) cy(i)];
    x = cx(i)+r*cos(th) ;
    y = cy(i)+r*sin(th) ;
    set(h,'XData',x,'YData',y) ;
    set(k,'XData',bar(:,1),'YData',bar(:,2));
    axis([-l l -l l]*2)    
    drawnow
    pause(0.0005)
end

%% Linearization and Integral State Space Feedback
% now adding in input to the system

% Linearized model:
m = 1; % kg
l = .1; % m
a = 9.8; % m/s^2
mu = 0.05; % friction factor (can also capture all losses in drive train)
A = [0 1; a/l 0];
B = [0; 1/(m*l^2)];
C = [1 0];
Aint = [A zeros(2,1); C 0];
Bint = [B; 0];
Cint = [C 0];
w = [0; 0]; % disturbances
yd = pi/2; % desired output of system
r = [w; -yd]; % reference point for integral feedback, coupled with disturbances

% State Feedback, placing poles at -1, -2, -3
poles = [-4; -6; -10];
Kint = acker(Aint,Bint,poles) % gains of controller u = -K*xtilde + r
Abar = Aint-Bint*Kint; % closed-loop A matrix
figure;
olsys = ss(A,B,C,0);
step(olsys)
figure;
clsys = ss(Abar,r,Cint,0);
step(clsys)

% simulating linearized integral CL system with initial conditions theta0
t = linspace(0,10,10^3)';
theta0 = [pi/2-pi/4; -pi; 0];
u = ones(length(t),1);
[y,t,theta] = lsim(clsys,u,t,theta0);
figure;
plot(t,theta)
legend('\theta',"\theta'","integral state")

%% Animation, Linearized and Integral State Feedback

th = linspace(0,2*pi) ; % the points for generating the plot of the circle mass
r = 0.01; % radius of circle
cx = l*cos(theta(:,1)); % x-axis location of the center of the circle
cy = l*sin(theta(:,1)); % y-axis location of the center of the circle
x = cx(1) + r*cos(th) ; 
y = cy(1) + r*sin(th) ;
bar = [0 0; cx(1) cy(1)];
figure;
h = plot(x,y,'r') ;
hold on
k = plot(bar(:,1),bar(:,2),'b');
axis equal 
%%animation
for i = 2:length(cx)
    bar = [0 0; cx(i) cy(i)];
    x = cx(i)+r*cos(th) ;
    y = cy(i)+r*sin(th) ;
    set(h,'XData',x,'YData',y) ;
    set(k,'XData',bar(:,1),'YData',bar(:,2));
    axis([-l l -l l]*2)    
    drawnow
    pause(0.0005)
end


