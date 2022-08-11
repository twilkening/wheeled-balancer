%% Cart-Pendulum Problem
% The equations of motion for the cart and the pendulum are as follows,
% where:
% g = gravity = 9.8 m/s^2
% I = inertial constant of the pendulum
% b = damping coefficient of friction of the cart
% x = cart position
% th = theta = pendulum angle w.r.t. vertical down
syms M m x_ddot x_dot x theta_ddot theta_dot theta l b I g F
eqn_cart = M*x_ddot == F - m*x_ddot - m*l*theta_ddot*cos(theta) +...
    m*l*theta_dot^2 - b*x_dot
eqn_pend = theta_ddot*(m*l^2 + I) == -m*g*l*sin(theta) -...
    m*l*cos(theta)*x_ddot

% Deriving the solutions for ONLY x_ddot and theta_ddot
[x_sol, theta_sol] = solve([eqn_pend,eqn_cart],[x_ddot,theta_ddot])


% % Setting up the State-space equations
syms beta_1 beta_2 beta_3 beta_4 u t
% f = [beta_1_dot; beta_2_dot; beta_3_dot; beta_4_dot] 
%   = [x_dot; x_ddot; theta_dot; theta_ddot]
f = [x_dot;
    x_sol;
    theta_dot;
    theta_sol]
f = subs(f,[x,x_dot,theta,theta_dot,F],[beta_1,beta_2,beta_3,beta_4,u])
beta = [beta_1; beta_2; beta_3; beta_4]
Asym = jacobian(f,beta)
Bsym = jacobian(f,u)
% system was linearized at equil. points (no movement), where theta=0,pi
% and x is not constrained. Thus, substituting in the equil point where
% theta = pi ;
Api = subs(Asym,[beta_1, beta_2, beta_3, beta_4], [0, 0, pi, 0])
Bpi = subs(Bsym,[beta_1, beta_2, beta_3, beta_4], [0, 0, pi, 0])

% Thus, our linearized state-space representation of the system can be
% represented as: diff(beta,t) = Api*beta + Bpi*u

% where beta_3 is the deviation from the equil point, phi


% % Inserting real values and creating State Space object
g_real = 9.8; % m/sec^2
M_real = 0.5; % kg - mass of cart
m_real = 0.2; % kg - mass of pendulum
b_real = 0.1; % N/m/sec
l_real = 0.3; % m - length to pendulum CoM
% kg*m^2 - mass moment of inertia of the pendulum (1/12*m*(length)^2):
I_real = 0.006; 
% for inertial calculation, assume uniform rod of negligible thickness,
% with inertia calculated about it's center. 
% reference: http://hyperphysics.phy-astr.gsu.edu/hbase/hframe.html
A = double(subs(Api,[g,M,m,b,l,I],[g_real,M_real,m_real,b_real,...
    l_real,I_real]));
B = double(subs(Bpi,[g,M,m,b,l,I],[g_real,M_real,m_real,b_real,...
    l_real,I_real]));
% outputting both the x-position and theta-angle
C = [1 0 0 0; 0 0 1 0];
D = [0;0];
states = {'x' 'x_dot' 'phi' 'phi_dot'};
inputs = {'u'};
outputs = {'x'; 'theta'};
olsys = ss(A,B,C,D,'statename',states,'inputname',inputs,...
    'outputname',outputs)
% initial(olsys,[0;0;7*pi/8;0],10) % this system blows up because we
% LINEARIZED around pi (vertical up), not 0 (vertical down)
olsys_tf = tf(olsys)


% % Stabilizing the System with SS Feedback gain
figure;
rlocus(olsys_tf(2))
% from the root locus plot theta/u, we see that there will always be a 
% pole in the RHP, so we have to do a state-space feedback to stabilize
% the system

% seeing if the system is controllable and observable (full rank = yes)
rank(ctrb(A,B))
rank(obsv(A,C))

p = [-3, -4, -5, -6];
K = acker(A,B,p)
Abar = A-B*K;
clsys = ss(Abar,B,C,D,'statename',states,'inputname',inputs,...
    'outputname',outputs)
[y_init,t_init,beta_state_init] = initial(clsys,[0;0;pi/8;0],10);
figure;
tiledlayout(2,1);
nexttile;
plot(t_init, y_init(:,1)); title('Response to Initial Conditions');...
    ylabel('To: x');
grid on
nexttile;
plot(t_init,y_init(:,2)); ylabel('To: \theta'); xlabel('Time (seconds)')
grid on


% % Desired output tracking
yd = [1; 0]; % x position: 1m, phi: 0rad deviation
clsys_gain = -C/Abar*B % this is the gain of the system with a step input
% from clsys_gain, we can see that phi will always return to 0, so we can
% only control the output x:
r1 = yd(1)/clsys_gain(1); % reference input

figure;
% step response with defined input so that output is at desired location
[y,t,beta_state] = step(clsys*r1);

%% Animating the desired output tracking - NEXT STEP
x = beta_state_init(:,1);
theta = beta_state_init(:,3);
slen = 0.2;
sq = [  -slen  slen   slen   -slen  -slen;
        0   0   slen   slen   0]; % the vertex locations of the cart
figure;
cart = plot(sq(1,:),sq(2,:),'r');
hold on

th = linspace(0,2*pi) ; % the points for generating the plot of a circle
r1 = 0.01; % radius of circle 1 (CoM)..?

cB = [0; l_real]; % initial x-axis location of the pendulum CoM, B-frame
eB = 2*cB; % initial x-axis location of the END of the pendulum 
% rotational matrix
R_AB = [cos(theta(1)) sin(theta(1)); -sin(theta(1)) cos(theta(1))]; 
d_AB = [0; slen]; % location of center of B-frame w.r.t. A-frame
g_AB = [R_AB d_AB; 0 0 1]; % transformation matrix
cA = g_AB*[cB; 1]; cA = cA(1:2);

eA = g_AB*[eB; 1]; eA = eA(1:2);

CoMx = cA(1) + r1*cos(th) ; 
CoMy = cA(2) + r1*sin(th) ;

CoM = plot(CoMx,CoMy,'bl') ;
bar = plot([d_AB(1) eA(1)], [d_AB(2) eA(2)],'b')
axis equal
minx = min(x)-slen;
maxx = max(x)+slen;
axis([minx maxx 0 1]) 

% animation loop
for i = 2:length(x)/2
    sq_draw = sq + [repmat(x(i),[1,size(sq,2)]); zeros(1,size(sq,2))];
    
    % rotational matrix
    R_AB = [cos(theta(i)) -sin(theta(i)); sin(theta(i)) cos(theta(i))]; 
    d_AB = [x(i); slen]; % location of center of B-frame w.r.t. A-frame
    g_AB = [R_AB d_AB; 0 0 1]; % transformation matrix
    cA = g_AB*[cB; 1]; cA = cA(1:2);
    eA = g_AB*[eB; 1]; eA = eA(1:2);
    
    CoMx = cA(1) + r1*cos(th) ; 
    CoMy = cA(2) + r1*sin(th) ;
    bar_x = [d_AB(1) eA(1)];
    bar_y = [d_AB(2) eA(2)];

    set(cart,'XData',sq_draw(1,:),'YData',sq_draw(2,:)) ;
    set(CoM,'XData',CoMx,'YData',CoMy) ;
    set(bar,'XData',bar_x,'YData',bar_y) ;
    axis([minx maxx 0 1])     
    drawnow
    pause(0.05)
end

% Done: Figure out why the angle is negative from the simulation...
% (it was my definition of the rotation matrix)

% NEXT STEPS:
% get the integral control working, 
% then do an optimal control solution