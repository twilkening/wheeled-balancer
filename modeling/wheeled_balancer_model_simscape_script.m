%% Define Mechanical Characteristics
g = 9.81; % [m/sec^2]
m_motor = 9.5*2/1000; % [kg] documentation
m_bump = 35/1000; % [kg] documentation
m_batt = 138/1000; % [kg] measured
m_wheels = 39.7/1000; % [kg] documentation
m_cb = (200/1000 - m_motor - m_bump - m_wheels); % [kg] measured value

m_c = (m_wheels); % [kg] - mass of cart (i.e. the wheels of the segway)
m = m_motor + m_bump + m_batt + m_cb; % [kg] - mass of "pole"
mu_c = 0.02; % bearing viscous friction coeff. for rotational speed

%% System Model

syms g_sym mu_c_sym m_sym m_c_sym l_sym J_sym r_sym
% Linearized Model about (x,theta) = (0,0)
Asym = [0, 1, 0, 0; 
    0, -(g_sym*mu_c_sym*l_sym^2*m_sym^2 + g_sym*m_c_sym*mu_c_sym*l_sym^2*m_sym + J_sym*g_sym*mu_c_sym*m_sym +...
    J_sym*g_sym*m_c_sym*mu_c_sym)/(r_sym*(m_sym*m_c_sym*l_sym^2 + J_sym*m_sym + J_sym*m_c_sym)), ...
    -(g_sym*l_sym^2*m_sym^2)/(m_sym*m_c_sym*l_sym^2 + J_sym*m_sym + J_sym*m_c_sym), 0;
    0, 0, 0, 1; 
    0, (l_sym*m_sym*(g_sym*m_sym*mu_c_sym + g_sym*m_c_sym*mu_c_sym))/(r_sym*(m_sym*m_c_sym*l_sym^2 + J_sym*m_sym + J_sym*m_c_sym)),...
    (l_sym*m_sym*(g_sym*m_sym*r_sym + g_sym*m_c_sym*r_sym))/(r_sym*(m_sym*m_c_sym*l_sym^2 + J_sym*m_sym + J_sym*m_c_sym)), 0];
Bsym = [0;
    (m_sym*r_sym*l_sym^2 + J_sym*r_sym)/(r_sym*(m_sym*m_c_sym*l_sym^2 + J_sym*m_sym + J_sym*m_c_sym));
    0;
    -(l_sym*m_sym)/(m_sym*m_c_sym*l_sym^2 + J_sym*m_sym + J_sym*m_c_sym)];

%%

% Calculating CoM of the pole w.r.t. location of axis of rotation
x_rot = 31/1000; % [m] axis of rotation is NOT at the bottom of the pole
x_motor = 0/1000; % [m] motor w.r.t. rotation axis
l_cb = 109.2/1000; % [m] lengh of control board
x_cb = l_cb/2 - x_rot; % [m] control board w.r.t. rotation axis
l_batt = 95.7/1000; % [m] length of batteries
x_batt = l_batt/2 - x_rot + 2.5/1000; % [m] batteries w.r.t. rotation axis
x_bump = l_cb - 30/1000 + 40/1000 - x_rot; % [m] bumpers w.r.t. rotation axis
y_motor = 16.416/1000; y_cb = 0; y_batt = 0; y_bump = 0; % [m] w.r.t. rotation axis
CM_x_rot = (m_motor*x_motor + m_cb*x_cb + m_batt*x_batt + m_bump*x_bump)/(m); % [m]
CM_y_rot = y_motor*m_motor/m; % [m]
CM_rot = sqrt(CM_x_rot^2 + CM_y_rot^2); % [m] length from CM of pole to axis of rotation

l = CM_rot; % [m] - length to pole CoM

% Calculating Inertia of the Pole:
I_motor = m_motor*y_motor^2; % [kg m^2] point mass inertia of motors
I_bump = m_bump*x_bump^2; % [kg m^2] point mass inertia of bumpers
% Parallel axis theorem: I = Icm + m h^2
% Inertia of a pole of equally distributed mass about one end: I = m l^2/3 
% Inertia of a pole of equally distributed mass about the CoM: I = m l^2/12
I_batt = m_batt*(l_batt/2)^2/12 + m_batt*x_batt^2; % [kg m^2] inertia of batteries
I_cb = m_cb*(l_cb/2)^2/12 + m_cb*x_cb^2; % [kg m^2] inertia of control board

J = (I_motor + I_bump + I_batt + I_cb); % kg*m^2 - mass moment of inertia of the pole (1/12*m*(length)^2)

r = 40/1000; % [m]

A = double(subs(Asym,[g_sym,m_c_sym,m_sym,mu_c_sym,l_sym,J_sym,r_sym],...
    [g,m_c,m,mu_c,l,J,r]));
B = double(subs(Bsym,[g_sym,m_c_sym,m_sym,mu_c_sym,l_sym,J_sym,r_sym],...
    [g,m_c,m,mu_c,l,J,r]));
% outputting both the x-position and theta_dot angle-rate
C = [1 0 0 0; 0 0 0 1];
D = [0;0];
states = {'x' 'x_dot' 'theta' 'theta_dot'};
inputs = {'u'};
outputs = {'x'; 'theta'};
plant_continuous = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);

%% Discretize System
T = 1/100; % 0.01s, 100Hz
plant = c2d(plant_continuous,T);
[G, H, Cd, Dd] = ssdata(plant);

%% Define Control Gains

pc_control = [-3, -4, -5, -6, -7]; % add one pole for integral state
pd_control = exp(T*pc_control);
pc_observer = [-20,-21,-22,-23]*2;
pd_observer = exp(T*pc_observer);

C2 = [1 0];
Cprime = C2*Cd;
Gbar = [G, zeros(4,1); -Cprime, 1];
Hbar = [H; 0];
K = acker(Gbar, Hbar, pd_control); Ks = K(1:end-1), Ki = K(end),
L = place(G',Cd',pd_observer)'

%% Plot after simulation

% x - xhat = xtilde (error)
xsim = [out.x.data, out.v.data, out.x.data - out.xhat.data];
figure; t2 = myplot2(out.x.Time, xsim);
title(t2, {"$\bf I.C.\ Response$"; "From: u"}, 'Interpreter', 'Latex')

usim = -Ks*out.xhat.data' -Ki*out.v.data';
figure; stairs(out.x.Time,usim); ylabel('Newtons');
title('Control Effort');