function dtheta_dt=pend_nofriction(t,theta,a,m,l)
%First order representation of frictionless pendulum

dtheta_dt= zeros(2,1);
dtheta_dt(1) = theta(2);
dtheta_dt(2) = -a/l*cos(theta(1));

