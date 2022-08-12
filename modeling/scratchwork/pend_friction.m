function dtheta_dt=pend_friction(t,theta,a,m,l,mu)
%First order representation of frictionless pendulum

dtheta_dt= zeros(2,1);
dtheta_dt(1) = theta(2);
dtheta_dt(2) = -(1)*(a/l*cos(theta(1))) - mu*a/l*sign(theta(2));

