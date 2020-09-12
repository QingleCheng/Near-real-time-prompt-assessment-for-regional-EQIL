function [u] = duhamel(wn,zi,m,dt,nt,u0,v0,f)


% this function was developed by Dr. Luis Suarez, Professor at the
% department of Civil Engineering at the University of Puerto Rico at
% Mayaguez

% Program to perform the numerical integration of the equation of motion
% of a viscous single dof linear system using the recursive equation based 
% on the Duhamel's integral.
% The system is subjected to an excitation f(t) with arbitrary time 
% variation sampled at constant intervals.
% The excitation is assumed to vary linearly on each time interval 

%------------------------ Input data ----------------------------------
% wn  = natural frequency in rad/sec
% zi  = damping ratio for underdamped case
% m   = mass of oscillator
% dt  = constant time step
% nt  = number of time steps
% u0  = initial displacement at time t = 0
% v0  = initial velocity at time t = 0
% f   = vector with force time history: nt points

zk1 = [0 ; 0];
zk  = [u0 ; v0];

wd  = wn * sqrt(1-zi^2) ;
ex  = exp(-zi*wn*dt) ;
co  = cos(wd*dt) ;
si  = sin(wd*dt) ;
S   = sqrt(1-zi^2);
k   = 1 / (m * wn^3 * S * dt);

a11 = (co + zi/S * si) * ex ;
a12 = 1 / wd * si * ex ;
a21 = -wn/S * si * ex ;
a22 = (co - zi/S * si) * ex ;

b11 = 2*zi*S + ( (1-2*zi^2-zi*wn*dt)*si - (2*zi*S+wd*dt)*co ) * ex; 
b12 = wd*dt - 2*zi*S + ( (2*zi^2-1)*si + (2*zi*S)*co ) * ex;
b21 = -wd + ( (wn*zi+wn^2*dt)*si + (wd)*co ) * ex;
b22 = wd - ( (wn*zi)*si + (wd)*co ) * ex;

A1 =     [a11, a12 ; a21, a22] ;
B1 = k * [b11, b12 ; b21, b22] ;

u(1) = u0;
for k = 2:nt 
    zk1  = A1 * zk + B1 * [f(k-1) ; f(k)] ;
    u(k) = zk1(1); 
    zk   = zk1;
end

