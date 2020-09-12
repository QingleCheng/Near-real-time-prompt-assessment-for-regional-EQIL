function wv = zumontw(t,omega,zeta)

% Generates the Suarez-Montejo Wavelet function
% Ref. Suarez, L.E. & Montejo, L.A. Generation of artificial earthquakes 
% via the wavelet transform, Int. Journal of Solids and Structures, 42,
% 2005
% wv = exp(-zeta*omega*abs(t))*sin(omega*t)
wv = exp(-zeta*omega*abs(t)).*sin(omega*t);
