function PSA = responsespectrum(T,s,zi,dt)

na    = length(T);
nb    = length(s);
SD    = zeros(na,1);

for j = 1 : na
    om      = 2*pi/T(j);
    ub      = duhamel(om,zi,1,dt,nb,0,0,-s);
    SD(j)   = max(abs(ub));
end

PSA = (2*pi./T').^2 .* SD;
