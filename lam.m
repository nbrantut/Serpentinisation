function out=lam(ac)

out = 1/2*(0.13448e2 / pi * (asin(ac) + ac .* acos(ac)) - 0.1568e1 * ac + 0.56e0 * ac.^2 + 0.100e0 * ac.^3 );
