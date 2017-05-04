function out = gam(ac)

out = asin(ac) .* (0.13e1 - 0.143e0 .* ac - 0.12e0 * ac .^ 2 + 0.83e-1 * ac .^ 3);