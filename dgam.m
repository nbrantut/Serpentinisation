function [dgda,dgdc] = dgam(a,c)


ac=a./c;

gammaprime = (1.3-0.143*ac-0.12*ac.^2+0.083*ac.^3)./sqrt(1-ac.^2) +...
    (-0.143-0.24*ac+0.249*ac.^2)*asin(ac);

dgda = 1./c .*gammaprime;
dgdc = -a./(c.^2) .*gammaprime;