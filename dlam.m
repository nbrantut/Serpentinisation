function [dlda,dldc] = dlam(a,c)


ac=a./c;

lambdaprime = 2*(-0.392+0.28*ac+0.075*ac.^2 + ...
    1.07016.*(1./sqrt(1-ac.^2) - ac./sqrt(1-ac.^2) + acos(ac)));

dlda = 1./c .*lambdaprime;
dldc = -a./(c.^2) .*lambdaprime;