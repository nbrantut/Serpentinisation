function df = f(x,param)
%df = f(x,param)
%
%This is the function governing the ODE for crack aand wedge sizes:
%  {a',c',w'} = f(t, {a,c,w})
%And we include its parameters as the structure 'param'.

% parameters
mu = param.mu;
theta = param.theta;
E = param.Enorm;
soo = param.sigmainf;

% extract variable to ease manipulation
a = x(1);
c = x(2);
w = x(3);

% compute useful parameters
g = gam(a/c);
l = lam(a/c);
[dga, dgc] = dgam(a,c);
[dla, dlc] = dlam(a,c);

% compute normal stress
sn = (1/l)*(E*w/c  + soo*theta);

% energy release rate
G = c*(g/mu*sn-soo)^2;

dw = 1 - sn;
da = 1;

if G<1
    if a<c
        dc = 0;
    else
        dc=da;
    end
else
    if a>=c
        %recompute the derivatives with c slightly larger than a to avoid
        %encountering lots of NaNs:
        [dga, dgc] = dgam(a,c*(1+1e-8));
        [dla, dlc] = dlam(a,c*(1+1e-8));
    end
    
    dc = 2*c*((-dla*g + dga*l)*(theta*soo*c + E*w)*da + E*g*l*dw)/...
            (l*((-theta*soo*c+E*w)*g + mu*soo*c*l ...
            -2*c*(theta*soo*c + E*w)*dgc) + 2*c*(theta*soo*c + E*w)*g*dlc);
end

% if a<c     
%     if G<1   % condition for fracture propagation
%         dc = 0;
%     else     % from dG = 0
%         dc = 2*c*((dla*g - dga*l)*(theta*soo*c - E*w)*da + E*g*l*dw)/...
%             (soo*c*(l*(g*theta+l*mu)+2*theta*(dgc*l-dlc*g)*c) + ...
%             E*(l*g+2*(dlc*g-dgc*l)*c)*w);
%     end
% else  % a = c
%     if G<1   % condition for fracture propagation
%         dc =da; 
%     else      
%         dc = (2*E*g*l*c*dw)/...
%             (soo*c*(l*(g*theta+l*mu)+2*theta*(-(dla+dlc)*g+(dga+dgc)*l)*c) + ...
%             E*(l*g+2*((dla+dlc)*g-(dga+dgc)*l)*c)*w);
%     end
% end

df=[da; dc; dw];