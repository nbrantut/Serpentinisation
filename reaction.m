function sol = reaction(param)
%sol = reaction(param)
%
%This function solves the reaction kinetics problem based on the input
%parameters contained in the structure 'param'.
%
%The output structure contains the solution for grain size, reaction
%progress and number of fragmentation events as a function of time, all in
%normalised quantities.

%notation:
Di = param.hsgrainnorm;

%solve the crack growth problem
options = odeset(...
    'Events',@(t,x) grainsplit(x(2),Di/param.grainsizeprop),...
    'abstol',1e-8,...
    'reltol',1e-8);

soln = ode45(@(t,x) f(x,param),[0 30], [param.a0;param.cminnorm*1.05;param.w0] ,options);

time=soln.x;
c=soln.y(2,:);

%speed of olivine interface during reaction to serpentine:
v = (2)*(2*param.Vm_for/(param.Vm_liz+param.Vm_brc))*param.k/param.l*param.tau; %normalised !
%factor 1/2 is because both sides of the grain are attacked at the same
%time !

%initialise:
Nk = 100; %guess number of timesteps
Dk = zeros(Nk,1);
fragk = Dk;
xik = Dk;
tk = Dk;

%initial condition:
Dk(1) = Di;
xik(1) = 0;

n=1;

check = 1;

while check
%     n=n+1;
%     
%     %grain size reduction threshold:
%     dD = min(param.cminnorm, Dk(n-1));
%     %increment time and grain size etc:
%     tk(n) = tk(n-1) + dD/v;
%     Dk(n) = Dk(n-1) - dD;
%     fragk(n) = fragk(n-1);
%     xik(n) = 1 - 8.^fragk(n).*(Dk(n)/Di)^3;

    %if we can fracture (grain size not too small):
    if Dk(n)>param.hsgrainminnorm && Dk(n)>(param.grainsizeprop)*param.cminnorm*1.05
        
        % look for the time at which we fracture:
        [cmind,i0,~] = unique(c - (1/param.grainsizeprop)*(Dk(n) - v*time));
        %plot(time,c, time, (1/param.grainsizeprop)*(Dk(n) - v*time));
        %drawnow;
        %pause;
        tcrack = interp1(cmind, time(i0), 0);

        %increment time, grain size etc vectors:
        n=n+1;
        tk(n) = tk(n-1) + tcrack;
        Dk(n) = Dk(n-1) - v*tcrack;
        fragk(n) = fragk(n-1);
        xik(n) = 1 - param.Nfold.^fragk(n).*(Dk(n)/Di)^3;

        %fracture happens and set D <-- D/2
        n=n+1;
        tk(n) = tk(n-1);
        Dk(n) = Dk(n-1)/(param.Nfold^(1/3));
        fragk(n) = fragk(n-1)+1;
        xik(n) = 1 - param.Nfold.^fragk(n).*(Dk(n)/Di)^3;
    else
        n=n+1;
        %grain size reduction threshold:
        dD = Dk(n-1);
        %increment time and grain size etc:
        tk(n) = tk(n-1) + dD/v;
        Dk(n) = Dk(n-1) - dD;
        fragk(n) = fragk(n-1);
        xik(n) = 1 - param.Nfold.^fragk(n).*(Dk(n)/Di)^3;
    end
    check = xik(n)<1;
end

%remove extra stuff
tk(n+1:end)  = [];
xik(n+1:end) = [];
Dk(n+1:end)  = [];
fragk(n+1:end) = [];

%assign results
sol.raw.t = tk;
sol.raw.xi = xik;
sol.raw.D = Dk;
sol.raw.frag = fragk;

%now make nicer curves by interpolating:
nsub = 64;
N = nsub*(length(tk) - 1);

sol.t = zeros(N,1);
sol.D = zeros(N,1);
sol.frag = zeros(N,1);
inc = 0;

for k=1:length(tk)-1
    if tk(k+1)>tk(k)
        t_in = linspace(tk(k), tk(k+1), nsub);
        D_in = interp1(tk(k:k+1), Dk(k:k+1), t_in);
        frag_in = t_in*0 + fragk(k);
        
        inc = inc+1;
        
        sol.t(1+(inc-1)*nsub:inc*nsub) = t_in;
        sol.D(1+(inc-1)*nsub:inc*nsub) = D_in;
        sol.frag(1+(inc-1)*nsub:inc*nsub) = frag_in;
    end
end

sol.t(inc*nsub+1:end)=[];
sol.D(inc*nsub+1:end)=[];
sol.frag(inc*nsub+1:end)=[];

sol.xi = 1 - param.Nfold.^sol.frag.*(sol.D/Di).^3;




    
    
