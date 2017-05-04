%% Crack growth during serpentinisation
%
% This script computes the dimensions of an edge crack that propagates due
% to a crystallisation force applied over a a short distance at its base.
% The whole crack is assumed to be under a compressive load which tends to
% stop its propagation as it grows.

%% Setup parameters and options
p = parameters('hsgrain',100e-6,...
    'grainsizeprop',2.9,...
    'supcrtfile','SUPCRT/data2_P',...
    'ac0',0.8);

options = odeset(...%'outputfcn', @odeplot,...
    'Events',@(t,x) grainsplit(x(2),p.hsgrainnorm/p.grainsizeprop),...
    'abstol',1e-8,...
    'reltol',1e-8);

%% Solve and extract derived variables

sol = ode45(@(t,x) f(x,p),[0 30], [p.a0;p.cminnorm*1.05;p.w0] ,options);

t=sol.x;
a=sol.y(1,:);
c=sol.y(2,:);
w=sol.y(3,:);

%derived variables:

sn = (1./lam(a./c)).*(p.Enorm*w./c  + p.sigmainf*p.theta);
G = c.*(gam(a./c)/p.mu.*sn-p.sigmainf).^2;
KI = sqrt(c).*(gam(a./c)/p.mu.*sn-p.sigmainf);

%detect when crack starts growing beyond a:

i0=find(c-c(1),1,'first');

%% plots
subplot 311
plot(t*p.tau/3600, a*p.l/1e-6, 'k');
hold on;
plot(t*p.tau/3600, c*p.l/1e-6, 'k');
plot(t*p.tau/3600, w*100*p.l/1e-6,'k');
xlim([0, p.tau*t(end)/3600]);
ylabel('dimensions (\mum)');

set(gca, 'xtick',[],...
    'position',get(gca,'position').*[1.1 1 1 1]);

text(8000, 29,'{\itc}',...
    'verticalalignment','bottom');
text(8000, 15,'{\ita}',...
    'verticalalignment','bottom');
text(8000, 5,'100 {\itw}',...
    'verticalalignment','bottom');

text(0, 40, ' (a)',...
    'verticalalignment','top');

subplot 312
plot(t*p.tau/3600,G*p.Gc,'k');
hold on
plot([0, t(end)*p.tau/3600], p.Gc*[1 1],'k:');
xlim([0, t(end)*p.tau/3600]);
ylim([0 3]);
ylabel({'energy release rate (J m^{-2})'});

text(8000, 2,'{\itG}_c',...
    'verticalalignment','bottom');
text(120, 1.3, '{\itG}');

text(0, 3, ' (b)',...
    'verticalalignment','top');

set(gca, 'xtick',[],...
    'position',get(gca,'position').*[1.1 1.1 1 1]);

subplot 313
plot(t(1:end)*p.tau/3600, sn(1:end)*p.sc/1e6,'k');
hold on
plot([0, t(end)*p.tau/3600], p.sc*[1 1]/1e6,'k:');
xlim([0, t(end)*p.tau/3600]);
ylim([sn(end)*0.99  max(sn)*1.01]*p.sc/1e6);
xlabel('time {\itt} (hours)');
ylabel('normal stress (MPa)');
text(8500, p.sc/1e6,'{\sigma}_c',...
    'verticalalignment','bottom');
text(1000, 165.25, '\sigma_n',...
    'verticalalignment','top');

text(0, max(sn)*1.01*p.sc/1e6, ' (c)',...
    'verticalalignment','top');


set(gca, 'position',get(gca,'position').*[1.1 1.6 1 1]);

exportfig('crackgrowth','ySize',13, 'font','Helvetica','fontsize',8)