%% Kinetics of serpentinisation
%
% This script computes the extent of serpentinisation as a function of
% time, including regular grain fragmentation events based on a crack
% growth model.

%% Setup paramters and solve

p = parameters('hsgrain',55e-6,...
    'hsgrainmin',10e-6,...
    'grainsizeprop',2.9,...
    'supcrtfile','SUPCRT/data2_P',...
    'ac0',0.8);

p.sigmainf = 0;

sol = reaction(p);

%% Compute reaction progress without fracturing
time_wf = linspace(0,max(sol.t*p.tau),1000);
k_eff   = 2*p.Vm_for/(p.Vm_liz+p.Vm_brc)*p.k;
xi_wf   = 1 - (p.hsgrain-2*k_eff*time_wf).^3./(p.hsgrain).^3;

%% plots

plot_expe;

figure;
subplot 211
plot(sol.t*p.tau/3600, sol.xi,'-');
hold on
plot(time_wf/3600, xi_wf,'k--');
plot(data_time,data_taux/100,'k.','MarkerSize',20)
xlim([0 16e3]);
ylim([0 1]);
xlabel('time (hours)');
ylabel('reaction progress (\xi)');
text(16000,0.9,' (a) ',...
    'verticalalignment','top',...
    'horizontalalignment','right');

subplot 212
plot(sol.t*p.tau/3600, sol.D*p.l*1e6,'-')
xlim([0 16e3]);
ylim([0 sol.D(1)*p.l*1e6]);
xlabel('time (hours)');
ylabel('grain size ({\mu}m)');
text(16000,55,' (b) ',...
    'verticalalignment','top',...
    'horizontalalignment','right');

exportfig('reaction', 'ySize',12, 'font','Helvetica','fontsize',8);