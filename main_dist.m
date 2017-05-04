%% evolution of grain size with an initial grain distribution from statistics 
clear all
close all
load K3a_zone3 % load of statistic results on sample K3a_zone4 
dm = (dgrand+dpetit)/2; % initial grain size
ngg = length(dm);%length(dm);     % number of grain for the modelling (<length(dm)) 
gs = datasample(dm,ngg);% random sampling of the data set
mol_frac = gs.^3./(sum(gs.^3)); % volume/molar fraction of each grain for xitot


%% run the microcracks model for determining the maximum size of xi and time
p = parameters('hsgrain',max(gs)*10^-6,...
        'hsgrainmin',10e-6,...
        'grainsizeprop',2.9,...
        'supcrtfile','SUPCRT/data2_P',...
        'ac0',0.8);
p.sigmainf = 0; 
sol = reaction(p);

max_t = sol.t(end)*p.tau/3600; % maximum time for the longest simulation
dt = max_t/1000; % timestep
xitot = zeros(1,length([0:dt:max_t])); % initialize xitot
timetot = [0:dt:max_t]; % vector for the time

for i = 1:ngg
    disp(i)
    clear p sol
    %% run the microcracks code for each grain
    p = parameters('hsgrain',gs(i)*10^-6,...
        'hsgrainmin',10e-6,...
        'grainsizeprop',2.9,...
        'supcrtfile','SUPCRT/data2_P',...
        'ac0',0.8);
    p.sigmainf = 0;
    sol = reaction(p);
    
    %% save the results
    [time{i},ix,ib] = unique(sol.t*p.tau/(3600)); % time
    xi{i} = sol.xi(ix);% reaction progress
    d{i} = sol.D(ix)*p.l*1e6; % grain size
    ng{i} = p.Nfold.^sol.frag(ix); % number of grains
    
    %% make the total parameters
    % boundaries for having the same time scales
    if max(time{i})<max_t
        lr = 1001-length([0:dt:max(time{i})]);
        lrm = max(time{i});
    elseif max(time{i})>max_t
        lr = 0;
        lrm = max_t;
    else
        lrm = max(time{i});
        lr = 0;
    end
    xitot = xitot + [interp1(time{i},xi{i},[0:dt:lrm]) ones(1,lr)]*mol_frac(i); % total reaction progress
    dtot(i,:) = [interp1(time{i},d{i},[0:dt:lrm]) zeros(1,lr)]; % matrix with all the grain size
    ngtot(i,:) = [interp1(time{i},ng{i},[0:dt:lrm]) zeros(1,lr)]; % matrix with the number of grains
    
end

%% make the weighted mean as in 2D so surfaces 
dtot2 = dtot;
dtot2(dtot2<5)=NaN; % excluding grains smaller than 5 mum
pmean = nansum(ngtot.*dtot2.^4,1)./nansum(ngtot.*dtot2.^3,1);
%% total number of grains
ngtot2 = ngtot;
ngtot2(dtot<5) = NaN;% excluding grains smaller than 5 mum
ngraint = nansum(ngtot2,1);

%% make the figure
figure(1)
subplot(131)
plot(timetot/24,xitot,'-', 'color',0.9*[1 1 1]);
hold on
plot(timetot/24,xitot,'k.');
xlabel('time (days)');
ylabel('reaction progress (\xi)'); 
text(0,1,' (a)',...
    'horizontalalignment','left',...
    'verticalalignment','top');
set(gca,'xlim',[0 4500])

subplot(132)
plot(xitot,pmean,'-', 'color',0.9*[1 1 1]);
hold on
plot(xitot,pmean,'k.')
xlabel('reaction progress (\xi)');
ylabel('weigted mean grain size (\mum)');
set(gca, 'xtick',[0:.25:1]);
ylim([0 1000])
text(0,1000,' (b)',...
    'horizontalalignment','left',...
    'verticalalignment','top');

subplot(133)
plot(xitot,ngraint/1000,'-', 'color',0.9*[1 1 1]);
hold on
plot(xitot,ngraint/1000,'k.')
xlabel('reaction progress (\xi)');
ylabel('number of grains ({\times}1000)');
set(gca, 'xtick',[0:.25:1]);
ylim([0 2100]);
text(0,2100,' (c)',...
    'horizontalalignment','left',...
    'verticalalignment','top');

exportfig('upscaling', 'xSize',19, 'ysize',6.5,'font','Helvetica','fontsize',8);