param = parameters();

A = textread(param.supcrtfile);
P_vec = A(:,1)/10*1e6;    % pressure in megapascal
T_vec = A(:,2)+273.15;    % temperature in K
omega_vec = 10.^(A(:,4)); % supersaturation index
sc_vec = 8.314.*T_vec.*log(omega_vec)./param.Vm_liz; % crystallisation pressure

%make interpolation table:
T_vec  = unique(T_vec)-273.15;
P_vec  = unique(P_vec);
sc_tab = reshape(sc_vec,length(T_vec),length(P_vec));

cmin = @(sc,spinf) (gam(1)*sc-param.mu*spinf>0).*param.Gc*pi*param.Eprime./(4*(gam(1)*sc-param.mu*spinf).^2);

%% plots
figure;
plot(T_vec, sc_tab(:,[10  end])/1e6,'k');
xlim([T_vec(1)  T_vec(end)]);
xlabel('temperature ({\circ}C)');
ylabel('crystallisation stress \sigma_c (MPa)')
text(T_vec(85), sc_tab(85,10)/1e6,...
    ['{\itP}='  num2str(P_vec(10)/1e6) ' MPa'],...
    'horizontalalignment','right',...
    'verticalalignment','top');
text(T_vec(65), sc_tab(65,end)/1e6,...
    ['{\itP}='  num2str(P_vec(end)/1e6) ' MPa'],...
    'horizontalalignment','left',...
    'verticalalignment','bottom');

exportfig('stress', 'font','Helvetica','fontsize',8);

figure;
spinf_tab = [0 100 200]*1e6;
i0 = [90 60 50];
for k=1:3
    semilogy(T_vec,  cmin(sc_tab(:,37),spinf_tab(k))*1e6,'k');
    hold on;
end
text(T_vec(80),  cmin(sc_tab(80,37),spinf_tab(1))*1e6,...
    {[num2str(spinf_tab(1)/1e6)]},...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top');
text(T_vec(65),  cmin(sc_tab(65,37),spinf_tab(2))*1e6,...
    {[num2str(spinf_tab(2)/1e6)]},...
    'HorizontalAlignment','right',...
    'VerticalAlignment','bottom');
text(T_vec(50),  cmin(sc_tab(50,37),spinf_tab(3))*1e6,...
    {['\sigma''_\infty='  num2str(spinf_tab(3)/1e6) ' MPa']},...
    'HorizontalAlignment','right',...
    'VerticalAlignment','bottom');

xlim([T_vec(1)  T_vec(end)]);
ylim([1 1e3])
text(T_vec(1), 1e3, [' {\itP}='  num2str(P_vec(37)/1e6) ' MPa'],...
    'horizontalalignment','left',...
    'verticalalignment','top');
xlabel('temperature ({\circ}C)');
ylabel('minimum crack size {\itc}_{min} (\mum)');

exportfig('cmin', 'font','Helvetica','fontsize',8)