% data from malvoisin et al (2012) for comparison Grain size of 50 - 63
% microns
data_taux = [0
    3.0103
    3.4347
    4.4335
    5.4676
    6.1780
    7.7234
    11.3395
    14.5853
    18.0742
    18.5277
    20.8452
    28.0875
    34.6370
    39.6751
    45.7209
    57.6864
    62.4726
    71.9190
    74.5640
    76.8312
    79.9800];

data_time = [0
    0.1659
    0.1826
    0.2157
    0.2582
    0.2910
    0.3486
    0.5192
    0.5788
    0.6721
    0.6795
    0.7330
    0.8701
    0.9300
    0.9679
    1.0088
    1.0639
    1.0803
    1.1285
    1.1696
    1.1867
    1.2324].*10^4;
% close all,
% figure(1),hold on
% plot(data_time,data_taux,'k.','MarkerSize',20)
% grid on
% xlabel('Time (hr)')
% ylabel('Reaction progress (%)')
% xlim([0 16000])