%% Program to create and analyze auralization
% Created Feb-May 2021
% Author: shreejay
% shreejayshrestha@gmail.com
%%
clc
clear variables
%Load force signals
load('D:\NTNU\MCT2018_2020\MasterThesis2021\LABMeasurements\Clean_measurement_set_for_Analysis\Clean_Output_Good_For_Analysis\data_auralization\force_aur.mat')
%Load soundpressure signals
load('D:\NTNU\MCT2018_2020\MasterThesis2021\LABMeasurements\Clean_measurement_set_for_Analysis\Clean_Output_Good_For_Analysis\data_auralization\sp_aur.mat')
%Load impulse response EQulization filter
load('D:\NTNU\MCT2018_2020\MasterThesis2021\LABMeasurements\Clean_measurement_set_for_Analysis\Clean_Output_Good_For_Analysis\data_auralization\irEQ.mat')

%% Calculate impulse response of the coupled system
irsys = contwo(irEQ,sp_aur{1}{1,1});

%% selecting force and sound pressure signals
% From footstep_with_shoes, neglect : P1 rep 2, P3-Rep 2,5 and 6
% From footstep_without_shoes, neglect : P2 rep 3
% Reason for neglecting them: They were excluded in the main correcpoiding
% analysis scripts. mic_footstep_without_shoes_nofilt.m and
% mic_footstep_withshoes_nofilt.m
%% For simplicity select Rep # 7 for all signals.
% use the following function to generate auralization in time and freq.
% domain. 
% function [p_aur,F_orig,F_aur,F_aur_shift,ecf] = auralise_calc_empfact_v2(irsys,forcesig,soundpressure,refp)
% refp = select reference repetition # for estimating the
% auralization.

repnum = 2;
[p_aur,p_aur_concat,F_orig,F_aur,F_aur_shift,ecf] = auralise_calc_empfact_v3(irsys,force_aur,sp_aur,repnum);

%% Play auralized sounds
fs = 48000;
% sound(p_aur_concat{1,3},fs);

%% Plot Auralization Spectrum vs original sound pressure spectrum
% function fvec = plot_aura_sps(F_orig,F_aur_shift,refsource); % refsource = [1,2,3]
fig = plot_aura_sps(F_orig,F_aur_shift,1,repnum); % refsource = [1,2,3]

%% Perform Octaveband analysis: using octavebandanalysis.m function from PS.
fs = 48000;
nfft = 2^19;
iv = 1:nfft/2;
fvec= fs/nfft*(0:nfft/2-1);
F_aur_3oct = cell(3,5);
F_orig_3oct = cell(3,5);
fvecout_orig = cell(3,5);
fvecout_shift = cell(3,5);
for i = 1:3
    for j = 1:5
        [F_orig_3oct{i,j},fvecout_orig{i,j}] = octavebandanalysis(fvec,abs(F_orig{i,j}(iv)).^2,3,20,500);
        [F_aur_3oct{i,j},fvecout_shift{i,j}] = octavebandanalysis(fvec,abs(F_aur_shift{i,j}(iv)).^2,3,20,500);
    end
end

figure()
plot(fvecout_orig{1},10*log10([F_aur_3oct{1} F_aur_3oct{2} F_aur_3oct{3}]))
% plot(fvecout{1},20*log10([F_aur_spec{1} F_aur_spec{2} F_aur_spec{3}])./p0)
legend('Medisin ball','Footstep with shoes','Footstep without shoes')
xlim([20 200])
xticks([20	25	32	40	50	63	79	100	126	159	200])
xticklabels({'20', '25.1', '31.6', '39.8', '50.1', '63.1', '79.4', '100','125.9','158.5','199.5'})
grid on
g = xlabel('Frequency [ Hz ]');
set(g,'FontSize',16)
g = ylabel('Amplitude dB [ re uncal. ]');
set(g,'FontSize',16)
g = title('Auralized Sound Pressure in 1/3rd Octave Band');
set(g,'FontSize',16)
ax = gca;
ax.FontSize = 12;
%5 plot option 2
%% PLOT Uncalibrated auralised Vs original sps
linS = {'-','--','-','--','-','--'};
markS = {'d','*','o','s','x','p'};
col = [{'b'},{'g'},{'m'}];
legentry = [{'Medisin ball original'},{'Medisin ball auralised'},{'Footstep with shoes original'},{'Footstep with shoes auralised'},{'Footstep without shoes original'},{'Footstep without shoes auralised'}];
ep = 4;
figure()
hold on
for c = 1:3
    plot(fvecout_orig{1},8*log10(F_orig_3oct{c,ep}),'-','Marker',markS{c},'MarkerSize',7,'LineWidth',1,'Color',col{c})
    plot(fvecout_orig{1},8*log10(F_aur_3oct{c,ep}),'--','Marker',markS{c},'MarkerSize',7,'LineWidth',1,'Color',col{c})

end
grid on
ax = gca;
ax.FontSize = 12;
legend([legentry])
g = xlabel('Frequency [ Hz ]');
set(g,'FontSize',16)
g = ylabel('Amplitude dB [ uncal. ]');
set(g,'FontSize',16)
g = title('Auralized vs original Sound Pressure in 1/3rd Octave Band');
set(g,'FontSize',16)
ylim ([30 90])
xlim([20 200])
xticks([20.0	25.1	31.6	39.8	50.1	63.1	79.4	100.0	125.9	158.5	199.5])
xticklabels({'20', '25', '32', '40', '50', '63', '79', '100','126','159','200'})
yticks ([30	32	34	36	38	40	42	44	46	48	50	52	54	56	58	60	62	64	66	68	70	72	74	76	78	80	82	84	86	88	90])
yticklabels({'30', '', '', '', '', '40', '', '','','','50', '', '', '','','60', '', '', '','','70', '', '', '','','80', '', '', '','','90'})

%% PLOT Calibrated auralised Vs original sps
p0 = 2e5;
linS = {'-','--','-','--','-','--'};
markS = {'d','*','o','s','x','p'};
col = [{'b'},{'g'},{'m'}];
legentry = [{'Medisin ball orig'},{'Medisin ball aur'},{'Footstep with shoes orig'},{'Footstep with shoes aur'},{'Footstep without shoes orig'},{'Footstep without shoes aur'}];
ep = 1;
figure()
hold on
for c = 1:3
    plot(fvecout_orig{1},10*log10(F_orig_3oct{c,ep}./p0),'-','Marker',markS{c},'MarkerSize',7,'LineWidth',1,'Color',col{c})
    plot(fvecout_orig{1},10*log10(F_aur_3oct{c,ep}./p0),'--','Marker',markS{c},'MarkerSize',7,'LineWidth',1,'Color',col{c})

end
grid on
legend([legentry])
g = xlabel('Frequency [ Hz ]');
set(g,'FontSize',14)
g = ylabel('Amplitude dB [ 20 \mu pa. ]');
set(g,'FontSize',14)
g = title('Auralized vs original Sound Pressure in 1/3rd Octave Band');
set(g,'FontSize',16)
xlim([20 200])
xticks([20.0	25.1	31.6	39.8	50.1	63.1	79.4	100.0	125.9	158.5	199.5])
xticklabels({'20', '25.1', '31.6', '39.8', '50.1', '63.1', '79.4', '100','125.9','158.5','199.5'})
%


F_aur_3oct{2} F_aur_3oct{3}]))
% plot(fvecout{1},20*log10([F_aur_spec{1} F_aur_spec{2} F_aur_spec{3}])./p0)
legend('Medisin ball','Footstep with shoes','Footstep without shoes')
xlim([20 200])
xticks([20.0	25.1	31.6	39.8	50.1	63.1	79.4	100.0	125.9	158.5	199.5])
xticklabels({'20', '25.1', '31.6', '39.8', '50.1', '63.1', '79.4', '100','125.9','158.5','199.5'})
grid on
g = xlabel('Frequency [ Hz ]');
set(g,'FontSize',14)
g = ylabel('Amplitude dB [ uncal. ]');
set(g,'FontSize',14)
g = title('Auralized Sound Pressure in 1/3rd Octave Band');
set(g,'FontSize',16)

