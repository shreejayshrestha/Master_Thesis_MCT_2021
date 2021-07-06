%% Program to read and proccess acceleration signals from medisin ball
% Created Feb-May 2021
% Author: shreejay
% shreejayshrestha@gmail.com

%% Stay in the root folder: Medisinball_clean, to run this script 
clc
clear variables
%% Part1: Reading wave files with decay signals and loading decaystarttime matrix
tic
source_title = cell(1,4);
source_title{1,1} = ['Medisinball on Main Floor Ch1'];
source_title{1,2} = ['Medisinball on Main Floor Ch2'];
source_title{1,3} = ['Medisinball on Floating Floor Ch1'];
source_title{1,4} = ['Medisinball on Floating Floor Ch2'];

reload = 0; 
if reload == 1
    path_audiofiles = 'D:\NTNU\MCT2018_2020\MasterThesis2021\LABMeasurements\Clean_measurement_set_for_Analysis\Medisinball_clean\data_acceleration';
    audiofiles = dir(fullfile(path_audiofiles,'*.wav'));
    filenames = natsort({audiofiles.name}');
   
    % associating the decaysignals of main floor and floating floor
    rawacc_mf = cell(1,6);
    rawacc_ff = cell(1,6);
    for i = 1:numel(filenames)/2
        rawacc_mf{1,i} = audioread(fullfile(path_audiofiles,filenames{i,1}));
        rawacc_ff{1,i} = audioread(fullfile(path_audiofiles,filenames{i+6,1}));
    end
    
    % Decaystarttime matrix
    struct_decaytimes = load(strcat(path_audiofiles,'\acc_decaystarttimes_mb.mat'));
    decaystarttimes = struct_decaytimes.decaystarttimes_mfff;
       
    %% Acceleration signals Calibration
    % Main floor
    acc_scaled = cell(1,4);
    calibrationfactor = 10; % ref signal from VOC calibrator = 10mv RMS eqv 10m/s^2
    [acc_scaled{1,1},check_rms_precal_ch1_mf,check_rms_postcal_ch1_mf] = calibrate_acceleration(rawacc_mf{1,1},rawacc_mf{1,3},rawacc_mf{1,5},calibrationfactor); 
    [acc_scaled{1,2},check_rms_precal_ch2_mf,check_rms_postcal_ch2_mf] = calibrate_acceleration(rawacc_mf{1,2},rawacc_mf{1,4},rawacc_mf{1,6},calibrationfactor); 

    % Floating floor
    [acc_scaled{1,3},check_rms_precal_ch1_ff,check_rms_postcal_ch1_ff] = calibrate_acceleration(rawacc_ff{1,1},rawacc_ff{1,3},rawacc_ff{1,5},calibrationfactor); 
    [acc_scaled{1,4},check_rms_precal_ch2_ff,check_rms_postcal_ch2_ff] = calibrate_acceleration(rawacc_ff{1,2},rawacc_ff{1,4},rawacc_ff{1,6},calibrationfactor); 

    %% Calculating decay, delay and aligned signals in time domain
    %Note: The delay time between each signals will allow to identify the early
    %and the most lagging signals which are helpful to align all the signals.
    decay = cell(1,4);
    delay = cell(1,4);
    sigalignrefnum = cell(1,2);
    decay_aligned = cell(1,2);
    
    fs = 48000;
    decaylength = 3;
    
    % Decay Main floor 
%     function [decay,delay,sigalignrefnum,sig_align_full,sig_align_trim] = func_decaydelayalignsig(fs,longsigwithdecays,decaystarttime,decaylen)

    [~,~,~,~,decay{1,1}] = func_decaydelayalignsig(fs,acc_scaled{1,1},decaystarttimes{1,1},decaylength);
    [~,~,~,~,decay{1,2}] = func_decaydelayalignsig(fs,acc_scaled{1,2},decaystarttimes{1,1},decaylength);
       
    % Decay Floating floor
    decay{1,3} = func_chopdecay(fs,acc_scaled{1,3},decaystarttimes{1,2},decaylength);
    decay{1,4} = func_chopdecay(fs,acc_scaled{1,4},decaystarttimes{1,2},decaylength);
    
     % check decay
    indx = 2;
    ch1_mf = func_plot_allrep_allp(decay{1,indx},source_title{1,indx});
    
    noiseselecttvec = cell(1,4);
    noisesamp = cell(1,4);
 
    % Noise selection Main floor
    noiseselecttvec{1,1}    =(round(51.5*fs):round(53.5*fs));
    noisesamp{1,1}          = acc_scaled{1,1}(noiseselecttvec{1,1});
    noiseselecttvec{1,2}    =(round(51.5*fs):round(53.5*fs));
    noisesamp{1,2}          = acc_scaled{1,2}(noiseselecttvec{1,2});
    
    % Noise selection Floating floor
    noiseselecttvec{1,3}   =(round(53*fs):round(56*fs));
    noisesamp{1,3}         = acc_scaled{1,3}(noiseselecttvec{1,3});
    noiseselecttvec{1,4}   =(round(53*fs):round(56*fs));
    noisesamp{1,4}         = acc_scaled{1,4}(noiseselecttvec{1,4});
    
    %% Exporting 
%     save data_acceleration\acc_sig_aligned_mb.mat decay_aligned
%     save data_acceleration\chopped_decay_mb.mat decay
%     save data_acceleration\acc_noiseselection_mb.mat noisesamp
%     save data_acceleration\acc_decaystarttimes_mb.mat decaystarttimes
end 

%% Part 2: This part can be run without running Part 1.
path_audiofiles = 'D:\NTNU\MCT2018_2020\MasterThesis2021\LABMeasurements\Clean_measurement_set_for_Analysis\Medisinball_clean\data_acceleration';
% Reset/Load following parameters/matrices
fs = 48000;
% struct_decay = load(strcat(path_audiofiles,'\chopped_decay_mb_48000Hz.mat'));
struct_sigaligned = load(strcat(path_audiofiles,'\acc_sig_aligned_mb.mat'));
struct_noise = load(strcat(path_audiofiles,'\acc_noiseselection_mb.mat'));

%% Assign individual decay signals 
% acc_sig = struct_decay.decay;
acc_sig = struct_sigaligned.decay_aligned;
noise_sig = struct_noise.noisesamp;

%% Filtering wit Lp filter
%% Filtering at 200 Hz (previously tried with 100. While no filtering at all, signals are very noisy)
[Blp,Alp]= butter(4,200/(fs/2));
% acc_temp_f = filter(Blp,Alp,acc_temp);
acc_filt = cell(1,4);
for k = 1:4
    for i = 1:5
        for j = 1:10
            acc_filt{1,k}{i,j}= filter(Blp,Alp,acc_sig{1,k}{i,j});
        end
    end
end
% % check
% % ch1_mf = func_acc_plot_allrep_allp(acc_filt{1,1},source_title{1,1});
% % acc_filt_temp = acc_filt{1,1};
% % plot(acc_filt_temp{1,1})

%% Final Cut: cuttiing the acc decay
acc_finalcut = cell(1,4);
acc_cut_win = cell(1,4);
% acc_filt = acc_sig;
nwin = 70000;
hanwin = hanning(nwin*2);
halfhanwin = hanwin(nwin+1:end);
figure()
plot(halfhanwin)

%% Filtering at 100 Hz
% for k = 1:4
%     for i = 1:5
%         for j = 1:10
%             if k ==1 && i==5
%             [npeak,ninit] = findinit(acc_filt{1,k}{i,j}.^2,0.03);
%             else
%             [npeak,ninit] = findinit(acc_filt{1,k}{i,j}.^2,0.006);
%             end
%             acc_finalcut{1,k}{i,j}= [0;acc_filt{1,k}{i,j}(ninit:end)];
%             acc_cut_win{1,k}{i,j}  = acc_finalcut{1,k}{i,j};
%             acc_cut_win{1,k}{i,j}(end-nwin+1:end) = acc_cut_win{1,k}{i,j}(end-nwin+1:end).*halfhanwin;
%         end
%     end
% end

%% use this while filtering at 200 Hz

for k = 1:4
    for i = 1:5
        for j = 1:10
            if k ==1 && i==5
            [npeak,ninit] = findinit(acc_filt{1,k}{i,j}.^2,0.03);
            elseif k==1 && i == 3 && j==4
            [npeak,ninit] = findinit(acc_filt{1,k}{i,j}.^2,0.008);
            elseif k==1 && i == 4 && j==5
            [npeak,ninit] = findinit(acc_filt{1,k}{i,j}.^2,0.008);
            elseif k==2 && i == 5 && j==10
            [npeak,ninit] = findinit(acc_filt{1,k}{i,j}.^2,0.008);
            else
            [npeak,ninit] = findinit(acc_filt{1,k}{i,j}.^2,0.006);
            end
            acc_finalcut{1,k}{i,j}= [0;acc_filt{1,k}{i,j}(ninit:end)];
            acc_cut_win{1,k}{i,j}  = acc_finalcut{1,k}{i,j};
            acc_cut_win{1,k}{i,j}(end-nwin+1:end) = acc_cut_win{1,k}{i,j}(end-nwin+1:end).*halfhanwin;
        end
    end
end
%% Export acceleration of the floor in time domain 
acc_floor_mb = acc_cut_win;
% save Final_Output_from_allscripts_23_04_2021\medisinball\acc_floor_mb.mat acc_floor_mb % cell1 is mf and cel2 is ff.
% save Final_Output_from_allscripts_23_04_2021\medisinball\acc_filtered_at_200Hz\acc_floor_mb.mat acc_floor_mb % cell1 is mf and cel2 is ff.

%% checking Final cut
checkfinalcut = 0;
if checkfinalcut ==1
    ind =1; 
    p = 4;  % ind1,p3 r4. p4 rep 5, r7 to be excluded. ind2 only p5 r 10
    for k = 1:10
        subplot(2,5,k)
    %     plot(acc_sig{1,ind}{p,k})
        plot(acc_cut_win{1,ind}{p,k})
        title (['P\_',int2str(p),'rep\_',int2str(k)])
    %       xlim([0 5000])
    %     ylim([0 0.3])
    end
end
   
%% PLOT: 
plotthissection = 0;
if plotthissection ==1
    %% CHECK  all reps at all points for 
    % function tvec = func_acc_plot_allrep_allp(decaysig,sigtitle)
    ch1_mf = func_plot_allrep_allp(acc_cut_win{1,1},source_title{1,1});
    ch2_mf = func_plot_allrep_allp(acc_cut_win{1,2},source_title{1,2});
    ch1_ff = func_plot_allrep_allp(acc_cut_win{1,3},source_title{1,3});
    ch2_ff = func_plot_allrep_allp(acc_cut_win{1,4},source_title{1,4});

    %% Check all rep overlappiing at one selected point.
    % function tvec = func_acc_plot_allrep_onep(decaysig,selectp,sigtitle)
    p = 1; % select point
    ch1_mf = func_plot_allrep_onep(acc_filt{1,1},p,source_title{1,1});
    ch2_mf = func_plot_allrep_onep(acc_filt{1,2},p,source_title{1,2});
    ch1_ff = func_plot_allrep_onep(acc_filt{1,3},p,source_title{1,3});
    ch2_ff = func_plot_allrep_onep(acc_filt{1,4},p,source_title{1,4});
end

%% Frequency Domain
fs=48000;
nfft = cell(1,4);
ivf = cell (1,4);
fvec = cell (1,4);
% decaylength = zeros(5,10,4); % to define nfft
% 
% for i = 1:4
%     for j = 1:5
%         for k = 1:10
%             decaylength(j,k,i)= length(acc_cut_win{1,i}{j,k});
%         end
%     end
% end
% 
% maxdecaylength = max(max(max(decaylength)));

for i = 1:4
%     nfft{1,i} = 2^(nextpow2(maxdecaylength)+2);
%     nfft{1,i} = 8388608;% previous nfft of microphone
    nfft{1,i} = 2097152;% updated nfft of microphone 20.04.2021
    ivf{1,i} = [1:nfft{1,i}/20];
    fvec{1,i}= fs/nfft{1,i}*[0:nfft{1,i}/2-1];
end

Fbig_acc = cell(1,4);
Fbig_acc{1,1} = func_freqdom(acc_cut_win{1,1},nfft{1,1});
Fbig_acc{1,2} = func_freqdom(acc_cut_win{1,2},nfft{1,2});
Fbig_acc{1,3} = func_freqdom(acc_cut_win{1,3},nfft{1,3});
Fbig_acc{1,4} = func_freqdom(acc_cut_win{1,4},nfft{1,4});

%% CHECK  all reps at all points  
% function fvec = func_acc_plotfreq_allrep_allp(decaysig,noiseselection,nfft,sigtitle)
check_freqplot =0;
if check_freqplot ==1
    fvec_ch1mf = func_plotfreq_allrep_allp(Fbig_acc{1,1},nfft{1,1},source_title{1,1});
    fvec_ch2mf = func_plotfreq_allrep_allp(Fbig_acc{1,2},nfft{1,2},source_title{1,2});
    fvec_ch1ff = func_plotfreq_allrep_allp(Fbig_acc{1,3},nfft{1,3},source_title{1,3});
    fvec_ch2ff = func_plotfreq_allrep_allp(Fbig_acc{1,4},nfft{1,4},source_title{1,4});

    %% CHECK  all reps at one selected point
    p = 1;
    % function fvec = func_acc_plotfreq_allrep_onep(decaysig,nfft,sigtitle,selectp)
    fvec_ch1mf = func_plotfreq_allrep_onep(Fbig_acc{1,1},nfft{1,1},source_title{1,1},p);
    fvec_ch2mf = func_plotfreq_allrep_onep(Fbig_acc{1,2},nfft{1,2},source_title{1,2},p);
    fvec_ch1ff = func_plotfreq_allrep_onep(Fbig_acc{1,3},nfft{1,3},source_title{1,3},p);
    fvec_ch2ff = func_plotfreq_allrep_onep(Fbig_acc{1,4},nfft{1,4},source_title{1,4},p);
end

%% Averaging signals
Favgcomp = cell(1,4);
Favgreal_acc = cell(1,4);

for i = 1:4
    for j = 1:5
        Favgcomp{1,i}= squeeze(Fbig_acc{1,i}(:,:,j));
        if i==1 && j==4
           Favgcomp{1,i}(:,7)= [];
        end
        Favgreal_acc{1,i}(:,j)=  mean(abs(Favgcomp{1,i}).' ).'; 
    end
end

%% Export average acceleration spectrum in frequency domain
acc_Spectrum_Avg_mb = Favgreal_acc;
% save Final_Output_from_allscripts_23_04_2021\medisinball\acc_filtered_at_200Hz\acc_Spectrum_Avg_mb.mat acc_Spectrum_Avg_mb % cell1 is mf and cel2 is ff.

%% Plotting avg acc with observed and analytical natural freq of the conc floor
checkavg =1;
if checkavg==1
%     plot(Favgreal{1,1}(ivf{1,1},1))
    % plotavg = 0;
    % if plotavg ==1
    f0_analyt = [12.80 29.22 34.79 51.21 56.59 71.43 78.58 87.85 115.22];
    f0_observed = [29 38 49 56 70 79];
    
   figure()
   for i = 1:2
        subplot (1,2,i)
        h(1:5) = semilogx(fvec{1,i}(ivf{1,i}),20*log10(Favgreal_acc{1,i}(ivf{1,i},:)));
        c=6;
        for k = 1:length(f0_observed)
            h(c) = xline(f0_observed(k),'-.r');
            c=c+1;
        end
        dd=c+1;
         for k = 1:length(f0_analyt)
            h(dd) = xline(f0_analyt(k),'-.b');
            dd=dd+1;
         end
        xlim([20 100])
        ylim([30 75])
        xlabel('Frequency(Hz)')
        grid on
        title(['Avg of acc spectrum around P1: Mb main floor Ch:',int2str(i),' (using LP filter)']);
   end
        legend(h([1,2,3,4,5,6,21]),{'P1','P2','P3','P4','P5','f_0 measurement','f_0 analytical'});
        legend('boxoff')
        legend('Location','EastOutside')
        ylabel('dB re 1m/s^2')
end

%% Finding relation between acce. spektrum amplitude vs different positions
figure()
semilogx(fvec{1,1}(ivf{1,1}),20*log10(Favgreal_acc{1,1}(ivf{1,1},1)))
xlim([20 100])
xlabel('Frequency (Hz)')
ylabel('dB (re 1m/s^2)')
legend('P1','P2','P3','P4','P5');
title('Average of acceleration spectrum CH1 at each excitation Point: Mb main floor');

%% Finding Peaks, locs in each excitation point within 20-100 Hz
pks = cell(2,5);
lcs = cell(2,5);
Favgrealdb = cell(1,2);

for k = 1:2
    for j = 1:5
        Favgrealdb{1,k}(:,j) = 20*log10(Favgreal_acc{1,k}(:,j));
    end
end

for k = 1:2
        for j = 1:5
            target = 20*log10(Favgreal_acc{1,k}((round(20*nfft{1,1}/fs):round(100*nfft{1,1}/fs)),j));
             [p,l] = findpeaks(target,'MinPeakDistance',65);
%             [p,l] = findpeaks(target,'MinPeakDistance',1000); % old value 65
            pks{k,j} = p;
            lcs{k,j}= l;
            clear p l
        end
end

seg = cell(1,2);
ll = 20;
ul = 300;
% seg{1,1} = (round(ll*nfft{1,1}/fs):round(ul*nfft{1,1}/fs))';
% seg{1,2} = (round(ll*nfft{1,2}/fs):round(ul*nfft{1,2}/fs))';

seg{1,1} = (round(ll*nfft{1,1}/fs):round(ul*nfft{1,1}/fs))';
seg{1,2} = (round(ll*nfft{1,2}/fs):round(ul*nfft{1,2}/fs))';
figure()
ch = 1; % input 1 for channel 1 and 2 for channel 2
for j = 1:5
    subplot (3,2,j)
%     semilogx(fvec{1,ch}(ivf{1,ch}(lim{1,1})),20*log10(Favgreal{1,ch}(ivf{1,ch}(lim{1,1}),j)),fvec{1,ch}(ivf{1,ch}(lcs{1,j}+218-1)),pks{1,j},'or')
    semilogx(fvec{1,ch}(seg{1,ch}),20*log10(Favgreal_acc{1,ch}(seg{1,ch},j)),fvec{1,ch}(lcs{ch,j}+(seg{1,ch}(1))-1),pks{ch,j},'or')
    xlim([20 300])
    legend('Ch1 mf')
    title(['Average of acceleration spectrum for excitation Point1\_',int2str(j),' Mb main floor Ch:',int2str(ch)])
    grid on
%     pause
end

%% locating pks and lcs
lcslength = zeros(2,5);
for ch = 1:2
    for j = 1:5
        lcslength(ch,j)= length(lcs{ch,j});
    end
end

maxlcslength = max(max(lcslength));
lcs_Hz = cell(2,5);
lcs_Hz2 = cell(2,5);
lcs_adj = cell(2,5);
lcs_Hz_adj_acc = zeros(maxlcslength,5,2);
lcs_adj_acc = zeros(maxlcslength,5,2);
% 

for ch = 1:2
    for j = 1:5
       lcs_Hz{ch,j}= (lcs{ch,j}+(seg{1,ch}(1))-1)*(fs/nfft{1,ch});
       lcs_Hz2{ch,j} =  lcs_Hz{ch,j};
       lcs_adj{ch,j} = lcs{ch,j};
       
       if length(lcs_Hz2{ch,j}) <maxlcslength
            lcs_Hz2{ch,j}(end+1:maxlcslength)=0;
             lcs_adj{ch,j}(end+1:maxlcslength)=0;
%             lcs_Hz_adj(:,j,ch)=  round(lcs_Hz2{ch,j});
            lcs_Hz_adj_acc(:,j,ch)=  (lcs_Hz2{ch,j});
            lcs_adj_acc(:,j,ch)=  lcs_adj{ch,j};

       elseif length(lcs_Hz2{ch,j})== maxlcslength
%             lcs_Hz_adj(:,j,ch)=  round(lcs_Hz2{ch,j}); 
            lcs_Hz_adj_acc(:,j,ch)=  (lcs_Hz2{ch,j}); 
            lcs_adj_acc(:,j,ch)=  lcs_adj{ch,j};
       end
    end
end

%% Export peakindex
acc_Peakindex_Hz_mb = lcs_Hz_adj_acc;
% save Final_Output_from_allscripts_23_04_2021\medisinball\acc_filtered_at_200Hz\acc_Peakindex_Hz_mb.mat acc_Peakindex_Hz_mb % cell1 is mf and cel2 is ff.

% this matrix contains the index of peak of avg acc in Hz
%  save data_acceleration\peakindexin_Hz_accavg.mat lcs_Hz_adj_acc; %

%% locating f0 observed
% It is important to obtain amplitude of mic at each point at the same
% natural freqneucy. So define a fixed natural freq. for all points.

f0pk_Hzloc = zeros(5,7,2);
f0pk_indloc = zeros(5,7,2);
f0pk = zeros(5,7,2);
lookupf0 = [29 37 40 49 56 71 79]; % <28-29>,<37-38>, 
for j = 1:2
    for c = 1:7
        for r = 1:5
            [row,~] = find(lcs_Hz_adj_acc(:,r,j)>lookupf0(c) & lcs_Hz_adj_acc(:,r,j)<lookupf0(c)+1);
%             if  r ==1 && isempty(row)
%                 row = 0;
%             elseif isempty(row)
%                 row = f0pk_Hzloc(r-1,c,j);
%             end
            if  isempty(row)
               [~,closest_index] = min(abs(lcs_Hz_adj_acc(:,r,j)-lookupf0(c)));
               f0pk_Hzloc(r,c,j) = closest_index;
            else
                f0pk_Hzloc(r,c,j) = row;
            end
            f0pk_indloc(r,c,j) = lcs_adj_acc(f0pk_Hzloc(r,c,j),r,j);
%             
        end
    end
end

% [~,closest_index] = min(abs(lcs_Hz_adj_acc(:,3,1)-lookupf0(7)));

%% Defining fixed natural frequency
f0fix = zeros(4,7);
for ii = 1:7
    f0fix(2,ii) = median(f0pk_indloc(:,ii,1)); % freq index ch 1
    f0fix(1,ii) = fvec{1,1}(f0fix(2,ii)+(seg{1,1}(1))-1); % freq in Hz ch 1
    f0fix(4,ii) = median(f0pk_indloc(:,ii,2)); % freq index ch 2
    f0fix(3,ii) = fvec{1,1}(f0fix(4,ii)+(seg{1,1}(1))-1); % freq in Hz ch 2
end

%% Export fixed observed natural frequency of the floor
acc_floor_Observed_Resfreq_mb = f0fix;
% save Final_Output_from_allscripts_23_04_2021\medisinball\acc_filtered_at_200Hz\acc_floor_Observed_Resfreq_mb.mat acc_floor_Observed_Resfreq_mb % cell1 is mf and cel2 is ff.

%% creating matrix with amplitude for fixed f0 at all 5 points
%% This is the latest version for this part.

%% f0fix(4,:) is chosen as the f0 for the floor which is acc from ch 2
amp_acc_vs_pts = zeros(5,7,2);
for p = 1:5
    for ii = 1:7
        amp_acc_vs_pts(p,ii,1) = Favgrealdb{1,1}((f0fix(4,ii)+(seg{1,1}(1))-1),p); 
        amp_acc_vs_pts(p,ii,2) = Favgrealdb{1,2}((f0fix(4,ii)+(seg{1,1}(1))-1),p); 

    end
end
%% Export average amplitude at observed resonance frequency
acc_amp_at_Obs_Resfreq_vs_pts_mb = amp_acc_vs_pts ;
% save Final_Output_from_allscripts_23_04_2021\medisinball\acc_filtered_at_200Hz\acc_amp_at_Obs_Resfreq_vs_pts_mb.mat acc_amp_at_Obs_Resfreq_vs_pts_mb % cell1 is mf and cel2 is ff.
toc


