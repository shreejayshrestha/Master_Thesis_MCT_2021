%% Program to read raw microphone recordings of transmitted sound pressure in the receiving room
%% Sound source: Footstep with and without shoes on Main and Floating Floor

% Created Feb-May 2021
% Author: shreejay
% shreejayshrestha@gmail.com
clc
clear variables
fs = 48000;

source_title_time = cell(1,2);
source_title_freq = cell(1,2);
source_title_time{1,1} = ['Footstep with shoes on Main Floor'];
source_title_time{1,2} = ['Footstep with shoes on Floating Floor'];
source_title_freq{1,1} = ['Sound pressure spectrum: Footstep with shoes on Main Floor'];
source_title_freq{1,2} = ['Sound pressure spectrum: Footstep with shoes on Floating Floor'];

reload = 0;
if reload ==1
    path_audiofiles = 'D:\NTNU\MCT2018_2020\MasterThesis2021\LABMeasurements\Clean_measurement_set_for_Analysis\Footstep_clean\data_microphone\with_shoes';
    audiofiles = dir(fullfile(path_audiofiles,'*.wav'));
    filenames = natsort({audiofiles.name}');

    rawmic_mf = cell(1,3);
    rawmic_ff = cell(1,3);
    for ii = 1:numel(filenames)/2
        rawmic_mf{1,ii} = audioread(fullfile(path_audiofiles,filenames{ii,1}));
        rawmic_ff{1,ii} = audioread(fullfile(path_audiofiles,filenames{ii+3,1}));
    end
    %     
    %% Microphone Calibration
    [micscaled_mf,xprecal_db_mf,xpostcal_db_mf] = func_miccalibration(rawmic_mf{1,1},rawmic_mf{1,2},rawmic_mf{1,3});
    [micscaled_ff,xprecal_db_ff,xpostcal_db_ff] = func_miccalibration(rawmic_ff{1,1},rawmic_ff{1,2},rawmic_ff{1,3});

    micscaled = cell(1,2);
    micscaled{1,1} = micscaled_mf;
    micscaled{1,2} = micscaled_ff;

    %% Decaystarttime matrix
    decaystarttimes = cell(1,2);
    decaystarttimes{1,1} = [1.249	197.9	392.86	574.53	783.56
    17.18	213.9	408.86	595.43	799.96
    34.98	230.6	425	612.56	816.29
    51.28	262.4	441.3	628.98	832.51
    67.52	279.09	458.29	645.05	851.42
    84.08	295.33	474.42	661.3	867.54
    100.23	311.36	490.78	677.4	884.54
    116.12	327.36	507.61	693.35	900.66
    132.33	343.52	523.89	709.31	922.79
    148.59	359.54	540	725.45	939.57]; % mainfloor (defined by manually reading the signal in Reaper(DAW))

    decaystarttimes{1,2} = [0.15	242.36	465.73	700.6	943.15
    24.47	264.92	488.51	725.7	966.53
    47.29	290.42	511.68	748.16	992.34
    71.45	312.34	535.19	771.06	1015.22
    95.51	334.1	558.47	793.56	1038.99
    118.57	356.18	582.91	816.45	1062.21
    141.62	377.99	606.63	838.51	1084.303
    164.6	399.8	629.58	861.67	1106.51
    187.59	421.5	652.76	884.97	1128.51
    210.21	444.56	676.02	909.51	1150.68];% floating floor (defined by manually reading the signal in Reaper(DAW))

    %  save Footstep_clean\data_microphone\with_shoes\decaystarttimes_withshoes.mat decaystarttimes

    %% Calculating decay, delay and aligned signals in time domain
    %Note: The delay time between each signals will allow to identify the early
    %and the most lagging signals which are helpful to align all the signals.
    decay = cell(1,2);
    delay = cell(1,2);
    noisesamp = cell(1,2);

    sigalignrefnum = cell(1,2);
    decay_aligned = cell(1,2);
    fs = 48000;
    decaylength = 7;

    % function [decay,delay,sigalignrefnum,sig_align_full,sig_align_trim] = func_decaydelayalignsig(fs,longsigwithdecays,decaystarttime,decaylen)
    [decay{1,1},delay{1,1},sigalignrefnum{1,1},~,decay_aligned{1,1}] = func_decaydelayalignsig(fs,micscaled{1,1},decaystarttimes{1,1},decaylength);
    [decay{1,2},delay{1,2},sigalignrefnum{1,2},~,decay_aligned{1,2}] = func_decaydelayalignsig(fs,micscaled{1,2},decaystarttimes{1,2},decaylength);

    noiseselect_mf =(round(12.9*fs):round(15.8*fs));
    noisesamp{1,1} = micscaled{1,1}(noiseselect_mf);

    noiseselect_ff =(round((19.4)*fs):round((21.4)*fs));
    noisesamp{1,2} = micscaled{1,2}(noiseselect_ff); 

    %% Check decays if the chosen decaystarttimes are ok
    check_decaystarttimes = 0;
    if check_decaystarttimes==1
        %% CHECK  all reps at all points for 
            % function tvec = func_acc_plot_allrep_allp(decaysig,sigtitle)
        ch = 1; % 1= mainfloor 2= floating floor
        figure()
        mic_check_decay = func_acc_plot_allrep_allp(decay_aligned{1,ch},source_title_time{1,ch}); % Looks beautiful
               %% Check all rep overlappiing at one selected point before filtering and windowing
            % function tvec = func_acc_plot_allrep_onep(decaysig,selectp,sigtitle)
        figure()
        ch = 1; % 1= mainfloor 2= floating floor
        for p = 1:5 % select point
        %     mic_mf = func_acc_plot_allrep_onep(mic_cut_win{1,1},p,source_title{1,1});
            mic_decay_check = func_acc_plot_allrep_onep(decay_aligned{1,ch},p,source_title_time{1,ch});
            pause
        end
    end

        
    %% Export decay aligned
    decay_aligned_footstep_with_shoes = decay_aligned;
%     save Footstep_clean\data_microphone\with_shoes\decay_aligned_footstep_with_shoes.mat decay_aligned_footstep_with_shoes
end
 
load 'D:\NTNU\MCT2018_2020\MasterThesis2021\LABMeasurements\Clean_measurement_set_for_Analysis\Footstep_clean\data_microphone\with_shoes\decay_aligned_footstep_with_shoes.mat';

%% Filtering wit Lp filter
[Blp,Alp]= butter(4,300/(fs/2));
mic_filt = cell(1,2);
for k = 1:2
    for ii = 1:size(decay_aligned_footstep_with_shoes{1,k},1)
        for j = 1:size(decay_aligned_footstep_with_shoes{1,k},2)
            mic_filt{1,k}{ii,j}= filter(Blp,Alp,decay_aligned_footstep_with_shoes{1,k}{ii,j});
        end
    end
end

%% Using findit function to cut the signal and windowing to slightly fade the signal towards the end
mic_cut = cell(1,2);
mic_cut_win = cell(1,2);
nwin = 70000; % maximum possible window in this case is 23500, varies from signal to signal
hanwin = hanning(nwin*2);
halfhanwin = hanwin(nwin+1:end);
reltrig = cell(1,2);% this may have to be adjusted depending upon the signal. KEEP AN EYE ON IT. DO CHECK!!
reltrig{1,1} = [0.003 0.003 0.006 0.005 0.009]; % for main floor
reltrig{1,2} = [0.002 0.004 0.008 0.003 0.009]; % for floating floor

for k = 1:2
    for ii = 1:size(mic_filt{1,k},1)
        for j = 1:size(mic_filt{1,k},2)
             reltrigval = reltrig{1,k}(ii);
              if k ==1 && ii== 1 && j==9
                   reltrigval = 0.004;
               elseif  k ==1 && ii == 3 && j==7 
                   reltrigval = 0.07;
               elseif  k ==1 && ii == 5 && j==7
                  reltrigval = 0.01;
               elseif  k ==2 && ii == 2 && j==1
                  reltrigval = 0.006;
               elseif  k ==2 && ii == 3 && ismember(j,[1,2 3])
                  reltrigval = 0.02;
               elseif  k ==2 && ii == 5 && ismember(j,[1,2])
                  reltrigval = 0.02;
              end
            
           [npeak,ninit] = findinit(mic_filt{1,k}{ii,j}.^2,reltrigval);
            mic_cut{1,k}{ii,j}= [0;mic_filt{1,k}{ii,j}(ninit:end)];
            mic_cut_win{1,k}{ii,j}  = mic_cut{1,k}{ii,j};
            mic_cut_win{1,k}{ii,j}(end-nwin+1:end) = mic_cut_win{1,k}{ii,j}(end-nwin+1:end).*halfhanwin;
        end
    end
end   
%% PLOT: 
% check_signal_cut = 0;
% if check_signal_cut ==1
%% Check all rep overlappiing at one selected point after filtering and windowing    % function tvec = func_acc_plot_allrep_allp(decaysig,sigtitle)
%     mic_decay_check = func_plot_allrep_allp(mic_cut_win{1,1},source_title{1,1});% Looks good
%     mic_decay_check = func_plot_allrep_allp(mic_cut_win{1,2},source_title{1,2});
% function tvec = func_plot_allrep_onep(decaysig,selectp,sigtitle)
   
ch = 2; % 1= mainfloor 2= floating floor
    
for p = 3 % select point
%     mic_mf = func_plot_allrep_onep(mic_cut_win{1,1},p,source_title{1,1});
    mic_decay_check = func_plot_allrep_onep(mic_cut_win{1,ch},p,source_title_time{1,ch});
end
%% Export decay singals in time domain filtered with LP at 300 Hz
soundpressure_footstep_with_shoes = mic_cut_win;
% save Final_Output_from_allscripts_23_04_2021\footstep\soundpressure_footstep_with_shoes.mat soundpressure_footstep_with_shoes % cell1 is mf and cel2 is ff.
% save data_microphone\with_shoes\mic_cut_win_with_shoes.mat mic_cut_win
   
%% Frequency Domain
fs=48000;
nfft = cell(1,2);
ivf = cell (1,2);
fvec = cell (1,2);

% decaylength = zeros(5,10,2); % to define nfft. 
% Use this if its not desired to take nfft as 2097152 for consistency with medisin ball analysis 
% for i = 1:2
%     for j = 1:5
%         for k = 1:10
%             decaylength(j,k,i)= length(mic_cut_win{1,i}{j,k});
%         end
%     end
% end
% 
% maxdecaylength = max(max(max(decaylength)));

for ii = 1:2
%      maxdecaylength = length(mic_cut_win{{1,ii}{1});    
%      nfft{1,ii} = 2^(nextpow2(maxdecaylength)+4); % it is slightly higher. Lets check this again
    nfft{1,ii} = 2097152; % for consistency with medisin ball analysis
    ivf{1,ii} = [1:nfft{1,ii}/20];
    fvec{1,ii}= fs/nfft{1,ii}*[0:nfft{1,ii}/2-1];
end

Fbig_mic_with_shoes = cell(1,2);
Fbig_mic_with_shoes{1,1} = func_freqdom(mic_cut_win{1,1},nfft{1,1});
Fbig_mic_with_shoes{1,2} = func_freqdom(mic_cut_win{1,2},nfft{1,2});
  
%% CHECK  all reps at all points  
% function fig = func_plotfreq_allrep_allp(decaysig,nfft,sigtitle)
checkplot =0;
if checkplot ==1
    fig_mf = func_plotfreq_allrep_allp(Fbig_mic_with_shoes{1,1},nfft{1,1},source_title_freq{1,1});
    fig_ff = func_plotfreq_allrep_allp(Fbig_mic_with_shoes{1,2},nfft{1,2},source_title_freq{1,2});
 %% CHECK  all reps at one selected point
    p = 2;
    % function fig = func_plotfreq_allrep_onep(decaysig,nfft,sigtitle,selectp)
    fig_ff = func_plotfreq_allrep_onep(Fbig_mic_with_shoes{1,2},nfft{1,2},source_title_freq{1,2},p);
    
     %% CHECK  all reps at one selected point WITH PAUSE
    p = 1;
    ch = 1;
    % function fig = func_plotfreq_allrep_onep(decaysig,nfft,sigtitle,selectp)
    fig_ff = func_plotfreq_allrep_onep_withpause(Fbig_mic_with_shoes{1,ch},nfft{1,ch},source_title_freq{1,ch},p);
end

%% Averaging signals
Favgcomp = cell(1,2);% (1,1)= main floor and (1,2)= floating floor
Favgreal_mic = cell(1,2); % (1,1)= main floor and (1,2)= floating floor

for i = 1:2
    for j = 1:5
        Favgcomp{1,i}= squeeze(Fbig_mic_with_shoes{1,i}(:,:,j));
        if i==1 && j == 1
           Favgcomp{1,i}(:,2)= [];% removing P1 rep 2 main floor
        elseif i==1 && j == 3
            Favgcomp{1,i}(:,2)= [];% removing P3 rep 2 main floor
            Favgcomp{1,i}(:,4)= [];% removing P3 rep 5 main  floor
            Favgcomp{1,i}(:,4)= [];% removing P3 rep 6 main floor
        elseif i==2 && j == 1
            Favgcomp{1,i}(:,4)= [];% removing P1 rep 4 floating floor
        elseif i==2 && j == 2
            Favgcomp{1,i}(:,5)= []; % removing P2 rep 5 floating floor
            Favgcomp{1,i}(:,5)= []; % removing P2 rep 6 floating floor
            Favgcomp{1,i}(:,5)= []; % removing P2 rep 7 floating floor
            Favgcomp{1,i}(:,6)= []; % removing P2 rep 9 floating floor
            % the place of the deleted column is taken by following column
         elseif i==2 && j == 5
            Favgcomp{1,i}(:,3)= [];% removing P2 rep 3 floating floor
            Favgcomp{1,i}(:,3)= [];% removing P2 rep 4 floating floor
            Favgcomp{1,i}(:,7)= [];% removing P2 rep 9 floating floor
        end
        Favgreal_mic{1,i}(:,j)=  mean(abs(Favgcomp{1,i}).' ).'; 
    end
end

%% Export average sound pressure spectrum in frequency domain
SoundPressureAvg_footstep_with_shoes = Favgreal_mic;
% save Final_Output_from_allscripts_23_04_2021\footstep\SoundPressureAvg_footstep_with_shoes.mat SoundPressureAvg_footstep_with_shoes  % cell1 is mf and cel2 is ff.
checkavg = 0;
if checkavg ==1
% check avgfreq    
%% Plotting all 5 average in one plot
% fvec = fs/nfft{1,1}*[0:nfft{1,1}/2-1];
% ivf = [1:nfft{1,1}/20];
  ch = 2; % main floor
    figure()
    semilogx(fvec{1,1}(ivf{1,1}),20*log10((Favgreal_mic{1,ch}(ivf{1,1},:))))
    xlim([20 100])
%     xl = xline(55.54,'--r');
%     xl.LabelVerticalAlignment = 'bottom';
%     xl.LabelHorizontalAlignment = 'right';
%     for ii = 1:26
%         xl2 = xline(f0room_theo(ii),'--');
%     end
%     set(xl2,'linewidth',1);
    xlabel('Frequency (Hz)')
    ylabel('dB (arbritrary)')
    legend('P1','P2','P3','P4','P5','Natural freq.floor','Natural freq.room','Location','NorthWest');
    title([source_title_freq{1,ch},' (Average at each excitation Point)']);
    grid on
end

%% Finding peaks and location of the average
pks = cell(1,2);
lcs = cell(1,2);
Favgrealdb = cell(1,2);
for k = 1:2
    for j = 1:5
        Favgrealdb{1,k}(:,j) = 20*log10(Favgreal_mic{1,k}(:,j));
    end
end

for k = 1:2
        for j = 1:5
            target = 20*log10(Favgreal_mic{1,k}((round(20*nfft{1,k}/fs):round(100*nfft{1,k}/fs)),j));
%             [p,l] = findpeaks(target,'MinPeakDistance',70);
%             [p,l] = findpeaks(target,'MinPeakDistance',300); % old value 65
            [p,l] = findpeaks(target,'MinPeakDistance',70); % This works good for nfft = 2097152.00

            pks{k,j} = p;
            lcs{k,j}= l;
            clear p l
        end
end
seg = cell(1,2);
seg{1,1} = (round(20*nfft{1,1}/fs):round(100*nfft{1,1}/fs))';
seg{1,2} = (round(20*nfft{1,2}/fs):round(100*nfft{1,2}/fs))';

figure()
ch = 1; % input 1 for channel 1 and 2 for channel 2
for j = 1:5
    subplot (3,2,j)
%     semilogx(fvec{1,ch}(ivf{1,ch}(lim{1,1})),20*log10(Favgreal{1,ch}(ivf{1,ch}(lim{1,1}),j)),fvec{1,ch}(ivf{1,ch}(lcs{1,j}+218-1)),pks{1,j},'or')
    semilogx(fvec{1,ch}(seg{1,ch}),20*log10(Favgreal_mic{1,ch}(seg{1,ch},j)),fvec{1,ch}(lcs{ch,j}+(seg{1,ch}(1))-1),pks{ch,j},'or')
%     for ii = 1:26
%          xl2 = xline(f0room_theo(ii),'--');
%     end
    xlim([20 100])
    legend('Ch1 mf')
    title(['Average of Measured SPL spectrum for excitation Point1\_',int2str(j),' Mb main floor',int2str(ch)])
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
lcs_Hz_adj_mic = zeros(maxlcslength,5,2);
lcs_adj_mic = zeros(maxlcslength,5,2);
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
            lcs_Hz_adj_mic(:,j,ch)=  (lcs_Hz2{ch,j});
            lcs_adj_mic(:,j,ch)=  lcs_adj{ch,j};

       elseif length(lcs_Hz2{ch,j})== maxlcslength
%             lcs_Hz_adj(:,j,ch)=  round(lcs_Hz2{ch,j}); 
            lcs_Hz_adj_mic(:,j,ch)=  (lcs_Hz2{ch,j}); 
            lcs_adj_mic(:,j,ch)=  lcs_adj{ch,j};
       end
    end
end

% maxlcslength = max(max(lcslength));
% lcs_Hz = cell(2,5);
% lcs_Hz2 = cell(2,5);
% lcs_Hz_adj_mic = zeros(maxlcslength,5,2);
% lcs_ajd = cell(1,2);
% 
% for ch = 1:2
%     for j = 1:5
%        lcs_ajd {1,ch}(:,j)= lcs_ajd {1,ch}(:,j);
%     end
% end
% save data_microphone\peakindexin_Hz_micavg.mat lcs_Hz_adj_mic

%% locating f0 observed
% It is important to obtain amplitude of mic at each point at the same
% natural freqneucy. So define a fixed natural freq. for all points.

f0pk_Hzloc = zeros(5,8,2);
f0pk_indloc = zeros(5,8,2);
f0pk = zeros(5,8,2);
lookupf0 = [29 40 51 54 67 71 79 82];

for j = 1:2
    for c = 1:8
        for r = 1:5
            [row,~] = find(lcs_Hz_adj_mic(:,r,j)>lookupf0(c) & lcs_Hz_adj_mic(:,r,j)<lookupf0(c)+1);
            if  isempty(row)
               [~,closest_index] = min(abs(lcs_Hz_adj_mic(:,r,j)-lookupf0(c)));
               f0pk_Hzloc(r,c,j) = closest_index;
            else
                f0pk_Hzloc(r,c,j) = row;
            end
            f0pk_indloc(r,c,j) = lcs_adj_mic(f0pk_Hzloc(r,c,j),r,j);
        end
    end
end

% for j = 1:2
%     for c = 1:8
%         for r = 1:5
%             [row,~] = find(lcs_Hz_adj_mic(:,r,j)>lookupf0(c) & lcs_Hz_adj_mic(:,r,j)<lookupf0(c)+1);
%             if isempty(row)
%                 row = f0pk_Hzloc(r-1,c,j);
%             end
%             f0pk_Hzloc(r,c,j) = row;
%             f0pk_indloc(r,c,j) = lcs_adj_mic(f0pk_Hzloc(r,c,j),r,j);
%         end
%     end
% end
%% Defining fixed natural frequency
f0fix = zeros(2,8);
for j = 1:2
    for ii = 1:8
        f0fix(2,ii) = median(f0pk_indloc(:,ii,1));% freq index
        f0fix(1,ii) = fvec{1,1}(f0fix(2,ii)+(seg{1,1}(1))-1); % freq in Hz.
    end
end

room_Obs_Resfreq_fstep_with_shoes = f0fix;

%% Export the observed resonance frequency of the sound pressure from the measurements
% save Final_Output_from_allscripts_23_04_2021\footstep\room_Obs_Resfreq_fstep_with_shoes.mat room_Obs_Resfreq_fstep_with_shoes % cell1 is mf and cel2 is ff.
 
%% Creating matrix with sound pressure amplitude at the observed resonance frequency at all 5 excitation points

amp_mic_vs_pts = zeros(5,8);
for p = 1:5
    for ii = 1:8
        amp_mic_vs_pts(p,ii) = Favgrealdb{1,1}((f0fix(2,ii)+(seg{1,1}(1))-1),p); 
    end
end

soundpres_amp_at_Obs_Resfreq_vs_pts_ftep_with_shoes = amp_mic_vs_pts ;
%% Export
% save Final_Output_from_allscripts_23_04_2021\footstep\soundpres_amp_at_Obs_Resfreq_vs_pts_ftep_with_shoes.mat soundpres_amp_at_Obs_Resfreq_vs_pts_ftep_with_shoes% cell1 is mf and cel2 is ff.
