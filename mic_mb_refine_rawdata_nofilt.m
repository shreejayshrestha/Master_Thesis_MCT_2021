%% Program to read and process microphone signals for medisin ball 
%% This version is revised with reduced nfft size from 8388608 to 2097152
% revesion is also done for mic amplitude vs excitation points.

% Created Feb-May 2021
% Author: shreejay
% shreejayshrestha@gmail.com

clc
clear variables
tic
fs = 48000;
source_title_time = cell(1,2);
source_title_freq = cell(1,2);
source_title_time{1,1} = ['Medisinball on Main Floor'];
source_title_time{1,2} = ['Medisinball on Floating Floor'];
source_title_freq{1,1} = ['Sound pressure spectrum: Medisinball on Main Floor'];
source_title_freq{1,2} = ['Sound pressure spectrum: Medisinball on Floating Floor'];

reload =0 ;
if reload == 1
    path_audiofiles = 'D:\NTNU\MCT2018_2020\MasterThesis2021\LABMeasurements\Clean_measurement_set_for_Analysis\Medisinball_clean\data_microphone';
    audiofiles = dir(fullfile(path_audiofiles,'*.wav'));
    filenames = natsort({audiofiles.name}');

    rawmic_mf = cell(1,3);
    rawmic_ff = cell(1,3);
    for ii = 1:numel(filenames)/2
        rawmic_mf{1,ii} = audioread(fullfile(path_audiofiles,filenames{ii,1}));
        rawmic_ff{1,ii} = audioread(fullfile(path_audiofiles,filenames{ii+3,1}));
    end

    %% Decaystarttime matrix
    struct_decaytimes = load(strcat(path_audiofiles,'\mic_decaystarttimes_mb.mat'));
    decaystarttimes = struct_decaytimes.decaystarttimes_mfff;

    %% Microphone Calibration
    [micscaled_mf,xprecal_db_mf,xpostcal_db_mf] = func_miccalibration(rawmic_mf{1,1},rawmic_mf{1,2},rawmic_mf{1,3});
    [micscaled_ff,xprecal_db_ff,xpostcal_db_ff] = func_miccalibration(rawmic_ff{1,1},rawmic_ff{1,2},rawmic_ff{1,3});

    micscaled = cell(1,2);
    micscaled{1,1} = micscaled_mf;
    micscaled{1,2} = micscaled_ff;

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

    noiseselect_mf =(round(63.03*fs):round(69.65*fs));
    noisesamp{1,1} = micscaled{1,1}(noiseselect_mf);

    noiseselect_ff =(round((15*60+48)*fs):round((15*60+54)*fs));
    noisesamp{1,2} = micscaled{1,2}(noiseselect_ff); 
    
    %% Export decay aligned
    decay_aligned_mb = decay_aligned;
%     save Medisinball_clean\data_microphone\decay_aligned_mb.mat decay_aligned_mb
end

load 'D:\NTNU\MCT2018_2020\MasterThesis2021\LABMeasurements\Clean_measurement_set_for_Analysis\Medisinball_clean\data_microphone\decay_aligned_mb.mat';

%% Filtering wit Lp filter
%% Not Filtering this time

% [Blp,Alp]= butter(4,100/(fs/2));
% mic_filt = cell(1,2);
% for k = 1:2
%     for ii = 1:5
%         for j = 1:10
%             mic_filt{1,k}{ii,j}= filter(Blp,Alp,decay_aligned_mb{1,k}{ii,j});
%         end
%     end
% end

%% Using findit function to cut the signal and windowing to slightly fade the signal towards the end
mic_finalcut = cell(1,2);
mic_cut_win = cell(1,2);
nwin = 70000;
hanwin = hanning(nwin*2);
halfhanwin = hanwin(nwin+1:end);
reltrig = 0.001;

mic_filt = decay_aligned_mb;

for k = 1:2
    for ii = 1:5
        for j = 1:10
             if k ==2 && ii==5
                 reltrig = 0.003;
                 [npeak,ninit] = findinit(mic_filt{1,k}{ii,j}.^2,reltrig);
             else
             [npeak,ninit] = findinit(mic_filt{1,k}{ii,j}.^2,reltrig);
             end
            mic_finalcut{1,k}{ii,j}= [0;mic_filt{1,k}{ii,j}(ninit:end)];
            mic_cut_win{1,k}{ii,j}  = mic_finalcut{1,k}{ii,j};
            mic_cut_win{1,k}{ii,j}(end-nwin+1:end) = mic_cut_win{1,k}{ii,j}(end-nwin+1:end).*halfhanwin;
        end
    end
end   

soundpressure_mb = mic_cut_win;

%% Export total sound pressure signals in time domain
%  save Final_Output_from_allscripts_23_04_2021\medisinball\signals_not_filtered\soundpressure_mb.mat soundpressure_mb % cell1 is mf and cel2 is ff.

%% Checking  
plotthissection = 0;
if plotthissection ==1
    %% CHECK  all reps at all points for 
    % function tvec = func_acc_plot_allrep_allp(decaysig,sigtitle)
    mic_mf = func_plot_allrep_allp(mic_cut_win{1,1},source_title_time{1,1});
    mic_ff = func_plot_allrep_allp(mic_cut_win{1,2},source_title_time{1,2});
    %% Check all rep overlappiing at one selected point.
    % function tvec = func_acc_plot_allrep_onep(decaysig,selectp,sigtitle)
    figure()
    p = 5; % select point
%     mic_mf = func_acc_plot_allrep_onep(mic_cut_win{1,1},p,source_title{1,1});
    mic_ff = func_acc_plot_allrep_onep(mic_cut_win{1,2},p,source_title{1,2});
end

%% Frequency Domain
fs=48000;
nfft = cell(1,2);
ivf = cell (1,2);
fvec = cell (1,2);

for ii = 1:2
    maxdecaylength = length(mic_cut_win{1,ii}{1});    
%     nfft{1,ii} = 2^(nextpow2(maxdecaylength)+4); % it is slightly more
    nfft{1,ii} = 2^(nextpow2(maxdecaylength)+2);
    ivf{1,ii} = [1:nfft{1,ii}/20];
    fvec{1,ii}= fs/nfft{1,ii}*[0:nfft{1,ii}/2-1];
end

fvec_2097152_nfft_mb = fvec;

%% Exporting frequency vector
% save Final_Output_from_allscripts_23_04_2021\medisinball\fvec_2097152_nfft_mb.mat fvec_2097152_nfft_mb % cell1 is mf and cel2 is ff.

Fbig_mic = cell(1,2);
Fbig_mic{1,1} = func_freqdom(mic_cut_win{1,1},nfft{1,1});
Fbig_mic{1,2} = func_freqdom(mic_cut_win{1,2},nfft{1,2});
   
%% CHECK  all reps at all points  
% function fvec = func_acc_plotfreq_allrep_allp(decaysig,noiseselection,nfft,sigtitle)
checkplot =0;
if checkplot ==1
    fvec_mf = func_plotfreq_allrep_allp(Fbig_mic{1,1},nfft{1,1},source_title_freq{1,1});
    fvec_ff = func_plotfreq_allrep_allp(Fbig_mic{1,2},nfft{1,2},source_title_freq{1,2});
 %% CHECK  all reps at one selected point
    p = 1;
    % function fvec = func_acc_plotfreq_allrep_onep(decaysig,nfft,sigtitle,selectp)
    fvec_mf = func_plotfreq_allrep_onep(Fbig_mic{1,1},nfft{1,1},source_title_freq{1,1},p);
    fvec_ff = func_plotfreq_allrep_onep(Fbig_mic{1,2},nfft{1,2},source_title_freq{1,2},p);
end

%% Averaging signals
Favgcomp = cell(1,2);
Favgreal_mic = cell(1,2);

for ii = 1:2
    for j = 1:5
        Favgcomp{1,ii}= squeeze(Fbig_mic{1,ii}(:,:,j));
        Favgreal_mic{1,ii}(:,j)=  mean(abs(Favgcomp{1,ii}).' ).'; 
    end
end

SoundPressureAvg_mb = Favgreal_mic;

%% Export average sound pressure spectrum in frequency domain
%  save Final_Output_from_allscripts_23_04_2021\medisinball\signals_not_filtered\SoundPressureAvg_mb.mat SoundPressureAvg_mb % cell1 is mf and cel2 is ff.

%% Plotting all 5 average in one plot
% fvec = fs/nfft{1,1}*[0:nfft{1,1}/2-1];
% ivf = [1:nfft{1,1}/20];
p0 = 2e-5;
ch = 2;
plotavg = 0;
    figure()
    semilogx(fvec{1,ch}(ivf{1,ch}),20*log10((Favgreal_mic{1,ch}(ivf{1,ch},:))))
%     xlim([20 100])
%     ylim([10 95])

%     xl = xline(55.54,'--r');
%     xl.LabelVerticalAlignment = 'bottom';
%     xl.LabelHorizontalAlignment = 'right';
%     for ii = 1:26
%         xl2 = xline(f0room_theo(ii),'--');
%     end
%     set(xl2,'linewidth',1);
    g = xlabel('Frequency [ Hz ]');
    set(g,'FontSize',14);
    g = ylabel('Amplitude dB [ uncal. ]');
    set(g,'FontSize',14);
    g = legend('P1','P2','P3','P4','P5','Natural freq.floor','Natural freq.room','Location','NorthWest');
    set(g,'FontSize',12); 
    legend show
    legend('boxoff')
    legend('Location','SouthOutside','Orientation','Horizontal')
%     title('Sound pressure spectrum for exciat each excitation Point: Mb main floor');
    ax = gca;
    ax.FontSize = 12; 
%     g = subtitle();
%     set(g,'FontSize',14);
    grid on
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
ch = 2; % input 1 for channel 1 and 2 for channel 2
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

% f0pk_Hzloc = zeros(5,8,2);
% f0pk_indloc = zeros(5,8,2);
% f0pk = zeros(5,8,2);
% lookupf0 = [29 40 51 54 67 71 79 82]; 
% %  lookupf0 = [29 37 40 49 56 70 79]; % <28-29>,<37-38>, This is for
% %  acceleration
% for j = 1:2
%     for c = 1:8
%         for r = 1:5
%             [row,~] = find(lcs_Hz_adj_mic(:,r,j)>lookupf0(c) & lcs_Hz_adj_mic(:,r,j)<lookupf0(c)+1);
%             if  isempty(row)
%                [~,closest_index] = min(abs(lcs_Hz_adj_mic(:,r,j)-lookupf0(c)));
%                f0pk_Hzloc(r,c,j) = closest_index;
%             else
%                 f0pk_Hzloc(r,c,j) = row;
%             end
%             f0pk_indloc(r,c,j) = lcs_adj_mic(f0pk_Hzloc(r,c,j),r,j);
%         end
%     end
% end

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

%% locating f0 observed
% It is important to obtain amplitude of mic at each point at the same
% natural freqneucy. So define a fixed natural freq. for all points.

f0pk_Hzloc_full = zeros(5,25,2);
f0pk_indloc_full = zeros(5,25,2);
f0pk_full = zeros(5,25,2);
lookupf0_full = resfreqs_analytical(1:25); 
% valvec = lookupf0_full ./(fs/nfft{1});
% valvec = zeros(25,1);
% for i = 1:25
%     valvec(i) = find(fvec,lookupf0_full(i));
% end


% lookupf0_full = [29
% 35
% 40
% 46
% 50
% 53
% 58
% 61
% 68
% 70
% 71
% 76
% 79
% 81
% 81
% 86
% 86
% 88
% 88
% 91
% 93
% 94
% 96
% 100
% 100];
%  lookupf0 = [29 37 40 49 56 70 79]; % <28-29>,<37-38>, This is for
%  acceleration
for j = 1:2
    for c = 1:25
        for r = 1:5
            [row,~] = find(lcs_Hz_adj_mic(:,r,j)>lookupf0_full(c) & lcs_Hz_adj_mic(:,r,j)<lookupf0_full(c)+1);
            if  isempty(row)
               [~,closest_index] = min(abs(lcs_Hz_adj_mic(:,r,j)-lookupf0_full(c)));
               f0pk_Hzloc_full(r,c,j) = closest_index;
            else
                f0pk_Hzloc_full(r,c,j) = row;
            end
            f0pk_indloc_full(r,c,j) = lcs_adj_mic(f0pk_Hzloc_full(r,c,j),r,j);
        end
    end
end
%% Defining fixed natural frequency
% f0fix = zeros(2,8);
% for j = 1:2
%     for ii = 1:8
%         f0fix(2,ii) = median(f0pk_indloc_full(:,ii,1));% freq index
%         f0fix(1,ii) = fvec{1,1}(f0fix(2,ii)+(seg{1,1}(1))-1); % freq in Hz.
%     end
% end
% 
% room_Observed_Resfreq_mb = f0fix;
%% Defining fixed natural frequeny: full
f0fix_full = zeros(2,25);
for j = 1:2
    for ii = 1:25
        f0fix_full(2,ii) = median(f0pk_indloc_full(:,ii,1));% freq index
        f0fix_full(1,ii) = fvec{1,1}(f0fix_full(2,ii)+(seg{1,1}(1))-1); % freq in Hz.
    end
end
resfreq_mb_obs = cell(1,2);
resfreq_mb_obs{1} = unique(f0fix_full(1,:));
resfreq_mb_obs{2} = unique(f0fix_full(2,:));

room_Observed_Resfreq_mb_full = f0fix_full;


%% Export the observed resonance frequency of the sound pressure from the measurements
%  save Final_Output_from_allscripts_23_04_2021\medisinball\signals_not_filtered\room_Observed_Resfreq_mb.mat room_Observed_Resfreq_mb % cell1 is mf and cel2 is ff.
 
%% Creating matrix with sound pressure amplitude at the observed resonance frequency at all 5 excitation points
%% with floor natura freq
load('D:\NTNU\MCT2018_2020\MasterThesis2021\LABMeasurements\Clean_measurement_set_for_Analysis\Final_Output_from_allscripts_08_05_2021\acc\floor_f0fix_mb.mat')

%% with room natura freq

% amp_mic_vs_pts_f0floor = zeros(5,6);
amp_mic_vs_pts_f0room = zeros(5,8,2);

% for p = 1:5
%     for ii = 1:6
%         amp_mic_vs_pts_f0floor(p,ii) = Favgrealdb{1,1}((floor_f0fix_mb(4,ii)+(seg{1,1}(1))-1),p); 
%     end
% end
for ch = 1:2
    for p = 1:5
        for ii = 1:8
            amp_mic_vs_pts_f0room(p,ii,ch) = Favgrealdb{1,ch}((f0fix_full(2,ii)+(seg{1,ch}(1))-1),p); 
        end
    end
end
% for p = 1:5
%     for ii = 1:8
%         amp_mic_vs_pts(p,ii) = Favgrealdb{1,1}((f0fix(2,ii)+(seg{1,1}(1))-1),p); 
%     end
% end

%% Export
%% updated 09-05-2021
% save Final_Output_from_allscripts_23_04_2021\medisinball\signals_not_filtered\amp_mic_vs_pts_f0room.mat amp_mic_vs_pts_f0room % cell1 is mf and cel2 is ff.
% save Final_Output_from_allscripts_23_04_2021\medisinball\signals_not_filtered\amp_mic_vs_pts_f0floor.mat amp_mic_vs_pts_f0floor % cell1 is mf and cel2 is ff.

toc