%% Program to calculate analytical sound pressure in the receiving room 
% while using modal sum theory applying the function modalsum_rigidroom.m
% by Peter Svennson

% Created Feb-May 2021
% Author: shreejay
% shreejayshrestha@gmail.com
clc
clear variables
%% Setting up input values
source = zeros(5,3);
source(1,:)= [2.44 2.94 4.25]; %P1
source(2,:)= [1.92 2.94 4.25]; %P2
source(3,:)= [1.4  2.94 4.25]; %P3
source(4,:)= [1.92 3.56 4.25]; %P4
source(5,:)= [1.4  4.19 4.25]; %P5

receiver = [0.2 0.24 0.25]; % corner mic position
roomsize = [4.88 5.87 4.25];


%% Calculating shortest distance between source and receiver position
base_PR = [3.5 3.2 2.95 2.69 1.88];
height_PR = 4.25 -0.25 ;
short_dist_PR = zeros (1,5);
for i = 1:5
    short_dist_PR(i) = sqrt((base_PR(i))^2 + (height_PR)^2);
end
%%
T60 = [5.2 2.5 1];

safetyfactor = 1.1;
rhoair = 1.2;
cair = 344;
nfft = 2097152;% corresponding to the medisinball microphone nfft
fs = 48000;
fvec= fs/nfft*[0:nfft/2-1];
seg = (round(20*nfft/fs):round(120*nfft/fs))';

freqvec = fvec(seg);

room_analytical_soundpres = cell(3,5); %% cell(a,b) where a1 = T60_5sec, a2= T60_2.5sec, a3= T60_1sec andb = excitation points P1:P5
nmodes_modsum = cell(3,5);
room_analytical_resfreq = cell(3,5);
room_analytical_resfreq_sorted = cell(3,5);
modeamp_modsum = cell(3,5);

%% Using the modal sum function from Peter. S.
for t = 1:3
    for ii = 1:5
        [room_analytical_soundpres{t,ii},nmodes_modsum{t,ii},room_analytical_resfreq{t,ii},modeamp_modsum{t,ii}] = modalsum_rigidroom(source(ii,:),receiver,freqvec,roomsize,T60(t),safetyfactor,rhoair,cair);
        room_analytical_resfreq_sorted{t,ii} = sort(room_analytical_resfreq{t,ii});
        %     pressureadj = pressure;
    %     pressureadj(isnan(pressureadj)) = 0;
    end
end
 room_analytical_resfreq_cut = room_analytical_resfreq_sorted{1,1}(2:26);

%% Export
% save Final_Output_from_allscripts_23_04_2021\room_modal_analysis_modsum\room_analytical_soundpres.mat room_analytical_soundpres;
% save Final_Output_from_allscripts_23_04_2021\room_modal_analysis_modsum\nmodes_modsum.mat nmodes_modsum;
% save Final_Output_from_allscripts_23_04_2021\room_modal_analysis_modsum\room_analytical_resfreq.mat room_analytical_resfreq;
% save Final_Output_from_allscripts_23_04_2021\room_modal_analysis_modsum\room_analytical_resfreq_sorted.mat room_analytical_resfreq_sorted;
% save Final_Output_from_allscripts_23_04_2021\room_modal_analysis_modsum\room_analytical_resfreq_cut.mat room_analytical_resfreq_cut;
% save Final_Output_from_allscripts_23_04_2021\room_modal_analysis_modsum\modeamp_modsum.mat modeamp_modsum;

%% Check Plots
check = 0;
if check ==1
    prs = 20*log10(abs(room_analytical_soundpres{1,1}));
    figure(2)
    semilogx(fvec(seg),prs);
end

check2 = 0;
if check2 ==1
    RT= 1; % select 1,2 or3
    figure()
    for ii = 1:5
        y1 = 20*log10(abs(room_analytical_soundpres{1,ii}));
        y2 = 20*log10(abs(room_analytical_soundpres{2,ii}));
        y3 = 20*log10(abs(room_analytical_soundpres{3,ii}));
        subplot (3,2,ii)
        semilogx(fvec(seg),y1,fvec(seg),y2,fvec(seg),y3)
        for j = 1:24
            xl2 = xline(resfreqs(j),'--');
        end
        xlim([20 100])
        ylim([0 70])
    %     set(xl2,'linewidth',1)
%         ylabel('dB (arbritrary)')
%         xlabel('Frequency (Hz)')
        g = title(['Sound source @ Point :P',int2str(ii)]);
        set(g,'FontSize',12);
        grid on
        ax = gca;
        ax.FontSize = 12;
    end
  g = legend( 'RT 5sec','RT 2.5sec','RT 1sec','f\_analytical');
% g = legend(h([1,2,3,18]),{'RT 5sec','RT 2.5sec','RT 1sec','f_0 analytical'});
 set(g,'FontSize',12);
legend('boxoff')
legend('Location','SouthOutside')
end

%% Calculating room modes
% [resfreqs,modenumbers] = calcroommodefreqs(roomsize,maxmodenumber,cair);
 maxmodenumber = [30,30,30];
[resfreqs,room_modes_lydrom3] = calcroommodefreqs(roomsize,maxmodenumber,cair);
%% Rough
% Plotting mode amplitude
% figure()
% plot(20*log10(abs(modeamp_modsum{1,5}(1,:))))

%% Converting complex into real dB values
FrealdB = cell(1,3);

for ii = 1:3
    for j = 1:5
        FrealdB{ii}(:,j)=  20*log10(abs(room_analytical_soundpres{ii,j})); 
    end
end

modsum_SP_dB = FrealdB ;
% save Final_Output_from_allscripts_23_04_2021\room_modal_analysis_modsum\modsum_SP_dB.mat modsum_SP_dB;


xx1 = 20*log10(abs(room_analytical_soundpres{1,1}(1,:)));
xx11 = 20*log10(abs(room_analytical_soundpres{1,1}));
xx2 = 20*log10(abs(room_analytical_soundpres{1,2}(1,:)));


nfft = 2097152;
fvec = fs/nfft*[0:nfft/2-1];
ivf = find(fvec>= 20 & fvec<=100);
RT = 1;
figure()
    semilogx(fvec(seg),FrealdB{1}(:,3))
    xlim([20 100])
    ylim([10 95])
%% Finding peaks and location 
pks = cell(1,3);
lcs = cell(1,3);
for k = 1:3
        for j = 1:5
            target = FrealdB{1,k}(:,j);
            [p,l] = findpeaks(target,'MinPeakDistance',70); % Working good for all now
% %             [p,l] = findpeaks(target,'MinPeakDistance',300); % old value 65
%             if k==1            
%                 [p,l] = findpeaks(target,'MinPeakDistance',70); 
%             else
%                 [p,l] = findpeaks(target,'MinPeakDistance',70); 
%             end
            pks{k,j} = p;
            lcs{k,j}= l;
            clear p l
        end
end
figure()
RT = 1;
for j = 1:5
    subplot (3,2,j)
%     semilogx(fvec{1,ch}(ivf{1,ch}(lim{1,1})),20*log10(Favgreal{1,ch}(ivf{1,ch}(lim{1,1}),j)),fvec{1,ch}(ivf{1,ch}(lcs{1,j}+218-1)),pks{1,j},'or')
    semilogx(fvec(seg),(FrealdB{RT}(:,j)),fvec(lcs{RT,j}+(seg(1)-1)),pks{RT,j},'or')
%     for ii = 1:26
%          xl2 = xline(f0room_theo(ii),'--');
%     end
    xlim([20 100])
%     legend('Ch1 mf')
    title(['Average of Measured SPL spectrum for excitation Point1\_',int2str(j),' Mb main floor',int2str(RT)])
    grid on
%     pause
end
%% locating pks and lcs
lcslength = zeros(3,5);
for ch = 1:3
    for j = 1:5
        lcslength(ch,j)= length(lcs{ch,j});
    end
end

maxlcslength = max(max(lcslength));
lcs_Hz = cell(3,5);
lcs_Hz2 = cell(3,5);
lcs_adj = cell(3,5);
lcs_Hz_adj_mic = zeros(maxlcslength,5,3);
lcs_adj_mic = zeros(maxlcslength,5,3);
% 

for ch = 1:3
    for j = 1:5
       lcs_Hz{ch,j}= (lcs{ch,j}+(seg(1)-1))*(fs/nfft);
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

%% SPS vs Points
lookup_freq = [29.3;35.2;40.4;49.9;58.6;68.3;76.3;81.2]; 
ind = zeros(1,8);
for i = 1:8
    ind(i) = find(ceil(resfreqs(1:24))==ceil(lookup_freq(i)));
end

ss = find(ceil(resfreqs(1:24))==ceil(lookup_freq(8)));
indfvec_resfreqs = round(resfreqs(1:24).*(nfft/fs));
SP_modsum_pts = zeros(8,5,3);
for i = 1:3
    for j = 1:5
        for k = 1:8
            SP_modsum_pts(k,j,i) = FrealdB{1,i}(indfvec_resfreqs(ind(k))-seg(1),j);
        end
    end
end

sp_pts_fig = plot_SP_vs_pts_modsum(SP_modsum_pts,lookup_freq);

ll = FrealdB{1,1}(indfvec_resfreqs(1)-seg(1),5);
