
%% Program to read Impulse Response

% Created: 20.03.2021
% updated: 11.04.2021

% Author: shreejay
% shreejayshrestha@gmail.com
%%
tic
clc
clear variables

% path_txtfiles = 'D:\NTNU\MCT2018_2020\MasterThesis2021\LABMeasurements\Clean_measurement_set_for_Analysis\ImpulseResponse';
path_txtfiles = 'D:\NTNU\MCT2018_2020\MasterThesis2021\LABMeasurements\Clean_measurement_set_for_Analysis\Medisinball_clean\data_ir';
txtfiles = dir(fullfile(path_txtfiles,'*.etx'));
filenames = natsort({txtfiles.name}');
fs = 48000;
Mic_index   = cell(1,numel(filenames));

%% Checking IR signals 
reload = 1;
if reload == 1
    impres = cell(1,length(filenames));
    for i = 1:length(filenames)
        fid = fopen(fullfile(path_txtfiles,filenames{i,1}));
        TxtScan = textscan(fid,'%f%f','HeaderLines',22); % skipping first 21 lines
        impres{1,i} = TxtScan{2}; %retreiving force signal data from column 2
        fclose(fid);
        Mic_index{1,i} = extractBefore(filenames{i,1},'.');
    end
end

plotir = 0;
if plotir ==1
    tvec = cell(1,length(filenames));
    for i = 1:length(impres)
        tvec{1,i} = (0:length(impres{1,i})-1)*1/fs;
    %     subplot (6,2,i)
        plot(tvec{1,i},impres{1,i});
        xlabel('time (s)')
        ylabel('Amplitude');
        title(['Impulse Response: Mic @ P\_',int2str(i)]);
    %     xlim([-0.25 5])
        ylim ([-0.2 0.2])
        pause
    end
end

% plot([impres{1,1},impres{1,2}]) % P1_2 looks good
% legend ('P1\_1','P1\_2')

% plot([impres{1,3},impres{1,4}]) % P2_2 looks good
% legend ('P2\_1','P2\_2')% P2_2 looks good

% plot([impres{1,5},impres{1,6}]) % P3_2 looks good
% legend ('P3\_1','P3\_2')% P2_2 looks good

% plot([impres{1,7},impres{1,8},impres{1,9}]) % P4_1 looks good
% legend ('P4\_1','P4\_2','P4\_3')% P2_2 looks good

% plot([impres{1,10},impres{1,11},impres{1,12}]) % P5_1 looks good
% legend ('P5\_1','P5\_2','P5\_3')% P2_2 looks good

%% Selecting one good signal for each point and saving the matrix for further processing
impres_good = cell(1,5);
impres_good{1,1}= impres{1,2};
impres_good{1,2}= impres{1,4};
impres_good{1,3}= impres{1,6};
impres_good{1,4}= impres{1,7};
impres_good{1,5}= impres{1,10};
% save data_ir\impulseresponse.mat impres_good

%% re-check selected impulse response
checkir =0;
if checkir ==1
    figure()
    tvec = cell(1,length(filenames));
    for i = 1:length(impres_good)
        tvec{1,i} = (0:length(impres_good{1,i})-1)*1/fs;
         subplot (3,2,i)
        plot(tvec{1,i},impres_good{1,i});
        xlabel('time (s)')
        ylabel('Amplitude');
        title(['Impulse Response: Mic @ P\_',int2str(i)]);
    %     xlim([-0.25 5])
        %ylim ([-0.2 0.2])
    %     pause
    end
end

% semilogy(abs(impres_good{1,1}));

%% Frequency Domain
plotirspect = 0;
L = cell(1,5);
Fimpres = cell(1,5);
Frealir = cell(1,1);
for k = 1:5
    ires = impres_good{1,k}(1:3e5);
    noise = impres_good{1,k}(1e6:1.5e6);
%     nfft = 2^nextpow2(length(ires));
%     nfft = 8388608; % nfft of microphone signals
     nfft = 2097152; % nfft of microphone and acc signals readjusted

    fs = 48000;
    fvec = fs/nfft*([0:nfft/2-1]);
    figure(1)
    Fimpres{1,k} = fft(ires,nfft);
    Frealir{1,1}(:,k) = abs(Fimpres{1,k}(1:nfft/2));
    L{1,k} = 20*log10(abs(Fimpres{1,k}(1:nfft/2)));
    Fn = fft(noise,nfft);
    Ln = 20*log10(abs(Fn(1:nfft/2)));
    if plotirspect ==1
        subplot(3,2,k)
        semilogx(fvec,L{1,k},fvec,Ln);
        grid on
        xlim([20 200])
        title(['Impulse Response of the room on mainfloor at P:',int2str(k)])
    end
end
legend ('signal','noise','Location','EastOutside')
% save data_ir\Frealir.mat Frealir

pks = cell(1,5);
lcs = cell(1,5);
seg = (round(20*nfft/fs):round(120*nfft/fs))';

for j = 1:5
    [p,l] = findpeaks(L{1,j}(seg),'MinPeakDistance',75);
    pks{1,j} = p;
    lcs{1,j}= l;
    clear p l
end
pkslcs = 1;
if pkslcs ==1
figure()
    for j = 1:5
        subplot (3,2,j)
    %     semilogx(fvec{1,ch}(ivf{1,ch}(lim{1,1})),20*log10(Favgreal{1,ch}(ivf{1,ch}(lim{1,1}),j)),fvec{1,ch}(ivf{1,ch}(lcs{1,j}+218-1)),pks{1,j},'or')
        semilogx(fvec(seg),L{1,j}(seg),fvec(lcs{1,j}+seg(1)-1),pks{1,j},'or')
    %     for k = 1:length(f0_observed)
    %         xline(f0_observed(k),'-.o');
    %     end
    %     semilogx(fvec(seg),L{1,j}(seg),fvec(lcs{1,j}),pks{1,j},'or')
        xlim([20 100])
    %     legend('Ch1 mf')
        title(['Resonance Frequency of the room for excitation Point1\_',int2str(j)])
        grid on
    %     pause
    end
end
roomsize = [4.88 5.87 4.25];
maxmodenumber = [10,10,10];
cair = 344;
[f0room_theo,f0room_modenumbers] = calcroommodefreqs(roomsize,maxmodenumber,cair);
% save data_ir\f0room_theo.mat f0room_theo
% save data_ir\f0room_modenumbers.mat f0room_modenumbers

lcslength = zeros(1,5);
    for j = 1:5
        lcslength(1,j)= length(lcs{1,j});
    end

maxlcslength = max(max(lcslength));
lcs_Hz = cell(1,5);
lcs_Hz2 = cell(1,5);
lcs_adj = cell(1,5);
lcs_Hz_adj_ir = zeros(maxlcslength,5);
lcs_adj_ir = zeros(maxlcslength,5);

    for j = 1:5
       lcs_Hz{1,j}= (lcs{1,j}+seg(1)-1)*(fs/nfft);
       lcs_Hz2{1,j} =  lcs_Hz{1,j};
       lcs_adj {1,j}= (lcs{1,j}+seg(1)-1);
       
       if length(lcs_Hz2{1,j}) <maxlcslength
            lcs_Hz2{1,j}(end+1:maxlcslength)=0;
            lcs_adj{1,j}(end+1:maxlcslength)=0;
            lcs_Hz_adj_ir(:,j)= (lcs_Hz2{1,j});
            lcs_adj_ir(:,j)= (lcs_adj{1,j});

       elseif length(lcs_Hz2{1,j})== maxlcslength
            lcs_Hz_adj_ir(:,j)= (lcs_Hz2{1,j}); 
            lcs_adj_ir(:,j)= (lcs_adj{1,j}); 

       end
    end
   
    
sorted_lcs_Hz_adj_ir = sort(lcs_Hz_adj_ir);
lcs_Hz_adj_ir_floor = floor(lcs_Hz_adj_ir );

freq_order = zeros(100,5);

for i = 1:43
    for j = 1:5
        val = floor(lcs_Hz_adj_ir(i,j));
        if val == 0
            val = 1;
        end
        freq_order(val,j) = lcs_Hz_adj_ir(i,j);
%         min_val = floor(min(lcs_Hz_adj_ir(i,:)));
%         if floor(lcs_Hz_adj_ir(i,j)) == min_val
%             freq_order(min_val,j) = lcs_Hz_adj_ir(i,j);
%         else
%             diff = floor(lcs_Hz_adj_ir(i,j)) - min_val;
%             freq_order(i+diff,j) = lcs_Hz_adj_ir(i,j);
%         end
    end
end
freq_order_median3 = zeros(length(freq_order),1);
for i = 1:length(freq_order)
    freq_order_median3(i) = median(freq_order(i,:));
    if freq_order_median3(i)==0
         freq_order_median3(i) = max(freq_order(i,:));
    end
end

freq_order_median = nonzeros(freq_order_median3); 
freq_order_median_orig = freq_order_median;
freq_order_median(1:8) = [];
freq_order_median(49:end) = [];
freq_order_median(2:4) = [];
freq_order_median(3:4) = [];
freq_order_median(4:6) = [];
freq_order_median(5:6) = [];
freq_order_median(6:8) = [];
freq_order_median(7) = [];
freq_order_median(8) = [];
freq_order_median(8) = [];
freq_order_median(9:10) = [];
freq_order_median(12:13) = [];
freq_order_median(13) = [];
freq_order_median(16:17) = [];
freq_order_median(22) = [];
freq_order_median(23) = [];

freq_order_median(11) = [];
freq_order_median(15) = [];

room_measured_Resfreq_from_ir_updated = freq_order_median;
% save Final_Output_from_allscripts_23_04_2021\room_impulse_response\room_measured_Resfreq_from_ir_updated.mat room_measured_Resfreq_from_ir_updated  
 %% locating f0 observed
% It is important to obtain amplitude of mic at each point at the same
% natural freqneucy. So define a fixed natural freq. for all points.

% lookupf0 = f0room_theo(1:25) ;
% lookupf0 = 20:100 ;
% % f0pk_Hzloc = zeros(5,25,1);
% % f0pk_indloc = zeros(5,25,1);
% % f0pk = zeros(5,25,1);
% f0pk_Hzloc = zeros(5,81,1);
% f0pk_indloc = zeros(5,81,1);
% f0pk = zeros(5,81,1);
% 
% for j = 1:1
%     for c = 1:81
%         for r = 1:5
%             [row,~] = find(lcs_Hz_adj_ir(:,r,j)>lookupf0(c) & lcs_Hz_adj_ir(:,r,j)<lookupf0(c)+1);
%             if  isempty(row)
%                [~,closest_index] = min(abs(lcs_Hz_adj_ir(:,r,j)-lookupf0(c)));
%                f0pk_Hzloc(r,c,j) = closest_index;
%             else
%                 f0pk_Hzloc(r,c,j) = row;
%             end
%             f0pk_indloc(r,c,j) = lcs_adj_ir(f0pk_Hzloc(r,c,j),r,j);
%         end
%     end
% end   

%% Defining fixed natural frequency
% f0fix_room = zeros(2,81);
% for j = 1:1
%     for ii = 1:81
%         f0fix_room(2,ii) = median(f0pk_indloc(:,ii,1));% freq index
% %         f0fix_room(1,ii) = fvec(f0fix_room(2,ii)+(seg(1)-1)); % freq in Hz.
%         f0fix_room(1,ii) = fvec(f0fix_room(2,ii)); % freq in Hz.
%     end
% end
% room_measured_Resfreq_from_ir = f0fix_room;
% room_Obs_Resfreq = f0fix_room;
% 
% save room_measured_Resfreq_from_ir.mat room_measured_Resfreq_from_ir  
% save data_ir\peakindexin_Hz_ir.mat lcs_Hz_adj_ir

% % commonHz = cell (1,1);
% f0obs = lcs_Hz_adj_ir(:,1);
% 
% for ii = 2:4
%  f0obs = intersect(f0obs,lcs_Hz_adj_ir(:,ii));
% end
% f0obs(1)=[];
% f0obs(end+1)=[71];
% f0obs(end+1)=[51];
% f0obs(end+1)=[82];
% 
% f0floor = [29 38 49 56 70 79];
% f0obs_adj = repmat(f0obs,[1,5]);
% 


% plotcomp = 0;
% if plotcomp ==1
%     figure()
%     for j = 1:5
%            subplot (3,2,j)
%     %     semilogx(fvec{1,ch}(ivf{1,ch}(lim{1,1})),20*log10(Favgreal{1,ch}(ivf{1,ch}(lim{1,1}),j)),fvec{1,ch}(ivf{1,ch}(lcs{1,j}+218-1)),pks{1,j},'or')
%         h(1:5)= semilogx(fvec(seg),L{1,j}(seg));
%     %     hold on
%     %     semilogx(fvec(lcs{1,j}+seg(1)-1),pks{1,j},'or');
%     %     hold off
%     %     hh.Color = '#ffa500';
%     %     for k = 1:length(f0floor)
%     %         xline(f0floor(k),'--k');
%     %     end
%         c=6;
%         for k = 1:26
%            h(c)= xline(f0room_theo(k),'-.o','Color','#ff7b00');
%            c=c+1;
%         end    
% 
%     %      for k = 1:length(f0obs)
%     %         xline(f0obs(k),'-.o','Color','#005aff');
%     %     end
%     %     semilogx(fvec(seg),L{1,j}(seg),fvec(lcs{1,j}),pks{1,j},'or')
%         xlim([20 100])
%     %     legend('Ch1 mf')
% 
%         title(['Resonance Frequency of the room for excitation Point: P',int2str(j)])
%         grid on
% %          pause
%     end
%         legend(h([1,6]),{'signal','f_0 analytical'},'Location','EastOutside');
%         xlabel('Frequency(Hz)')
%         ylabel('dB re ?')
% end
% % legend(h([1,2,3,4,5,6]),{'P1','P2','P3','P4','P5','f_0 analytical'});
% 
% % % NOTE commonHz excetions
% % % ch1 69-70: P1:P3 = 70 & P4:P5 = 69
% % f0obs_adj{1,1}(5,4:5) = 69;
% % % ch1 79-80: P3 = 76 & rest = 79
% % f0obs_adj{1,1}(6,3) = 76;
% % % ch2 69-70: P2 = 73 & rest = 70
% % f0obs_adj{1,2}(5,2) = 73;
% % % ch2 79-80: P1 & P5 = 79, P2 = 80 P3&P4 = 78
% % f0obs_adj{1,2}(6,2) = 80;
% % f0obs_adj{1,2}(6,3:4) = 78;
% % avg_amp = cell(1,2);
% 
% toc
