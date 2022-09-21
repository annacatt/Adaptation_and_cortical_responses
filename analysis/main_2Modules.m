%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% TWO-MODULE NETWORK: PLOT and ANALYZE SPIKING NETWORK DATA %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script shows how makes plots and analyses as in the following
% article:
% Title: "Adaptation shapes Adaptation shapes local cortical reactivity: from
% bifurcation diagram and simulations to human physiological and pathological responses"
% Authors: Anna Cattani, Andrea Galluzzi,  Matteo Fecchio,  Andrea Pigorini,  Maurizio Mattia,  Marcello Massimini
% https://www.biorxiv.org/content/10.1101/2022.06.11.493219v1

% More in detail, this script, which paralles main_1Module.m designed for the single-module
% network, shows how to plot and analyze firing rate time-series obtained
% through simulations of two-module spiking neuron networks having specific
% excitation and adaptation levels.
% This code is structure as follows:
% 1) set excitation and adaptation levels
% 2) load corresponding firing rate simulations and preprocess data
% 3) make figures of firing rate activity
% 4) compute probability of off periods in both spontaneous and stimulus-evoked activity
% 5) time-frequency analysis
%
% Help is available by contacting the creators, Anna Cattani (acattani@bu.edu), Andrea Galluzzi, Matteo Fecchio, and Andrea Pigorini.
%
% 2022 - Anna Cattani, Andrea Galluzzi, Matteo Fecchio, Andrea Pigorini
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
num_mod = 2;                        % number of modules in the network
realization = 1;                    % trial to be plotted

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set excitation and adaptation levels
C_ext = 3283.7;                     % excitation value used for Fig. 3
g_c = 0.078;                        % adaptation values ([0.048 0.078]) used for Fig. 3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load simulated data and preprocess data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an example of loading the firing rate of spiking neuron network
% characterized by a specific excitation and adaptatation level set
% above.
path_data = '';                     % set here the correct path
load([path_data,'R_1_2M_rate_Cext_',num2str(C_ext),'_Gc_',num2str(g_c),'_NuExtStim_4_stimlen_2ms/rates.mat'])
t=rates(:,1);                       % firing rates
rates = rates(:,[2 4]);             % firing rate of the excitatory populations of module 1 and 2
for k=1:num_mod 
    Ca(:,k)=g_c*computeCaDynamicsFromNu(t/1000, rates(:,k), 1, 0.15); % fatigue
end

% ---------- stimulus artifact removal (median filtering) -----------
time_trial = [0;time_trial];
for i=1:Nstim
    art_index(i) = RelaxTime+(i-1)*CycleTime+sum(time_trial(1:i));
end

p = 19;                             % length of median filter
m = 10;                             % length of on-to-off and off-to-on transition for median filter
k = 28;                             % length of median filtered segments - must be odd
x = rates';
N = size(x,1);                      % number of channels

mdn = medfilt1(x',p)';              % note medfilt1 operates along columns
window = [ (1-cos([1:m-1]*pi/m))/2, ones(1,k), (1+cos([1:m-1]*pi/m))/2]; % define window applied to each segment

art_wind = zeros(1,size(x,2));      % art_wind is series of windows which are "on" for the median filtered segments
for i = 1:length(art_index)
    art_wind(floor(art_index(i)-(k+2*m-1)/2+1):floor(art_index(i)+(k+2*m-1)/2-1)) = window;
end

sig_wind = ones(size(art_wind))-art_wind; % sig_wind is "on" for the parts of the record that are not median filtered
clean_data = mdn.*(ones(N,1)*art_wind) + x.*(ones(N,1)*sig_wind); % Now add up windowed original and median filtered data

rates = clean_data';

%---------- spontaneous activity ------------------------------------------
rates_spont = rates(5000:10000,:);

%---------- epoching ------------------------------------------------------
xint = 1:3000;                      % signal length   
timesignal = -5000:4999;            % stimulus set to 0

Nstim = 250;                        % number of trials

rates_concatenated = [];
rates_concatenated_filt = [];
Ca_concatenated = [];

for i=1:Nstim
    t_stim = RelaxTime+(i-1)*CycleTime+1+sum(time_trial(1:i));
    interval = [t_stim+timesignal(1) t_stim+timesignal(end)];
    Ca_concatenated = [Ca_concatenated; Ca(interval(1):interval(2),:)];
    rates_concatenated = [rates_concatenated; rates(interval(1):interval(2),:)];
end

for k=1:num_mod
    rates_epoched(realization,k,:,:) = reshape(rates_concatenated(:,k),length(timesignal),Nstim);
    Ca_epoched(realization,k,:,:) = reshape(Ca_concatenated(:,k),length(timesignal),Nstim);
end
clear Ca rates_ rates_filt Data Data_filt idx_down_in idx_down_en

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make figures of firing rate activity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

rates_epoched_ar = squeeze(mean(rates_epoched,4)); rates_epoched_ar = squeeze(mean(rates_epoched_ar,1));
Ca_epoched_ar = squeeze(mean(Ca_epoched,4)); Ca_epoched_ar = squeeze(mean(Ca_epoched_ar,1));

% all the realizations together

rates_epoched_1_realiz = squeeze(rates_epoched(:,1,:,:));
rates_epoched_2_realiz = squeeze(rates_epoched(:,2,:,:));

figure;
plot([1:10000]-5000,rates_epoched_1_realiz(:,1:8))
%title('spontaneous')
xlabel('t')
ylabel('FR [Hz]')
set(gca, 'fontsize', 16)
xlim([-1000 1500])
ylim([0 200])


figure;
plot(rates_spont)
title('spontaneous')
xlabel('t')
ylabel('FR [Hz]')
set(gca, 'fontsize', 16)
legend('module 1','module 2')
xlim([0 5000])
ylim([0 140])

figure;
hold on
plot(timesignal,mean(rates_epoched_1_realiz,2),'LineWidth',2)%,'k','LineWidth',2)
xlabel('t')
ylabel('FR [Hz]')
title('stimulus-evoked')
xlim([-1000 1500])
ylim([0 140])
hold on
plot(timesignal,mean(rates_epoched_2_realiz,2),'LineWidth',2)%,'k','LineWidth',2)
legend('module 1','module 2')
set(gca, 'fontsize', 16)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time-frequency analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for k = 1:num_mod
    if k == 1
        tmp_signal = rates_epoched_1_realiz;
        tmp = mean(squeeze(rates_epoched_2_realiz(1:900,:)),1);
    elseif k==2
        tmp_signal = rates_epoched_2_realiz;
        tmp = mean(squeeze(rates_epoched_2_realiz(1:900,:)),1);
    end
    Data(k,:,:) = tmp_signal-repmat(tmp,size(tmp_signal,1),1); % subtract the mean
end

% paramteers used for timef.m (implemented in EEGLAB (Delorme and Makeig,
% 2004))
rate=1000;
cycles=3.5;
baseline=-400;
lf=5;                           % the lowest frequency is given by (1000*cycles)/baseline length
hf=45;

allerspboot=[];
allersp=[];

numsamples=size(Data,2); 

for k = 1:num_mod               % only one module
    signal=squeeze(Data(k,:,:))';

    X=[];                       % X is the time-serie made of concatenated trials 
    for i=1:size(signal,1)
        X=[X signal(i,:)];
    end

    %figure,
    [ersp,itc,~,timesTF,freqsTF,erspboot,itcboot,itcphase]=timef(X,size(signal,2),...
        [timesignal(1) timesignal(end)],rate,cycles,'alpha',0.005,'naccu',1000,'padratio',32,'plotersp','off',...
        'plotitc','off','plotphase','off','maxfreq',hf,'winsize',floor(cycles*rate/lf),'baseline',baseline,'timesout',1000);   
    
    erspboot = erspboot';
    itcboot = itcboot';

    allerspboot(k,:,:)=erspboot;
    allersp(k,:,:)=ersp;

    allitcboot(k,:,:)=itcboot;
    allitc(k,:,:)=itc;

    ERSPstat=squeeze(allersp(k,:,:));
    ERSPstat(find((ERSPstat > repmat(erspboot(:,1),[1 length(timesTF)])) & ...
        (ERSPstat < repmat(erspboot(:,2),[1 length(timesTF)])))) = 0;

    ITCstat=abs(squeeze(allitc(k,:,:)));
    ITCstat(ITCstat < repmat(itcboot(:,1),[1 length(timesTF)])) = 0;

    ERPstat_(k,:,:) = ERSPstat;
    ITCstat_(k,:,:) = ITCstat;
end

if g_c == 0.048
    bounds = [-15 15];
else
    bounds = [-20 20];
end

figure; 
imagesc(timesTF,freqsTF,squeeze(ERPstat_(1,:,:)),bounds); axis 'xy'
colorbar
xlabel('t [ms]')
ylabel('freq [Hz]')
title('module 1')
set(gca, 'fontsize', 16)
xlim([-1000 2000])
colormap('jet')

figure; 
imagesc(timesTF,freqsTF,squeeze(ERPstat_(2,:,:)),bounds); axis 'xy'
colorbar
xlabel('t [ms]')
ylabel('freq [Hz]')
title('module 2')
set(gca, 'fontsize', 16)
xlim([-1000 1500])
colormap('jet')

ERPstat_1 = squeeze(ERPstat_(1,:,:));
ERPstat_1(ERPstat_1>0) = 1;
ERPstat_1(ERPstat_1<=0) = 0;
ITC_nosup(1,:,:) = squeeze(abs(ITCstat_(1,:,:))) .* ERPstat_1;

ERPstat_2 = squeeze(ERPstat_(2,:,:));
ERPstat_2(ERPstat_2>0) = 1;
ERPstat_2(ERPstat_2<=0) = 0;
ITC_nosup(2,:,:) = squeeze(abs(ITCstat_(2,:,:))) .* ERPstat_2;

figure;
tmp_itc_mod1 = squeeze(mean(ITC_nosup(1,:,:)));
plot(timesTF,tmp_itc_mod1)
xlim([-1000 1500])
ylim([0 1])
xlabel('t [ms]')
ylabel('average ITC')
title('module 1')
set(gca, 'fontsize', 16)

h = figure,
tmp_itc_mod2 = squeeze(mean(ITC_nosup(2,:,:)));
plot(timesTF,tmp_itc_mod2)
xlim([-1000 1500])
ylim([0 1])
xlabel('t [ms]')
ylabel('average ITC')
title('module 2')
set(gca, 'fontsize', 16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
