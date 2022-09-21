%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PLOT and ANALYZE SINGLE-MODULE SPIKING NEURON NETWORK DATA %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script shows how makes plots and analyses as in the following
% article:
% Title: "Adaptation shapes Adaptation shapes local cortical reactivity: from
% bifurcation diagram and simulations to human physiological and pathological responses"
% Authors: Anna Cattani, Andrea Galluzzi,  Matteo Fecchio,  Andrea Pigorini,  Maurizio Mattia,  Marcello Massimini
% https://www.biorxiv.org/content/10.1101/2022.06.11.493219v1
% 
% More in detail, this code present how to plot and analyze firing rate time-series
% obtained through simulations of single-module spiking neuron networks having specific
% excitation and adaptation levels.  
% This code is structured as follows:
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set excitation and adaptation levels
tmp_C_ext = 3200:3.75:3350;         % all the excitation levels considered
ggc = 20:2.5:100;                   % all the adaptation levels considered

C_ext = 3297.5;                     % excitation value used for Fig. 1
g_c = 30;                           %[30 40 45 57.5 90]; adaptation values used for Fig. 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load simulated data and preprocess data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an example of loading the firing rate of spiking neuron network
% characterized by a specific excitation and adaptatation level set
% above.
path_data = '';                     % set here the correct path
load([path_data,'1M_rate_Cext_',num2str(C_ext),'_Gc_',num2str(g_c),'_NuExtStim_4_stimlen_2ms/rates.mat'])
t = rates(:,1);                     % time
rates = rates(:,2);                 % firing rate of the excitatory population
Ca=g_c*computeCaDynamicsFromNu(t/1000, rates, 1, 0.15); % fatigue


%---------- stimulus artifact removal -------------------------------------
time_trial = [0;time_trial];

for i=1:Nstim
    stim_time = RelaxTime+(i-1)*CycleTime+sum(time_trial(1:i));
    m = rates(stim_time:stim_time+8);
    [~,max_m] = max(m);
    [~,min_max_m] = min(m(max_m:end));
    max_min_m = 1;
    min_max_m = min_max_m + max_m-1;

    min_m_id = stim_time;
    max_m_id = stim_time+min_max_m;

    val_max_min_m = m(max_min_m);
    val_min_max_m = m(min_max_m);

    x = 0:length(max_min_m:min_max_m)-1;

    q = val_max_min_m;
    m = (val_max_min_m-val_min_max_m)/(x(1)-x(end));
    y = m.*x+q;

    rates(min_m_id:max_m_id-1) = y;
end

%---------- spontaneous activity ------------------------------------------
rates_spont = rates(1:RelaxTime-1); % spontaneous firing rate of the excitatory population
rates_spont(1:4999) = [];           % remove transient in the first 5 seconds


%---------- epoching ------------------------------------------------------
xint = 1:3000;                      % signal length
timesignal = -1000:1999;            % stimulus set to 0

rates_concatenated = [];
Ca_concatenated = [];
for i=1:Nstim
    t_stim = RelaxTime+(i-1)*CycleTime+1+sum(time_trial(1:i));       % t_stim takes into account the random jitter applied to the external stimulus at each trial
    interval = [t_stim+timesignal(1) t_stim+timesignal(end)];
    Ca_concatenated = [Ca_concatenated; Ca(interval(1):interval(2),:)];
    rates_concatenated = [rates_concatenated; rates(interval(1):interval(2),:)];
end

rates_epoched = reshape(rates_concatenated,xint(end),Nstim);
Ca_epoched = reshape(Ca_concatenated,xint(end),Nstim);
rates_epoched(2500:end,:) = []; rates_epoched(1:500,:) = [];         % remove the first and last 500 ms
Ca_epoched(2500:end,:) = []; Ca_epoched(1:500,:) = [];               % remove the first and last 500 ms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Making figures of firing rate activity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
in = 92999; en = 99999; trials_to_plot = 1:3;                        % initial and final time to plot spontaneous activity, and trials to be plotted
if g_c == 45                                                         % to better reveal off periods
    in = 85999; en = 92999;
end

tmp_rates_epoched = rates_epoched(:,trials_to_plot); tmp_rates_epoched = tmp_rates_epoched(:);
tmp_Ca_epoched = Ca_epoched(:,trials_to_plot); tmp_Ca_epoched = tmp_Ca_epoched(:);

rates_forfig = [rates(in:en); tmp_rates_epoched];                    % rates
t_forfig = 1:length(rates_forfig);
Ca_forfig = g_c*computeCaDynamicsFromNu(t_forfig/1000, rates_forfig, 1, 0.15); % fatigue
rates_forfig(1:1000) = [];                                           % remove transient
Ca_forfig(1:1000) = [];                                              % remove transient

fig1 = figure();                                                     % figure of spontaneous and stimulus-evoked firing rate as a function of time
plot(rates_forfig)
xlim([0 length(rates_forfig)])
ylim([0 60])
set(fig1,'Position',[0 1000 400 300])
xlabel('t [ms]')
ylabel('FR [Hz]')
set(gca, 'fontsize', 18)
ax = gca;
ax.XRuler.Exponent = 0;
title(['g = ',num2str(g_c)])


fig2 = figure();                                                     % figure of spontaneous firing rate as a function of fatigue
plot(Ca_forfig(1:round(length(rates_forfig)/2)),rates_forfig(1:round(length(rates_forfig)/2)),'LineWidth',2)
set(fig2,'Position',[0 1000 400 300])
xlabel('fatigue')
ylabel('FR [Hz]')
xlim([0 330])
ylim([0 200])
title(['g = ',num2str(g_c)])
set(gca, 'fontsize', 18)

fig3 = figure();                                                      % figure of stimulus-evoked firing rate as a function of fatigue
for kk = 1:3
    hold on
    plot(squeeze(Ca_epoched(:,kk)),squeeze(rates_epoched(:,kk)),'LineWidth',2)

end
xlabel('fatigue')
ylabel('FR [Hz]')
title(['g = ',num2str(g_c)])
set(fig3,'Position',[0 1000 400 300])
xlim([0 330])
ylim([0 200])
set(gca, 'fontsize', 18)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute probability of off periods in both spontaneous and stimulus-evoked activity for a given excitation and adaptation level
%
rates_spont_trials = reshape(rates_spont,2000,Nstim);   % spontaneous dynamic (100 seconds) is splitted in trials lasting 2 seconds each 

fr_th = 5;                                              % fr_th and duration_th are parameters defining the off-period, i.e., an off-period is here characterized by the firing rate below fr_th Hz for more than duration_th ms
duration_th = 70;

for i = 1:Nstim
    rat_spont = rates_spont_trials(:,i);
    idx_down = diff(find(rat_spont<fr_th))';
    start_ones = strfind([0,idx_down==1],[0 1]);
    end_ones = strfind([idx_down==1,0],[1 0]);
    down_duration = end_ones - start_ones + 1;
    clear start_ones end_ones idx_down

    idx = find(down_duration>duration_th);

    if isempty(idx)==1
        down_spont(i) = 0;
    else
        down_spont(i) = 1;
    end

    rat_evo = rates_epoched(1001:end,i);
    idx_down = diff(find(rat_evo<fr_th))';

    start_ones = strfind([0,idx_down==1],[0 1]);
    end_ones = strfind([idx_down==1,0],[1 0]);
    down_duration = end_ones - start_ones + 1;
    clear start_ones end_ones

    idx = find(down_duration>70);

    if isempty(idx)==1
        down_evoked(i) = 0;
    else
        down_evoked(i) = 1;
    end

end
frac_down_spont = sum(down_spont)/Nstim;                % probability of off-periods in the spontaneous activity
frac_down_evo = sum(down_evoked)/Nstim;                 % probability of off-periods in the evoked activity


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time-frequency analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
rates_epoched = reshape(rates_concatenated,xint(end),Nstim);

fr_baseline = mean(squeeze(rates_epoched(1:900,:)),1); % average firing rate between -1000 and -100 ms
Data = rates_epoched - repmat(fr_baseline,size(rates_epoched,1),1);
Data = Data';

% parameters used for timef.m (implemented in EEGLAB (Delorme and Makeig,
% 2004))
rate=1000;
cycles=3.5;
baseline=-400;
lf=5;                               % the lowest frequency is given by (1000*cycles)/baseline length
hf=45;

numsamples=size(Data,2);
X=[];                               % X is the time-serie made of concatenated trials 
for i=1:size(Data,1)
    X=[X Data(i,:)];
end

figure,
[ersp,itc,~,timesTF,freqsTF,erspboot,itcboot,itcphase]=timef(X,size(Data,2),...
    [timesignal(1) timesignal(end)],rate,cycles,'alpha',0.005,'naccu',1000,'padratio',32,'plotersp','on',...
    'plotitc','on','plotphase','off','maxfreq',hf,'winsize',floor(cycles*rate/lf),'baseline',baseline,'timesout',1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
