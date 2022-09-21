clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters to be set
%
Nstim = 50;                                       % number of stimulations (i.e., trials)
CycleTime = 8000;                                 % duration of each trial [ms]
RelaxTime = 105000;                               % duration of the spontaneous activity at the beginning of the simulation [ms]
StimTime = 2;                                     % duration of the stimulation [ms]
NuExt = 4;
PlotFigures = 1;                                  % 0: no plot
SaveFigures = 1;                                  % 0: no save

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Files necessary for Perseo to run
%
InModuleFile='modules.ini';
InConnectivityFile='connectivity.ini';
NeuronType='LIFCA';
Net = loadPerseusParams('modules.ini', 'connectivity.ini', 'LIFCA');
%%

for Cext = 3200:3.75:3350                           % set the excitation level interval here 
    for tmp_g_c = [0.02:0.0025:0.1].*1000           % set the adaptation level interval here
        g_c = tmp_g_c/1000;
        tempdir=sprintf('Data');
        if exist(tempdir,'dir') == 0
            mkdir(tempdir);
        end
        tempdir=sprintf('Data/1M_rate_Cext_%d_Gc_%f_NuExtStim_%f_stimlen_%dms',Cext,g_c,NuExt,StimTime);
        if exist(tempdir,'dir') == 0
            mkdir(tempdir);
        end
        [status,message,messageId] = copyfile('perseo.exe', tempdir, 'f');
        [status,message,messageId] = copyfile('perseo.ini', tempdir, 'f');
        %[status,message,messageId] = copyfile('computeCaDynamicsFromNu.m', tempdir, 'f');
        %[status,message,messageId] = copyfile('plotMultimodalHistogram.m', tempdir, 'f');


        OutModuleFile=[tempdir,'/modules.ini'];
        OutConnectivityFile=[tempdir,'/connectivity.ini'];
        OutCommandFile=[tempdir,'/protocol.ini'];

        Net.SNParam.NExt(Net.ndxE)=Cext;
        Net.SNParam.GC(Net.ndxE)=g_c;

        savePerseusParams(Net, OutModuleFile, OutConnectivityFile)
        [LifeTime time_trial] = savePerseusCommands(Nstim, NuExt,StimTime,CycleTime,RelaxTime, OutCommandFile);

        cd(tempdir);
        [~,~]=system(['perseo.exe Life=',num2str(LifeTime)]);

        load('rates.dat')
        save(['rates.mat'],'rates','time_trial','Nstim','CycleTime','RelaxTime','StimTime')

        cd ..\..
    end
end
