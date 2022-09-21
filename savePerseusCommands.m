function [LifeTime time_trial]=savePerseusCommands(Nstim, NuExt,StimTime,CycleTime,RelaxTime, OutCommandFile)
%
%  savePerseusCommands(Nstim, NuExt,StimTime, OutCommandFile)
%
%  ...
%


%
% Opens command output file...
%
[fid, message] = fopen(OutCommandFile, 'wt');
if fid == -1
   disp(message);
   Net = [];
   return
end

%
% Scans commands...
%
LifeTime = 0;
fprintf(fid, 'SET_PARAM %8d %4d 4 %8g\n', RelaxTime,0, NuExt);
LifeTime = LifeTime+RelaxTime;
time_trial = randi(300,Nstim,1);
for n = 1:Nstim
   time_trial_=time_trial(n);
   fprintf(fid, 'SET_PARAM %8d %4d 4 %8g\n', StimTime,0, 0.75);
   fprintf(fid, 'SET_PARAM %8d %4d 4 %8g\n', CycleTime+time_trial_-StimTime,0, NuExt);
   LifeTime = LifeTime+CycleTime+time_trial_;
end

fclose(fid);
