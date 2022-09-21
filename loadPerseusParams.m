function Net = loadPerseusParams(InModuleFile, InConnectivityFile, NeuronType, NetToUpdate)
%
%  Net = loadPerseusModules(InModuleFile, InConnectivityFile, NeuronType[, NetToUpdate])
%
%  Load the neuron and connectivity parameters of the different populations composing 
%  the module as defined in the text files <InModuleFile>, <InConnectivityFile>. 
%  The format of such files is the one used by Perseus 2.x to set up the network 
%  parameters: each row is associated to a population of homogeneous neurons or 
%  homogeneous synapses connecting population. The population parameters are 
%  specified in different columns.
%
%  <InModuleFile>: The input file name of population definition in Perseus 2.x
%     format.
%  <InConnectivityFile>: The input file name of the synaptic connectivity between
%     populations in Perseus 2.x format.
%  <NeuronType>: Specifies the type of neurons composing the population defined in 
%     the input module file. Is case insensitive. The string accepted are the 
%     followings: 'LIF', leaky integrate and fire (IF) neurons (Beta is the membrane 
%     time constant); 'LIF_LUT', is the same as 'LIF' but will be used as current 
%     to rate gain function a Look Up Table (faster and less precise); 'VIF', is the 
%     VLSI IF neurons; 'LIFCA', is the adaptive LIF, the LIF with spike frequency 
%     adaptation. Constants will be associated directly to the structure
%     Net as Net.Constants.NT_LIF, Net.Constants.NT_LIFCA,
%     Net.Constants.NT_LIF_LUT, Net.Constants.NT_VIF. 
%  <NetToUpdate>: this OPTIONAL parameter is used if the population definitions
%     read from the input file have to be appended to an existent network 
%     definition.
%
%  <Net>: Is a structure with two fields (SNParam,CParam,P) grouping respectively
%    the parameters of the populations (and neurons composing them) and of the
%    connectivity (the probability to have two populations connected, the synaptic 
%    efficacies, and their variability). <Net.P> is the number of populations 
%    composing the network. <Net.SNParam> itself is a stucture of arrays of size 
%    <Net.P> related to the population parameters (<Net.SNParam.N> is the number 
%    of neurons in the populations, for instance). <Net.CParam> is a structure of
%    matrixes as <Net.CParam.c> (the connectivities), <Net.CParam.J> (the synaptic 
%    efficacies), <Net.CParam.Delta> (the relative standard deviation of the 
%    synaptic efficacies). The row indexes are related to the receiving populations,
%    while the columns are associated to the transmitting populations.
%    
%   Copyright 2013 Maurizio Mattia 
%   Version: 1.1 - Apr. 16, 2013
%

%% Recognizes the neurons type...
%
if exist('NetToUpdate','var')
   Net = NetToUpdate;
else
   Net.P = 0;
end
pop0 = Net.P + 1;


%-----
%   POPULATION PARAMETERS...
%-----

%
% Opens module input file...
%
[fid, message] = fopen(InModuleFile, 'rt');
if fid == -1
   disp(['[loadPerseusParams] Error: ' message]);
   Net = [];
   return
end

%
% Recognizes the neurons type...
%
Net.Constants.NT_VIF = 0;
Net.Constants.NT_LIF_LUT = 1;
Net.Constants.NT_LIF = 2;
Net.Constants.NT_LIFCA = 3;
switch upper(NeuronType)
   case 'VIF'
      ntype = Net.Constants.NT_VIF;
      phifun = @PhiLin;
   case 'LIF_LUT'
      ntype = Net.Constants.NT_LIF_LUT;
      phifun = @PhiExpFromLUT;
   case 'LIF'
      ntype = Net.Constants.NT_LIF;
      phifun = @PhiExp;
   case 'LIFCA'
      ntype = Net.Constants.NT_LIFCA;
      phifun = @PhiExp;
%       phifun = @PhiExpFromLUT;
   otherwise
      disp('[loadPerseusParams] Error: Unknown neuron type. Accepted values: LIF, LIF_LUT, VIF, LIFCA.');
      Net = [];
      return;
end

%
% Scans file rows...
%
rownum = 0;
while feof(fid) == 0
   rownum = rownum + 1;
   str = fgetl(fid);
   
   %
   % Removes comments...
   %

   ndx = findstr(str, '#');
   if ~isempty(ndx)
      str = str(1:(ndx(1)-1));
   end

   %
   % Reads parameters, if any...
   %
   num = sscanf(str, '%f');
   if length(num) > 9
      Net.P = Net.P + 1;
      Net.SNParam.N(Net.P,1) = num(1);
      Net.SNParam.Beta(Net.P,1) = num(6)/1000; % from ms to s.
      Net.SNParam.Theta(Net.P,1) = num(7);
      Net.SNParam.H(Net.P,1) = num(8);
      Net.SNParam.Tarp(Net.P,1) = num(9)/1000; % from ms to s.
      Net.SNParam.NExt(Net.P,1) = num(4);
      Net.SNParam.NuExt(Net.P,1) = num(5);
      Net.SNParam.JExt(Net.P,1) = num(2);
      Net.SNParam.DeltaExt(Net.P,1) = num(3);
      Net.SNParam.Type(Net.P,1) = ntype;
      Net.SNParam.Phi{Net.P,1} = phifun;
      if ntype == Net.Constants.NT_LIFCA
         Net.SNParam.AlphaC(Net.P,1) = num(10);
         Net.SNParam.TauC(Net.P,1) = num(11)/1000; % from ms to s.
         Net.SNParam.GC(Net.P,1) = num(12);
      end
%
% The following fields should be set from outside... 
%
%      Net.SNParam.Nu = ...;
%
   else
      if ~isempty(num)
         disp(['[loadPerseusParams] Warning: insufficient number of columns at row ' num2str(rownum)]);
      end
   end
end

fclose(fid);


%-----
%   CONNECTIVITY PARAMETERS...
%-----

%
% Opens module input file...
%
[fid, message] = fopen(InConnectivityFile, 'rt');
if fid == -1
   disp(message);
   Net = [];
   return
end

%
% Scans file rows...
%
rownum = 0;
while feof(fid) == 0
   rownum = rownum + 1;
   str = fgetl(fid);
   
   %
   % Removes comments...
   %
   ndx = findstr(str, '#');
   if ~isempty(ndx)
      str = str(1:(ndx(1)-1));
   end

   %
   % Reads parameters, if any...
   %
   [num,count,~,nextindex] = sscanf(str, '%f');
   if count > 0
      str = str(nextindex:end);
      ndx = findstr(str, ' ');
      num = [num' sscanf(str(ndx(1):end), '%f')']';
   end
   if length(num) == 7
      presyn = num(2) + pop0;
      postsyn = num(1) + pop0;

      Net.CParam.c(postsyn,presyn) = num(3);
      Net.CParam.J(postsyn,presyn) = num(6);
      Net.CParam.Delta(postsyn,presyn) = num(7);
      Net.CParam.DMin(postsyn,presyn) = num(4)/1000; % from ms to s.
      Net.CParam.DMax(postsyn,presyn) = num(5)/1000; % from ms to s.
   else
      if length(num) == 18
         presyn = num(2) + pop0;
         postsyn = num(1) + pop0;

         R0 = num(end);
         
         Net.CParam.c(postsyn,presyn) = num(3);
         Net.CParam.J(postsyn,presyn) = num(6)*(1-R0) + num(7)*R0;
         Net.CParam.Delta(postsyn,presyn) = sqrt(num(8)^2*(1-R0) + num(9)^2*R0);
         Net.CParam.DMin(postsyn,presyn) = num(4)/1000; % from ms to s.
         Net.CParam.DMax(postsyn,presyn) = num(5)/1000; % from ms to s.
      else
         if ~isempty(num)
            disp(['[loadPerseusParams] Warning: incorrect number of columns at row ' num2str(rownum)]);
         end
      end
   end
end

fclose(fid);

%
% Estimates the characteristic times from the delay time and 
% the excitatory or inhibitory nature of the populations...
%
TD = 0.05;
for p = pop0:Net.P
   Net.SNParam.DMin(p,1) = Net.CParam.DMin(p,p);
   Net.SNParam.TauD(p,1) = (Net.CParam.DMax(p,p) * TD - Net.CParam.DMin(p,p))/(TD-1) ...
                         - (Net.CParam.DMax(p,p) - Net.CParam.DMin(p,p)) / log(TD);
   if min(Net.CParam.J(:,p)) >= 0
      Net.SNParam.IsExc(p,1) = 1;
   else
      if max(Net.CParam.J(:,p)) <= 0
         Net.SNParam.IsExc(p,1) = 0;
      else
         Net.SNParam.IsExc(p,1) = -1;
         disp(['[loadPerseusParams] Warning: population ' num2str(p-pop0+1) ' is neither excitatory nor inhibitory.']);
      end
   end
end

%
% Selects excitatory background and foreground populations...
%
Net.ndxE = find(Net.SNParam.IsExc>0);
Net.ndxI = find(Net.SNParam.IsExc==0);
[~,Net.ndxEBg] = max(Net.SNParam.N(Net.ndxE));
if length(Net.ndxEBg) > 1
   disp('[loadPerseusParams] Warning: more than 1 background excitatory populations.');
end
Net.ndxEBg = Net.ndxE(Net.ndxEBg(1));
Net.ndxEFg = setdiff(Net.ndxE, Net.ndxEBg);
