function savePerseusParams(Net, OutModuleFile, OutConnectivityFile, SynapseType)
%
%  savePerseusParams(Net, OutModuleFile, OutConnectivityFile[, SynapseType])
%
%  Saves the neuron and connectivity parameters of the different populations composing 
%  the module <Net> in the text files <OutModuleFile>, <OutConnectivityFile>. 
%  The format of such files is the one used by Perseus 2.x to set up the network 
%  parameters: each row is associated to a population of homogeneous neurons or 
%  homogeneous synapses connecting population. The population parameters are 
%  specified in different columns.
%
%  <OutModuleFile>: The output file name of population definition in Perseus 2.x
%     format.
%  <OutConnectivityFile>: The output file name of the synaptic connectivity between
%     populations in Perseus 2.x format.
%  <SynapseType>: An optional string parameter indicating the synapse type connecting
%     the populations. If not specified it is assumed to be 'Fixed'.
%
%  <Net>: Is a structure with three fields (SNParam,CParam,P) grouping respectively
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
%   Version: 1.2.1 - May 24, 2016
%   Version: 1.2 - May 17, 2013
%   Version: 1.1 - Apr. 16, 2013
%

if exist('SynapseType','var') == 0
   SynapseType = 'Fixed';
end

% 2016.05.24: Temporary...
if ~isfield(Net,'Constants')
   Net.Constants.NT_VIF = 0;
   Net.Constants.NT_LIF_LUT = 1;
   Net.Constants.NT_LIF = 2;
   Net.Constants.NT_LIFCA = 3;
   Net.Constants.NT_VIFCA = 4;
end


% -----
%   POPULATION PARAMETERS...
% -----

%
% Opens module output file...
%
[fid, message] = fopen(OutModuleFile, 'wt');
if fid == -1
   disp(message);
   Net = [];
   return
end

%
% Scans population parameters...
%
for n = 1:Net.P
   fprintf(fid, '%6d %8g %4g ', Net.SNParam.N(n), Net.SNParam.JExt(n), Net.SNParam.DeltaExt(n));
   fprintf(fid, '%6d %4g %8g ', Net.SNParam.NExt(n), Net.SNParam.NuExt(n), Net.SNParam.Beta(n)*1000);
   fprintf(fid, '%4g %4g %4g ', Net.SNParam.Theta(n), Net.SNParam.H(n), Net.SNParam.Tarp(n)*1000);
   if Net.SNParam(1).Type == Net.Constants(1).NT_LIFCA
      fprintf(fid, '%4g %4g %4g ', Net.SNParam.AlphaC(n), Net.SNParam.TauC(n)*1000, Net.SNParam.GC(n));
   end
   fprintf(fid, '0\n');
end

fclose(fid);


%-----
%   CONNECTIVITY PARAMETERS...
%-----

%
% Opens module input file...
%
[fid, message] = fopen(OutConnectivityFile, 'wt');
if fid == -1
   disp(message);
   Net = [];
   return
end

%
% Scans populations of synapses...
%
for postsyn = 1:Net.P
   for presyn = 1:Net.P
      if Net.CParam.c(postsyn,presyn) > 0 
         fprintf(fid, '%4d %4d %8g ', postsyn-1, presyn-1, Net.CParam.c(postsyn,presyn));
         fprintf(fid, '%4g %4g ''%s'' ', Net.CParam.DMin(postsyn,presyn)*1000, ...
            Net.CParam.DMax(postsyn,presyn)*1000, SynapseType);
         fprintf(fid, '%8g %4g \n', Net.CParam.J(postsyn,presyn), Net.CParam.Delta(postsyn,presyn));
%       else
%          fprintf(fid, '%4d %4d %8g ', postsyn-1, presyn-1, 0);
%          fprintf(fid, '%4g %4g ''%s'' ', 0, 0, SynapseType);
%          fprintf(fid, '%8g %4g \n', 0, 0);
      end
   end
end

fclose(fid);