%% Param script file
SAVE_VIDEO = 0;
RESERVOIR = 1;
LOGTXT = 0;
SYNAPSE_ORDER = 2;
RESET_VAR_PER_SAMPLE = 0; 
LOAD_EXISTING_INPUT = 0;
W_INIT_RANDOM = 0;
OPTIMAL_WEIGHTS = 1; % This parameter was updated to 1 in the revised script.
MODIFIED_WT_UPDATE = 0;
PROBABILISTIC_WT_UPDATE = 1;
REVERSE_PRUNE = 0;
RESERVOIR_2DGRID = 0; % Added based on the primary script for specifying the grid structure of the reservoir.
RNG.option = 'default'; % Updated to 'default' from 'shuffle' based on your latest script.
RNG.rng = rng; % store rng object for reference.
GET_DIST = 1; % Overall distribution of spikes in reservoir (stored in param).
DATASET = 'TI46-IF.mat';
PREPROCESSING = 1; % Added preprocessing parameter based on the new information.

%% save file names
savefilename = mfilename;

%% PARAM1 - Input Parameters
fs = 12E3; 
out_fs = 1E3; 
df = 12; 
earQ = 8; 
stepfactor = 0.25; 
differ = 1; 
agcf = 1; 
tauf = 32; 
max_AGC = 4E-4; 
BSAtau = 4E-3; 
BSAtau2 = 1E-3;
BSAfilterFac = 1;
appendS = 0;
alphaRes = [1E-10 1 0.5 0.8 1.2 1.4 1.6 1.8 2 3 5 0.25 0.65 0.4 0.6 0.7 0.9]; % Newly added to manage reservoir weights scaling.

%% PARAM1 - Neuron Parameters
tauV = 64E-3; 
tauC = 64E-3; 
tauv1e = 8E-3; 
tauv2e = 4E-3; 
tauv1i = 4E-3; 
tauv2i = 2E-3; 
RefracPeriod = 3E-3; 
Vth = 20;

%% PARAM1 - Reservoir Parameters
resSize = [5 5 5]; 
Wres = alphaRes(16) * [3 6; -2 -2]; % Utilizing alphaRes for Wres calculation.
r0 = 2; 
Kres = [0.45 0.3; 0.6 0.15]; 
f_inhibit = 0.15; 
GinMag = 8; 
InputFanout = 4;

%% PARAM1 - Classifier Parameters
Nout = 10; 
dW0 = 0.01; 
Wlim = 8; 
Cth = 10; 
DeltaC = 2;   
dC = 1; 
Iinf = 10000; 
Cfactor = 1; 
probN = 0.1; 
probP = 0.1;

%% PARAM1 - Training Parameters
dt = 1E-3;
Nfold = 5; 
Nepochs = 1; % Updated based on your latest script. Previously set to 200 in the detailed script.
samplesPerSpeakerPerClass = 10; 
numSpeakers = 5; 
samplesPerClass = numSpeakers * samplesPerSpeakerPerClass;
fracResIn = 1;
