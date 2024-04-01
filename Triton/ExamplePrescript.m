

%% This file is an example prescript which contains all the parameters
% do not modify any parameters in other files
% (unless core algorithm/other encoded parameters needs fix or update)
% Author : Ajinkya Gorad
%% Setup parameters
SAVE_VIDEO =0; % IF you wish to make a movie output for visual
RESERVOIR = 1;  % RESERVOIR = 1 (in case of LSM), RESERVOIR = 0 (normal linear classifier - NOT LSM)
LOGTXT = 0;     % LOG the output text  (set to 0, unnecessary)
SYNAPSE_ORDER= 2; % Order of dynamics synapse exhibits (=0 : Static synapse (impulse response), =1: Exponentially decaying synapse, =2 : Second order synapse)
RESET_VAR_PER_SAMPLE =0; % Variables are Reset after every sample (Do not change - look at the core code)
LOAD_EXISTING_INPUT = 0; % Set to 0, always compute response from the input (otherwise, loads previous file- filename in core code)
W_INIT_RANDOM = 1;    % INitialize Linear Classifier weights as random, set to 0 initializes it from all 0.
OPTIMAL_WEIGHTS =0;   % First compute approximate optimal weight from matrix calculations
MODIFIED_WT_UPDATE = 0; % Different rules for weight update (includ in core-code if you have your own, otherwise Do not change)
PROBABILISTIC_WT_UPDATE = 1; % Use probabilistic weight update. Suggested for better performance.
PREPROCESSING=1;            % Do preprocess the input, Load existing input overrides this preprocessing (Check in core code)
REVERSE_PRUNE = 2;% Prune based on different strategies 1- Reverse Prune, 0 :forward prune, 2 Random Prune % change FracResIn for control (default  FracResIn= 1 : no pruning)
RESERVOIR_GRID=0; % if yes, then reservoir size must be 2D % 2 : 3D grid reservoir % old term : RESERVOIR_2DGRID
RNG.option = 'default'; rng(RNG.option); % default/shuffle
RNG.rng = rng; % store rng object for reference and reproducibility
GET_DIST = 1; % overall distribution of spikes in reservoir(stored in parameter structure - check core code)
%% Load Dataset
DATASET= 'TI46-IF.mat'; % Choose DATASET = 'TI46-IF.mat'; for 5 Female speakers dataset from TI46 (Standard)
% Choose DATASET= 'PoissonJittered_jitter16_lambda40.mat'; for Poisson Dataset
% Choose DATASET= 'digitdata.mat'; for Preeti's dataset (more speakers, very few samples)
%% save file names
savefilename = mfilename;
%% PARAM1Parameters based on paper
%Zhang, Y., Li, P., Jin, Y., & Choe, Y. (2015). A digital liquid state machine with biologically inspired learning and its application to speech recognition. IEEE transactions on neural networks and learning systems, 26(11), 2635-2649.
%with some additional parameters
%% PARAM1  (Primary paramters)
%Input Parameters
fs = 12E3; out_fs = 1E3; df =12; earQ = 8; stepfactor =0.25; differ=1; agcf=1; tauf=32; max_AGC = 4E-4; BSAtau =4E-3; BSAtau2 = 1E-3;BSAfilterFac = 1;
appendS = 50;
% Neuron Parameters
tauV = 64E-3; tauC = 64E-3; tauv1e = 8E-3; tauv2e = 4E-3; tauv1i = 4E-3; tauv2i = 2E-3; RefracPeriod = 3E-3; Vth = 20;
% Reservoir Parameters %3D reservoir
resSize = [5 5 5]; Wres = [3 6;-2 -2]; r0 = 2; Kres = [0.45 0.3;0.6 0.15]; f_inhibit = 0.2; GinMag = 8; InputFanout = 4;
% Classifier Parameters
Nout = 10; dW0 = 0.01; Wlim = 8; Cth = 10; DeltaC= 2   ; dC =1; Iinf = 10000; Cfactor = 1; probN = 0.1; probP = 0.1;
dt = 1E-3;
% Training Parameters
Nfold = 5; Nepochs = 10; samplesPerSpeakerPerClass=10; numSpeakers=5; samplesPerClass = numSpeakers*samplesPerSpeakerPerClass;
 fracResIn=1;
% END PARAM1s


%% Sweep Parameters (Look at coreCode first)
Wds=[];
Nsim=10;
%alphaG=[1E-10 linspace(0.1,5,Nsim)];
alphaG=[0.5 0.65 0.8 0.9 1 1.2 1.4 1.6 1.8 2 3 5];
%alphaG=[0.1:0.2:1 1.2:0.2:5 6:0.5:15];
%alphaG = 1E-10;
Nsim = numel(alphaG);
ALPHAG=repmat(reshape(alphaG,1,1,Nsim),2,2);
Wds=repmat([3 6;-2 -2],1,1,Nsim).*ALPHAG;

