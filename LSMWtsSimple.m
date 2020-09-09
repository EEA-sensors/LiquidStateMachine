%% Liquid State Machine 
% Weights  defining the edge of chaos
clear all;
%% Reset random number generator
rng(0);
%% Default configuration settings (could be modified in GUI)
STDP = 0;       % Spike Timing Dependent Plasticity
SAVE_VIDEO=0;   % Save the simulation frame by frame for each param setting
NOISE = 0;      % noise in membrane potential of neurons
TWO_SPIKES = 0; % Only give two input spikes
POISSON =0;     % Poisson spike inputs
%% Configuration parameters
eta = 0;        % noise level (NOISE only)
lambdaIn = 10;  % Poisson input spike rate (POISSON only)
DeltaT = 0;     % Separation between two spike input (TWO_SPIKES only)
%% Reservoir parameters
W0 = [3 6;-2 -2];   % Initial weights in the reservoir (E->E, E->I; I->E, I->I) Excitatory(E), Inhibitory(I)
K0 = [0.45 0.3;0.6 0.15];   % Probability coefficients for reservoir connections
Wee  = W0(1,1); Wei= W0(1,2); Wie= W0(2,1); Wii = W0(2,2); % assign individual weights
Wm = 50;              % Maximum limit for reservoir weight synaptic scaling (SLIDER)
alpha_G = 1;            % default synaptic scaling parameter
resSize= [5 5 5];    % 3D reservoir size
r0 = 2;              % effective connectivity distance parameter for gaussian profile
f_inhibit=0.2;       % fraction of inhibitory neurons in reservoir
d_syn = 1E-3;        % Synaptic Delay in reservoir (FIXED)
%%


%% Create 3D reservoir network
% createNetworkM creates an multidimensional graph on N-D grid with
% dimensions(M) = lenght(resSize) and neuron at each ND-point (n1,...,nM)
% X : source node of network graph
% Xn : destination node of network graph
% T : synaptic delay for each neuron (FIXED)
% W : weight of each synapse
% R : coordinates of each neuron
% E : value for choosing excitatory(1)/inhibitory(-1) neuron

[X,Xn,T,W,R,E] = createNetworkM(resSize,W0,r0,K0,f_inhibit,d_syn); % gaussian probablistically connected network

%% Get secondary parameters
Nres = length(E);   % number of reservoir neurons
V = zeros(Nres,1);  % Membrane potential of reservoir Neurons
v = zeros(Nres,3);  % Synaptic traces
I = zeros(Nres,3);  % Neuron input current
spiked = logical(zeros(Nres,1)); % to store spiked neurons at an instant
spikedLog =logical(sparse([])); % to store spiked neurons in all simulated time <K
G = sparse(X,Xn,W); % Sparse matrix is better to store loosely connected networks
G0=G;    % save default configuration for network weights and connectivity

% STDP definitions
dGp = -1; dGn = 1; % weight/conductance (G) update rate
A = sparse(zeros(size(G)));
A(E>0,E>0) = 1;A(E>0,E<0) = 1; % Define coefficients for STDP
A(E<0,E>0) = 1;A(E<0,E<0) = 1; %  -;;-
A = A.*logical(G);  % get the graph connectivity matrix A

%% Input mapping of spikes
Gin=ones(Nres,1)*8;     % Defines input weights for each neuron
Iapp = zeros(Nres,1);   % Externally applied current for each neuron
inFrac=0.9;             % Fraction of reservoir  neurons which take input
%% Values from the Literature : Zhang et al. 2015
dt = 1E-3; % Define time step of simulation
K = 2000;  % Define time of simulation
alphaV = (1-dt/20E-3); % membrane decay constant (discrete)
alpha1E = (1-dt/40E-3);alpha2E =(1-dt/6E-3); % synaptic trace decay constants
alpha3stdp = (1-dt/10E-3); % used in STDP
alpha1I = (1-dt/30E-3); alpha2I = (1-dt/4E-3);
% create and assign matrices
alpha1 = zeros(Nres,1); alpha2 = alpha1; alpha3 = alpha1;
Kt0e=0;Kt0i=0;Kt0=zeros(Nres,1);
alpha1(E>0) = alpha1E; alpha2(E>0) = alpha2E; alpha3(:) = alpha3stdp;
alpha1(E<0) = alpha1I; alpha2(E<0) = alpha2I;
Vth = 20;       % Threshold voltage (mV)
% normalization for PSP
Kt0e = (1-alpha1E)*(1-alpha2E)/(alpha1E-alpha2E)/dt; % normalize by area
Kt0i = (1-alpha1I)*(1-alpha2I)/(alpha1I-alpha2I)/dt;
Kt0(E>0) = Kt0e;Kt0(E<0) = Kt0i;

%% Create figure handles and GUI
fig_handle = figure('name','simulation'); % Main figure window
set(fig_handle, 'Position', [0, 0, 1280, 600]); % Defines the size of main window
text_subplot = axes('Position',[0 0 0.2 1],'Visible','off');infoText=text(0.05,0.65,' ','FontSize',10); % Information region
raster_subplot = axes('Position',[.2 0.1 .55 .35]);     % reservoir spikes
input_subplot = axes('Position',[.2 0.6 .55 .35]);      % input Spikes (to each neuron)
psp_subplot = axes('Position',[.8 0.6 .15 .35]);        % Synaptic model/Reservoir distribution
activity_subplot = axes('Position',[.8 0.1 .15 0.35]);  % Spike activity with time
glog_subplot = axes('Position',[0.025 0.8 0.15 0.2],'Visible','off'); % Weight update (STDP only)
set(gcf,'Color','w');
%% Choose weight sample from a random density of weights
for iter = 1:10000
W0_sampled =  Wm.*rand(4,1);
WLog(:,iter) =W0_sampled;
Wee = W0_sampled(1); Wei = W0_sampled(2); Wie = -W0_sampled(3); Wii = -W0_sampled(4);
G0(E>0,E>0) = Wee*(A(E>0,E>0)~=0);  W0(1,1) = Wee;
G0(E>0,E<0) = Wei*(A(E>0,E<0)~=0);  W0(1,2) = Wei;
G0(E<0,E>0) = Wie*(A(E<0,E>0)~=0);  W0(2,1) = Wie;
G0(E<0,E<0) = Wii*(A(E<0,E<0)~=0);  W0(2,2) = Wii;

%% Run simulation
% return : total number of spikes after end of impulse, 
% and time of last spike [totalSpikes,tLastSpike] 
% Excitation/inhibition potential after last spike [Ve,Vi]
VLog =[];
GsumLog = zeros(0,2,2);

IappLog=sparse([]);
V = zeros(Nres,1); % reset valuess
v = zeros(Nres,3);
I = zeros(Nres,3);
G = G0;
spikedLog(:)=0;
k_spiked = -Inf*ones(Nres,1);
inhN = find(E<0);
for k =1:K
    Iapp(:)=0;
    
   
            if(k==10 || k==70||k==140||k==210)
                Iapp(1:floor(Nres*inFrac)) = 10000*(rand(1,floor(Nres*inFrac))>0);
                if(k==400) rng(0);end
            end
    
    Iapp(inhN) = 0; IappLog(:,end+1) = Iapp; % no current to inhibitory neurons, log applied current
    V = alphaV*V+alpha_G*G'*I(:,mod(k-1,2)+1)*dt+Gin.*Iapp*dt+NOISE*eta*randn(Nres,1)*sqrt(dt);
    spikedOld = spiked; spiked = V>Vth;             V(spiked) = 0;V(V<0) = 0;
    spiked((k-k_spiked)<=2)=0;
    k_spiked(spiked)=k; % refractory period

    v(:,1:2) = [alpha1 alpha2].*v(:,1:2)+spiked;
    I(:,mod(k,3)+1) = Kt0.*v(:,1:2)*[1 -1]'; % second order synaptic model
    
    VLog(:,end+1)=V;
    spikedLog(:,k)=spiked;
end


%
% plot(0,0);
% axes(raster_subplot);
% imagesc(VLog); colorbar;hold on;
% [s,ts]=find(spikedLog.*(E>0));   plot(ts,s,'.','MarkerSize',0.01,'Color','g');
% [s,ts]=find(spikedLog.*(E<0));    plot(ts,s,'.','MarkerSize',0.01,'Color','r'); hold off;
spikesK = full(sum(spikedLog,1));
% spikesKE= full(sum(spikedLog(E>0,:),1));
% spikesKI =full(sum(spikedLog(E<0,:),1));
% spikesK1=spikesK;
% spikesK = conv(spikesK,ones(1,20)/20); % plot spikes over 2ms window (=refractory period)
% spikesKE = conv(spikesKE,ones(1,20)/20); % plot spikes over 2ms window (=refractory period)
% spikesKI = conv(spikesKI,ones(1,20)/20); % plot spikes over 2ms window (=refractory period)
% xlim([0 k]); ylim([0 Nres]); title('Output Raster and Neuron Potential');
% xlabel('time(ms)'); ylabel('Neuron ID');
% axes(input_subplot);
% spikedLog=logical(spikedLog);
% [s,ts]=find(IappLog);   plot(ts,s,'.','MarkerSize',0.1,'Color','k'); hold on;
% [s,ts]=find(spikedLog.*(E>0));   plot(ts,s,'.','MarkerSize',0.01,'Color','b');
% [s,ts]=find(spikedLog.*(E<0));    plot(ts,s,'.','MarkerSize',0.01,'Color','r'); hold off;
% hold off; xlabel('time(ms)'); ylabel('Neuron ID');
% xlim([0 k]); ylim([0 Nres]);title('Input Spikes & Output Raster');
% axes(activity_subplot);
% plot(spikesKE(1:end),'-o','Color','b','MarkerSize',2); hold on;
% plot(spikesKI(1:end),'-o','Color','r','MarkerSize',2);
% plot(spikesK(1:end),'-o','Color','g','MarkerSize',2); hold off; xlim([1 K]);
% title('Network Activity'); legend('Excitatory','Inhibitory','Total');
% xlabel('time(ms)'); ylabel('# of firing neurons');
% set(gcf,'Color','w');set(findobj(gcf,'type','axes'),'FontName','Consolas','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);

 % return function objects
    totalSpikes(iter) = sum(spikesK(211:end));
    tLastSpike(iter) = find(spikesK,1,'last');
    Ve(iter) = sum(sum(VLog(E>0,211:end)));
    Vi(iter) = sum(sum(VLog(E<0,211:end)));
    fprintf('%i: W : [%.2f,%.2f,%.2f,%.2f] (%i,%i) [%.2f,%.2f]\r\n',iter,W0_sampled,totalSpikes(end),tLastSpike(end),Ve(end),Vi(end));
end

%% total spikes vs weights
subplot(141);plot(WLog(1,:),totalSpikes,'.');ylabel('#spikes');xlabel('Wee ');
subplot(142);plot(WLog(2,:),totalSpikes,'.');xlabel('Wei ');
subplot(143);plot(WLog(3,:),totalSpikes,'.');xlabel('Wie ');
subplot(144);plot(WLog(4,:),totalSpikes,'.');xlabel('Wii ');
title('Total Spikes');
%% total time last spike vs weights
subplot(141);plot(WLog(1,:),tLastSpike,'.');ylabel('#spikes');xlabel('Wee '); ylim([0 600]);
subplot(142);plot(WLog(2,:),tLastSpike,'.');xlabel('Wei '); ylim([0 600]);
subplot(143);plot(WLog(3,:),tLastSpike,'.');xlabel('Wie '); ylim([0 600]);
subplot(144);plot(WLog(4,:),tLastSpike,'.');xlabel('Wii '); ylim([0 600]);
title('Time Last Spike');

%% total spikes vs weights
subplot(241);plot(WLog(1,:),totalSpikes,'.');ylabel('#spikes');xlabel('Wee ');
subplot(242);plot(WLog(2,:),totalSpikes,'.');xlabel('Wei ');
subplot(243);plot(WLog(3,:),totalSpikes,'.');xlabel('Wie ');
subplot(244);plot(WLog(4,:),totalSpikes,'.');xlabel('Wii ');
%title('Total Spikes');
% total time last spike vs weights
subplot(245);plot(WLog(1,:),tLastSpike,'.');ylabel('\tau');xlabel('Wee '); ylim([200 600]);
subplot(246);plot(WLog(2,:),tLastSpike,'.');xlabel('Wei '); ylim([200 600]);
subplot(247);plot(WLog(3,:),tLastSpike,'.');xlabel('Wie '); ylim([200 600]);
subplot(248);plot(WLog(4,:),tLastSpike,'.');xlabel('Wii '); ylim([200 600]);
%title('Time Last Spike');

%% Combine above to 3d figures
subplot(141);plot3(WLog(1,:),WLog(4,:),totalSpikes,'.');zlabel('#spikes');xlabel('Wee ');ylabel('Wii');
subplot(142);plot3(WLog(2,:),WLog(3,:),totalSpikes,'.');zlabel('#spikes');xlabel('Wei ');ylabel('Wie');
subplot(143);plot3(WLog(1,:),WLog(4,:),tLastSpike,'.');zlabel('\tau');xlabel('Wee ');ylabel('Wii');
subplot(144);plot3(WLog(2,:),WLog(3,:),tLastSpike,'.');zlabel('\tau');xlabel('Wei ');ylabel('Wie');
title('Total Spikes');



%% Plot total Spikes vs tLast spikes
plot(totalSpikes,tLastSpike,'.');
%%
subplot(221);scatter(totalSpikes,tLastSpike,2+0.5*WLog(1,:),WLog(1,:)); xlabel('#spikes');ylabel('\tau ');
subplot(222);scatter(totalSpikes,tLastSpike,2+0.5*WLog(2,:),WLog(2,:));
subplot(223);scatter(totalSpikes,tLastSpike,2+0.5*WLog(3,:),WLog(3,:));
subplot(224);scatter(totalSpikes,tLastSpike,2+0.5*WLog(4,:),WLog(4,:));  
%% total spikes vs weights
subplot(141);plot(WLog(1,:),totalSpikes,'.');ylabel('#spikes');xlabel('Wee ');
subplot(142);plot(WLog(2,:),totalSpikes,'.');xlabel('Wei ');
subplot(143);plot(WLog(3,:),totalSpikes,'.');xlabel('Wie ');
subplot(144);plot(WLog(4,:),totalSpikes,'.');xlabel('Wii ');

%% Plot Ve vs Vi
plot(Ve,Vi,'.');
%% Plot in a weights space
% Wee = +ve x, Wii = -ve x, Wei = +ve y, Wii = -y
colormap hsv;
nVe = Ve/2E4; nVi = Vi/2E4; nV = [Ve;Vi];
scatter(WLog(1,:),WLog(2,:),2+4*nVe,nVe);  hold on;  %1st Qdrnt
scatter(-WLog(4,:),-WLog(3,:),2+4*nVe,nVe); hold on; %3nd Qdrnt
scatter(WLog(1,:),-WLog(3,:),2+4*nVi,nVi);  hold on; %4nd Qdrnt 
scatter(-WLog(4,:),WLog(2,:),2+4*nVi,nVi); hold on;  %2nd Qdrnt

%% Plot Ve vs Vi @W
subplot(221);scatter(Ve,Vi,2+0.5*WLog(1,:),WLog(1,:));
subplot(222);scatter(Ve,Vi,2+0.5*WLog(2,:),WLog(2,:));
subplot(223);scatter(Ve,Vi,2+0.5*WLog(3,:),WLog(3,:));
subplot(224);scatter(Ve,Vi,2+0.5*WLog(4,:),WLog(4,:));  
%% Plot Ve vs Vi @W
subplot(221);plot3(Ve,Vi,WLog(1,:),'.'); xlabel('V_e');ylabel('V_i');zlabel('W_{ee}');
subplot(222);plot3(Ve,Vi,WLog(2,:),'.'); xlabel('V_e');ylabel('V_i');zlabel('W_{ei}');
subplot(223);plot3(Ve,Vi,WLog(3,:),'.'); xlabel('V_e');ylabel('V_i');zlabel('W_{ie}');
subplot(224);plot3(Ve,Vi,WLog(4,:),'.'); xlabel('V_e');ylabel('V_i');zlabel('W_{ii}');

%% Plot Ve vs Vi @W
tS = totalSpikes/max(totalSpikes); tLS = tLastSpike/max(tLastSpike);
subplot(221);scatter3(Ve,Vi,WLog(1,:),tS,tLS,'.'); xlabel('V_e');ylabel('V_i');zlabel('W_{ee}');
subplot(222);scatter3(Ve,Vi,WLog(2,:),tS,tLS,'.'); xlabel('V_e');ylabel('V_i');zlabel('W_{ei}');
subplot(223);scatter3(Ve,Vi,WLog(3,:),tS,tLS,'.'); xlabel('V_e');ylabel('V_i');zlabel('W_{ie}');
subplot(224);scatter3(Ve,Vi,WLog(4,:),tS,tLS,'.'); xlabel('V_e');ylabel('V_i');zlabel('W_{ii}');


