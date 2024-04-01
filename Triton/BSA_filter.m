% Generate filter used for BSA
function filter = BSA_filter(Fs,tau,tau2)
if(nargin<2)
    tau = 2E-3; % decay time constant of membrane integration
end
dt = 1/Fs;
t = 0:dt:8*tau;
if(nargin<3)
    alpha = 4; % ratio between tau & tauS
    tauS = tau/alpha;% decay time constant of synaptic current
else
    tauS = tau2;
    alpha = tau/tauS;
end
K0 = 1/(exp(-log(alpha)/(alpha-1))-exp(-log(alpha)/(alpha-1))^alpha); % PSP norm constant
K = @(t) K0*(exp(-t/tau)-exp(-t/tauS)); % Normalised PSP
filter = K(t);
filter = filter/sum(filter); % normalize filter
end
