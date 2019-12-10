function [Y varSim] = simulate_variance(sim, nSim, t, Ifunc,SigmaIn, Area, NoiseModel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edited from StochasticHH_func.m
% Written by Joshua H. Goldwyn
% April 22, 2011
% Distributed with:

%   JHG and E-Shea-Brown, "The what and where of channel noise in the Hodgkin-Huxley equations", submitted to PLoS Computational Biology, 2011.

%%% Inputs
% sim is the flag to indicate whether it is a simulation or not
% nSim number of iteration for the simulation of variance/mean
% t is vector of time values (ms)
% Ifunc is a function @(t)f(t) that returns stimulus value as a function of time (in ms)
% SigmaIn: st. dev. of current noise
% Area: Membrane Area (mu m^2)
% Noise Model (String), possible values:
%   None: 'ODE'
%   Current Noise: 'Current', must also have value for SigmaIn
%   Subunit: 'Subunit', Fox and Lu subunit model, must also have value for Area
%   Voltage Clamp Conductance: 'VClamp', Linaro et al model, must also have value for Area
%   System size Conductance: 'FoxLuSystemSize', Fox and Lu system size expansion, must also have value for Area
%   Markov Chain: 'Markov Chain', must also have value for Area

%%% Outputs
% Y(:,1) : t
% Y(:,2) : V
% Y(:,3) : fraction open Na channels
% Y(:,4) : fraction open K channels
% Y(:,5) : m
% Y(:,6) : h
% Y(:,7) : n
%varSim  : struct that contains simulations of variance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize quantities needed to run solver


% time step size
dt = t(2)-t(1);

% Number of time steps
nt = length(t);  % total
nt1 = nt-1;  % at which to solve

% Initial Values
t0 = t(1);


%simulation for variance
%no of simulation

%turn subunits into arrays if it is a simulation
if sim ==1
    nSim = nSim;
    simSystem = zeros(nt, nSim);
    varSim = struct();
    V0= zeros(1,nSim);
    m0 = alpham(V0(1)) / (alpham(V0(1)) + betam(V0(1))); % m
    h0 = alphah(V0(1)) / (alphah(V0(1)) + betah(V0(1))); % h
    n0 = alphan(V0(1)) / (alphan(V0(1)) + betan(V0(1))); % n
else
    V0 = 0;
    m0 = alpham(V0) / (alpham(V0) + betam(V0)); % m
    h0 = alphah(V0) / (alphah(V0) + betah(V0)); % h
    n0 = alphan(V0) / (alphan(V0) + betan(V0)); % n
end



NaFraction = m0^3*h0;
KFraction = n0^4;

% Initialize Output
Y = zeros(nt,7);
%Y(1,:) = [t0, V0, m0^3*h0, n0^4, m0, h0, n0];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Values

% Number of Channels
NNa = round(Area*60); % Na
NK = round(Area*18); % K

% Capacitance
C = 1; % muF /cm^2

% Na Current
gNa = 120; % mS/cm^2
ENa = 120; % mV

% K Current
gK = 36; % mS/cm^2
EK = -12; % mV

% Passive Leak
gL = 0.3; % mS / cm^2
EL = 10.6; % mV


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine Which Noise Model and Do Some Necessary Setup

% No Noise
if strfind(NoiseModel,'ODE')
    % Nothing to do
end

% Current Noise
if strfind(NoiseModel,'Current')
    VNoise = SigmaIn*randn(nt1,1);
else
    VNoise = zeros(nt1,1);
end

% Subunit Noise (FL Model)
if strfind(NoiseModel,'Subunit')
    mNoiseVec = randn(nt1,1);
    % Imposing bounds on argument of sqrt functions, not directly altering dynamics of the subunits
    mNoise = @(V,m,i) sqrt((alpham(V)*(1-m) + betam(V)*m)/NNa) * mNoiseVec(i-1);
    hNoiseVec = randn(nt1,1);
    hNoise = @(V,h,i) sqrt((alphah(V)*(1-h) + betah(V)*h)/NNa) * hNoiseVec(i-1);
    
    nNoiseVec = randn(nt1,1);
    nNoise = @(V,n,i) sqrt((alphan(V)*(1-n) + betan(V)*n)/NK)  * nNoiseVec(i-1);
else
    mNoise = @(V,m,i) 0;
    hNoise = @(V,h,i) 0;
    nNoise = @(V,n,i) 0;
end

% Conductance Noise (Linaro et al Voltage Clamp)
if strfind(NoiseModel,'VClamp')
    ConductanceNoise = 1;
    NaWeiner = randn(nt1,7);
    KWeiner = randn(nt1,4);
    NaNoise =0;  % Initialize
    KNoise =0;  % Initialize
    
    taum = @(V) 1./ (alpham(V) + betam(V));
    tauh = @(V) 1./ (alphah(V) + betah(V));
    denomNa = @(V) NNa * (alphah(V) + betah(V)).^2 .*(alpham(V) + betam(V)).^6;
    TauNa = @(V) [taum(V)./[1 2 3] ...
        tauh(V) ...
        taum(V).*tauh(V)./(taum(V) + tauh(V)) ...
        taum(V).*tauh(V)./(taum(V) + 2*tauh(V)) ...
        taum(V).*tauh(V)./(taum(V) + 3*tauh(V))];
    CovNa = @(V) [3*alphah(V).^2.*alpham(V).^5.*betam(V) ...
        3*alphah(V).^2.*alpham(V).^4.*betam(V).^2 ...
        alphah(V).^2.*alpham(V).^3.*betam(V).^3 ...
        alphah(V).*betah(V).*alpham(V).^6 ...
        3*alphah(V).*betah(V).*alpham(V).^5.*betam(V) ...
        3*alphah(V).*betah(V).*alpham(V).^4.*betam(V).^2 ...
        alphah(V).*betah(V).*alpham(V).^3.*betam(V).^3] ./denomNa(V);
    
    taun = @(V) 1./ (alphan(V) + betan(V));
    TauK = @(V) taun(V) ./ [1 2 3 4];
    CovK = @(V) [4*alphan(V).^7.*betan(V) ...
        4*alphan(V).^6.*betan(V).^2 ...
        4*alphan(V).^5.*betan(V).^3 ...
        4*alphan(V).^4.*betan(V).^4]./(NK*(alphan(V)+betan(V)).^8);
    
    SigmaNa = @(V) sqrt(2*CovNa(V) ./ TauNa(V));
    SigmaK =  @(V) sqrt(2*CovK(V) ./ TauK(V));
    
end

% Conductance Noise (FL Channel Model)
if strfind(NoiseModel,'FoxLuSystemSize')
    if sim ==1
        NaHat =zeros(8, nSim);
        KHat = zeros(5,nSim);
        NaNoise = randn(8,nt1,nSim);
        KNoise = randn(5,nt1,nSim);
    else
        NaHat = zeros(8,1);  %Initial values set to 0
        KHat = zeros(5,1); %Initial values set to 0
        NaNoise = randn(8, nt1);
        KNoise = randn(5, nt1);
        
    end
    
    % Drift Na
    ANa = @(V) ...
        [ -3*alpham(V)-alphah(V)       , betam(V)                , 0              , 0                      ,  betah(V)                 , 0                       , 0          , 0 ;
        3*alpham(V)               ,-2*alpham(V)-betam(V)-alphah(V), 2*betam(V)        , 0                      ,  0                     , betah(V)                   , 0          , 0 ;
        0                      , 2*alpham(V)             , -alpham(V)-2*betam(V)-alphah(V),  3*betam(V)        ,  0                     ,  0                      , betah(V)      , 0 ;
        0                      , 0                    , alpham(V)         , -3*betam(V)-alphah(V)        ,  0                     ,  0                      , 0          , betah(V)    ;
        alphah(V)                 , 0                    , 0              , 0                      ,  -3*alpham(V) - betah(V)     , betam(V)                   , 0          , 0 ;
        0                      , alphah(V)               , 0              , 0                      ,  3*alpham(V)              ,  -2*alpham(V)-betam(V)-betah(V)  ,   2*betam(V)  , 0 ;
        0                      , 0                    , alphah(V)         , 0                      ,  0                     ,  2*alpham(V)               ,   -alpham(V)-2*betam(V)-betah(V) , 3*betam(V)  ;
        0                      , 0                    , 0              , alphah(V)                 ,           0            ,  0                      ,  alpham(V)    , -3*betam(V)-betah(V)];
    
    % Drift K
    AK = @(V) ...
        [-4*alphan(V), betan(V)             , 0                , 0                  , 0
        4*alphan(V), -3*alphan(V)-betan(V), 2*betan(V)       , 0,                   0;
        0,        3*alphan(V),        -2*alphan(V)-2*betan(V), 3*betan(V),          0;
        0,        0,               2*alphan(V),          -alphan(V)-3*betan(V), 4*betan(V);
        0,        0,               0,                 alphan(V),          -4*betan(V)];
    
    
    % Diffusion Na : Defined in a afunction below
    
    
    % Diffusion K
    DK = @(V,X) (1/(NK)) * ...
        [   (4*alphan(V)*X(1) + betan(V)*X(2)) ,  -(4*alphan(V)*X(1) + betan(V)*X(2))                                    ,   0                                                                  , 0                                                                , 0 ;
        -(4*alphan(V)*X(1) + betan(V)*X(2)),  (4.*alphan(V)*X(1) + (3*alphan(V)+ betan(V))*X(2) + 2.*betan(V)*X(3))  ,  -(2*betan(V)*X(3) + 3*alphan(V)*X(2) )                               , 0                                                                , 0 ;
        0                                  ,  -(2*betan(V)*X(3) + 3*alphan(V)*X(2))                                   , (3*alphan(V)*X(2) + (2*alphan(V)+2*betan(V))*X(3) + 3*betan(V)*X(4)) , -(3*betan(V)*X(4) + 2*alphan(V)*X(3))                            , 0 ;
        0                                  ,  0                                                                      , -(3*betan(V)*X(4) + 2*alphan(V)*X(3))                                , (2*alphan(V)*X(3) + (alphan(V)+3*betan(V))*X(4) +4*betan(V)*X(5)), -(4*betan(V)*X(5) + alphan(V)*X(4)) ;
        0                                  ,  0                                                                      , 0                                                                    , -(4*betan(V)*X(5) + alphan(V)*X(4))                              , (alphan(V)*X(4) + 4*betan(V)*X(5)) ];
    
    
    % Take Matrix square roots numerically using SVD
    SNa = @(V,Y,NNa) mysqrtm(DNa(V,Y,NNa));
    SK = @(V,X) mysqrtm(DK(V,X));
    
end


% Markov chain
if strfind(NoiseModel,'MarkovChain')
    % Initialize channel states
    if isnumeric(Area)
        MCNa(1,1) = floor(NNa*(1-m0)^3*(1-h0));
        MCNa(2,1) = floor(NNa*3*(1-m0)^2*m0*(1-h0));
        MCNa(3,1) = floor(NNa*3*(1-m0)^1*m0^2*(1-h0));
        MCNa(4,1) = floor(NNa*(1-m0)*m0^3*(1-h0));
        MCNa(1,2) = floor(NNa*(1-m0)^3*(h0));
        MCNa(2,2) = floor(NNa*3*(1-m0)^2*m0*(h0));
        MCNa(3,2) = floor(NNa*3*(1-m0)^1*m0^2*(h0));
        MCNa(4,2) = NNa - sum(sum(MCNa));
        MCK(1:4) = floor(NK*[(1-n0)^4 4*n0*(1-n0)^3 6*n0^2*(1-n0)^2 4*n0^3*(1-n0)^1 ]);
        MCK(5) = NK-sum(sum(MCK));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% HERE IS THE SOLVER            %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% USING EULER FOR ODEs,         %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% EULER-MARUYAMA FOR SDEs, and  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GILLESPIE FOR MARKOV CHAIN    %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=2:nt
    
    % Input Current
    I = Ifunc(t(i-1));
    
    %Update subunits
    % Noise terms are non-zero for Subunit Noise model
    if  sim ==1
        m = m0 + dt*(alpham(V0(1))*(1-m0) - betam(V0(1))*m0) + mNoise(V0(1),m0,i)*sqrt(dt);  % shifted to i-1 in function
        h = h0 + dt*(alphah(V0(1))*(1-h0) - betah(V0(1))*h0) + hNoise(V0(1),h0,i)*sqrt(dt);
        n = n0 + dt*(alphan(V0(1))*(1-n0) - betan(V0(1))*n0) + nNoise(V0(1),n0,i)*sqrt(dt);
        
        % Enforce boundary conditions (only necessary for subunit noise model
        m = max(0,min(1,m));
        h = max(0,min(1,h));
        n = max(0,min(1,n));
    else
        m = m0 + dt*(alpham(V0)*(1-m0) - betam(V0)*m0) + mNoise(V0,m0,i)*sqrt(dt);  % shifted to i-1 in function
        h = h0 + dt*(alphah(V0)*(1-h0) - betah(V0)*h0) + hNoise(V0,h0,i)*sqrt(dt);
        n = n0 + dt*(alphan(V0)*(1-n0) - betan(V0)*n0) + nNoise(V0,n0,i)*sqrt(dt);
        
        % Enforce boundary conditions (only necessary for subunit noise model
        m = max(0,min(1,m));
        h = max(0,min(1,h));
        n = max(0,min(1,n));
    end
    %
    % Update Fluctuations if using conductance noise model
    if ~isempty(strfind(NoiseModel,'VClamp')) || ~isempty(strfind(NoiseModel,'FoxLuSystemSize'))
        switch NoiseModel
            case 'VClamp'  % Voltage Clamp (Linaro et al)
                NaNoise = NaNoise + dt*(-NaNoise ./ TauNa(V0)) + sqrt(dt)*(SigmaNa(V0).*NaWeiner(i-1,:));
                KNoise = KNoise + dt*(-KNoise ./ TauK(V0)) + sqrt(dt)*(SigmaK(V0).*KWeiner(i-1,:));
                NaFluctuation = sum(NaNoise);
                KFluctuation = sum(KNoise);
            case 'FoxLuSystemSize'  % System Size (Fox and Lu)
                if sim ==1
                    %deterministic part
                    NaFluctuation = zeros(1,nSim);
                    KFluctuation = zeros(1,nSim);
                    NaBar = [(1-m0)^3*(1-h0) , 3*(1-m0)^2*m0*(1-h0) , 3*(1-m0)*m0^2*(1-h0) , m0^3*(1-h0) , (1-m0)^3*h0 , 3*(1-m0)^2*m0*h0 , 3*(1-m0)*m0^2*h0 , m0^3*h0];
                    KBar  = [(1-n0)^4 , 4*n0*(1-n0)^3 , 6*n0^2*(1-n0)^2  , 4*n0^3*(1-n0)  , n0^4];
                    
                    for j =1:nSim
                        NaHat(:,j) = NaHat(:,j) + dt*ANa(V0(j))*NaHat(:,j) + sqrt(dt)*SNa(V0(j),NaBar,NNa)*NaNoise(:,i-1,j);
                        KHat(:,j) =  KHat(:,j)  + dt*AK(V0(j))*KHat(:,j)  + sqrt(dt)*SK(V0(j),KBar)*KNoise(:,i-1,j);
                        
                        NaFluctuation(j) = NaHat(end) ;
                        KFluctuation(j)  = KHat(end) ;
                    end
                    
                else
                    NaBar = [(1-m0)^3*(1-h0) , 3*(1-m0)^2*m0*(1-h0) , 3*(1-m0)*m0^2*(1-h0) , m0^3*(1-h0) , (1-m0)^3*h0 , 3*(1-m0)^2*m0*h0 , 3*(1-m0)*m0^2*h0 , m0^3*h0];
                    KBar  = [(1-n0)^4 , 4*n0*(1-n0)^3 , 6*n0^2*(1-n0)^2  , 4*n0^3*(1-n0)  , n0^4];
                    NaHat = NaHat + dt*ANa(V0)*NaHat + sqrt(dt)*SNa(V0,NaBar,NNa)*NaNoise(:,i-1);
                    KHat =  KHat  + dt*AK(V0) *KHat  + sqrt(dt)*SK(V0,KBar)*KNoise(:,i-1);
                    NaFluctuation = NaHat(end) ;
                    KFluctuation  = KHat(end) ;
                    
                end
                
                
        end
    else
        NaFluctuation = 0;
        KFluctuation = 0;
    end
    
    % Compute Fraction of open channels
    if strfind(NoiseModel,'MarkovChain')
        [MCNa, MCK]= MarkovChainFraction(V0, MCNa, MCK, t0,dt);
        NaFraction = MCNa(4,2) / NNa;
        KFraction = MCK(5) / NK;
    else
        % Note: Impose bounds on fractions to avoid <0 or >1 in dV/dt equation, this doesn't directly alter the dynamics of the subunits or channels
        NaFraction = max(0, min(1, m0^3*h0 + NaFluctuation));  % Fluctuations are non-zero for Conductance Noise Models
        KFraction = max(0, min(1, n0^4 + KFluctuation));
    end
    
    % Update Voltage
    if sim ==1
        Vrhs = (-gNa*(NaFraction).*(V0 - ENa)-gK*(KFraction).*(V0 - EK) - gL*(V0-EL) + I)/C;
        V = V0 + dt*Vrhs + sqrt(dt)*VNoise(i-1)/C ;   % VNoise is non-zero for Current Noise Model
        
    else
        Vrhs = (-gNa*(NaFraction)*(V0 - ENa)-gK*(KFraction)*(V0 - EK) - gL*(V0-EL) + I)/C;
        V = V0 + dt*Vrhs + sqrt(dt)*VNoise(i-1)/C ;   % VNoise is non-zero for Current Noise Model
    end
    % Save Outputs
    if sim ==1
        Y(i,1) = t(i);
        Y(i,2) = V(:,1);
        Y(i,3) = NaFraction(:,1);
        Y(i,4) = KFraction(:,1);
        varSim.NaFraction(i,:) = NaFraction;
        varSim.KFraction(i,:) = KFraction;
        VarSim.V(i,:) = V;
        varSim.t(i) = t(i);
    else
        Y(i,2) = V;
        Y(i,3) = NaFraction;
        Y(i,4) = KFraction;
        Y(i,1) = t(i);
        
        
        
        
    end
    V0 = V;
    
    Y(i,5) = m;
    Y(i,6) = h;
    Y(i,7) = n;
    m0 = m;
    h0 = h;
    n0 = n;
    % Keep "old values" to use in next Euler time step
    
end  % End loop over time for SDE solver


end % End Function Definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% END OF SOLVER                 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define functions used above

% subunit kinetics (Hodgkin and Huxley parameters)
function out = alpham(V)
out =  0.1 * (25-V)/ (exp((25-V)/10)-1);
end

function out = betam(V)
out = 4 * exp(-V/18);
end

function out = alphah(V)
out =  0.07 * exp(-V/20);
end

function out = betah(V)
out = 1/ (exp((30-V)/10)+1);
end

function out = alphan(V)
out = 0.01 * (10-V) / (exp((10-V)/10)-1);
end

function out = betan(V)
out = 0.125 * exp(-V/80);
end

% Computing matrix square roots with SVD
function S = mysqrtm(D)
[u,s,v] = svd(D);
S = u*sqrt(s)*v';
end

% Diffusion matrix for Na
function   D = DNa(V,Y,N)
D = zeros(8,8);
y00 = Y(1);
y10 = Y(2);
y20 = Y(3);
y30 = Y(4);
y01 = Y(5);
y11 = Y(6);
y21 = Y(7);
y31 = Y(8);

D(1,1) = (  (3*alpham(V) + alphah(V))*y00 + betam(V)*y10 + betah(V)*y01 )  ;
D(1,2) = (-3*alpham(V)*y00 - betam(V)*y10);
D(1,3) = 0;
D(1,4) = 0;
D(1,5) = - (alphah(V)*y00 + betah(V)*y01);
D(1,6) = 0;
D(1,7) = 0;
D(1,8) = 0;

D(2,1) = D(1,2);
D(2,2) = ((betam(V)+2*alpham(V))*y10 + 2*betam(V)*y20 + 3*alpham(V)*y00 + alphah(V)*y10 + betah(V)*y11) ;
D(2,3) = -(2*alpham(V)*y10 + 2*betam(V)*y20);
D(2,4) = 0;
D(2,5) = 0;
D(2,6) = -(alphah(V)*y10 + betah(V)*y11);
D(2,7) = 0;
D(2,8) = 0;

D(3,1) = D(1,3);
D(3,2) = D(2,3);
D(3,3) = ((2*betam(V)+alpham(V))*y20 + 3*betam(V)*y30 + 2*alpham(V)*y10 + alphah(V)*y20 + betah(V)*y21) ;
D(3,4) = -(alpham(V)*y20+3*betam(V)*y30);
D(3,5) = 0;
D(3,6) = 0;
D(3,7) = -(alphah(V)*y20+betah(V)*y21);
D(3,8) = 0;

D(4,1) = D(1,4);
D(4,2) = D(2,4);
D(4,3) = D(3,4);
D(4,4) = (3*betam(V)*y30 + alpham(V)*y20 + alphah(V)*y30 + betah(V)*y31);
D(4,5) = 0;
D(4,6) = 0;
D(4,7) = 0;
D(4,8) = -(alphah(V)*y30 + betah(V)*y31);

D(5,1) = D(1,5);
D(5,2) = D(2,5);
D(5,3) = D(3,5);
D(5,4) = D(4,5);
D(5,5) = (3*alpham(V)*y01 + betam(V)*y11 + betah(V)*y01 + alphah(V)*y00) ;
D(5,6) = -(3*alpham(V)*y01 + betam(V)*y11) ;
D(5,7) = 0;
D(5,8) = 0;

D(6,1) = D(1,6);
D(6,2) = D(2,6);
D(6,3) = D(3,6);
D(6,4) = D(4,6);
D(6,5) = D(5,6);
D(6,6) = ((betam(V)+2*alpham(V))*y11 + 2*betam(V)*y21 + 3*alpham(V)*y01 + betah(V)*y11 + alphah(V)*y10);
D(6,7) = -(2*alpham(V)*y11+2*betam(V)*y21);
D(6,8) = 0;

D(7,1) = D(1,7);
D(7,2) = D(2,7);
D(7,3) = D(3,7);
D(7,4) = D(4,7);
D(7,5) = D(5,7);
D(7,6) = D(6,7);
D(7,7) = ((2*betam(V)+alpham(V))*y21+3*betam(V)*y31+2*alpham(V)*y11+betah(V)*y21+alphah(V)*y20);
D(7,8) = -(alpham(V)*y21+3*betam(V)*y31);

D(8,1) = D(1,8);
D(8,2) = D(2,8);
D(8,3) = D(3,8);
D(8,4) = D(4,8);
D(8,5) = D(5,8);
D(8,6) = D(6,8);
D(8,7) = D(7,8);
D(8,8) = (3*betam(V)*y31 + alpham(V)*y21 + betah(V)*y31 + alphah(V)*y30);

D = D/(N);
end


% Markov chain
function [NaStateOut, KStateOut]= MarkovChainFraction(V, NaStateIn, KStateIn, t,dt)

tswitch = t;
Nastate = NaStateIn;
Kstate = KStateIn;
% Update Channel States
while (tswitch < (t+dt))
    
    % Determine which state switches by partitioning total rate into its 28 components
    rate(1) = 3.*alpham(V) * Nastate(1,1);
    rate(2) = rate(1) + 2.*alpham(V) * Nastate(2,1);
    rate(3) = rate(2) + 1.*alpham(V) * Nastate(3,1);
    rate(4) = rate(3) + 3.*betam(V) * Nastate(4,1);
    rate(5) = rate(4) + 2.*betam(V) * Nastate(3,1);
    rate(6) = rate(5) + 1.*betam(V) * Nastate(2,1);
    rate(7) = rate(6) + alphah(V) * Nastate(1,1);
    rate(8) = rate(7) + alphah(V) * Nastate(2,1);
    rate(9) = rate(8) + alphah(V) * Nastate(3,1);
    rate(10) = rate(9) + alphah(V) * Nastate(4,1);
    rate(11) = rate(10) + betah(V) * Nastate(1,2);
    rate(12) = rate(11) + betah(V) * Nastate(2,2);
    rate(13) = rate(12) + betah(V) * Nastate(3,2);
    rate(14) = rate(13) + betah(V) * Nastate(4,2);
    rate(15) = rate(14) + 3.*alpham(V) * Nastate(1,2);
    rate(16) = rate(15) + 2.*alpham(V) * Nastate(2,2);
    rate(17) = rate(16) + 1.*alpham(V) * Nastate(3,2);
    rate(18) = rate(17) + 3.*betam(V) * Nastate(4,2);
    rate(19) = rate(18) + 2.*betam(V) * Nastate(3,2);
    rate(20) = rate(19) + 1.*betam(V) * Nastate(2,2);
    rate(21) = rate(20) + 4.*alphan(V) * Kstate(1);
    rate(22) = rate(21) + 3.*alphan(V) * Kstate(2);
    rate(23) = rate(22) + 2.*alphan(V) * Kstate(3);
    rate(24) = rate(23) + 1.*alphan(V) * Kstate(4);
    rate(25) = rate(24) + 4.*betan(V) * Kstate(5);
    rate(26) = rate(25) + 3.*betan(V) * Kstate(4);
    rate(27) = rate(26) + 2.*betan(V) * Kstate(3);
    rate(28) = rate(27) + 1.*betan(V) * Kstate(2);
    
    % Total Transition Rate
    totalrate = rate(28);
    
    % Exponential Waiting Time Distribution
    tupdate = -log(rand()) / totalrate;
    
    % Time of Next Switching Event (Exp Rand Var)
    tswitch = tswitch + tupdate;
    
    if (tswitch < (t+dt))
        
        % Scaled Uniform RV to determine which state to switch
        r = totalrate*rand();
        
        if (r < rate(1))
            Nastate(1,1) = Nastate(1,1)-1;
            Nastate(2,1) = Nastate(2,1)+1 ;
        elseif (r < rate(2))
            Nastate(2,1) = Nastate(2,1)-1;
            Nastate(3,1) = Nastate(3,1)+1 ;
        elseif (r < rate(3))
            Nastate(3,1) = Nastate(3,1)-1;
            Nastate(4,1) = Nastate(4,1)+1 ;
        elseif (r < rate(4))
            Nastate(4,1) = Nastate(4,1)-1;
            Nastate(3,1) = Nastate(3,1)+1 ;
        elseif (r < rate(5))
            Nastate(3,1) = Nastate(3,1)-1;
            Nastate(2,1) = Nastate(2,1)+1;
        elseif (r < rate(6))
            Nastate(2,1) = Nastate(2,1)-1;
            Nastate(1,1) = Nastate(1,1)+1;
        elseif (r < rate(7))
            Nastate(1,1) = Nastate(1,1)-1;
            Nastate(1,2) = Nastate(1,2)+1;
        elseif (r < rate(8))
            Nastate(2,1) = Nastate(2,1)-1;
            Nastate(2,2) = Nastate(2,2)+1;
        elseif (r < rate(9))
            Nastate(3,1) = Nastate(3,1)-1;
            Nastate(3,2) = Nastate(3,2)+1;
        elseif (r < rate(10))
            Nastate(4,1) = Nastate(4,1)-1;
            Nastate(4,2) = Nastate(4,2)+1;
        elseif (r < rate(11))
            Nastate(1,2) = Nastate(1,2)-1;
            Nastate(1,1) = Nastate(1,1)+1;
        elseif (r < rate(12))
            Nastate(2,2) = Nastate(2,2)-1;
            Nastate(2,1) = Nastate(2,1)+1;
        elseif (r < rate(13))
            Nastate(3,2) = Nastate(3,2)-1;
            Nastate(3,1) = Nastate(3,1)+1;
        elseif (r < rate(14))
            Nastate(4,2) = Nastate(4,2)-1;
            Nastate(4,1) = Nastate(4,1)+1;
        elseif (r < rate(15))
            Nastate(1,2) = Nastate(1,2)-1;
            Nastate(2,2) = Nastate(2,2)+1;
        elseif (r < rate(16))
            Nastate(2,2) = Nastate(2,2)-1;
            Nastate(3,2) = Nastate(3,2)+1;
        elseif (r < rate(17))
            Nastate(3,2) = Nastate(3,2)-1;
            Nastate(4,2) = Nastate(4,2)+1;
        elseif (r < rate(18))
            Nastate(4,2) = Nastate(4,2)-1;
            Nastate(3,2) = Nastate(3,2)+1;
        elseif (r < rate(19))
            Nastate(3,2) = Nastate(3,2)-1;
            Nastate(2,2) = Nastate(2,2)+1;
        elseif (r < rate(20))
            Nastate(2,2) = Nastate(2,2)-1;
            Nastate(1,2) = Nastate(1,2)+1;
        elseif (r < rate(21))
            Kstate(1) = Kstate(1)-1;
            Kstate(2) = Kstate(2)+1;
        elseif (r < rate(22))
            Kstate(2) = Kstate(2)-1;
            Kstate(3) = Kstate(3)+1;
        elseif (r < rate(23))
            Kstate(3) = Kstate(3)-1;
            Kstate(4) = Kstate(4)+1;
        elseif (r < rate(24))
            Kstate(4) = Kstate(4)-1;
            Kstate(5) = Kstate(5)+1;
        elseif (r < rate(25))
            Kstate(5) = Kstate(5)-1;
            Kstate(4) = Kstate(4)+1;
        elseif (r < rate(26))
            Kstate(4) = Kstate(4)-1;
            Kstate(3) = Kstate(3)+1;
        elseif (r < rate(27))
            Kstate(3) = Kstate(3)-1;
            Kstate(2) = Kstate(2)+1;
        else
            Kstate(2) = Kstate(2)-1;
            Kstate(1) = Kstate(1)+1;
        end % End if statement
        
    end % end if tswitch<dt
    
end % end while tswitch<dt

NaStateOut = Nastate;
KStateOut = Kstate;
end % end Markov chain Gillespie update function
