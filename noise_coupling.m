function YCoupled = noised_coupling( noOfNeurons, K, noise_intensity, t, Ifunc,SigmaIn, Area, NoiseModel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noised_coupling.m
% Written by AMATH422 group edited from Joshua Goldwyn
% December,2019
% April 22, 2011
% Distributed with:

%   JHG and E-Shea-Brown, "The what and where of channel noise in the Hodgkin-Huxley equations", submitted to PLoS Computational Biology, 2011.

%%% Inputs
% noOfNeurons is the number of neurons on a row to be simulated
% K is the coupling constant
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
% A struct with following variables
%    NaFraction: - fraction of open NA channels, each row 
%    is a timestep, each column is a neuron
%    KFraction: - fraction of open K channels, each row 
%    is a timestep, each column is a neuron
%      V: Voltage traces for each neuron
%      t: timesteps
%      m: subunits
%      h: 
%      n:
%      I: Input current (both coupled+input)
%      Icoupled: just coupled current
% 

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



simSystem = zeros(nt, noOfNeurons);
YCoupled = struct();
V0 = 20*rand(1,noOfNeurons); %random initial conditions
for i =1:noOfNeurons
    m0(i) = alpham(V0(i)) / (alpham(V0(i)) + betam(V0(i))); % m
    h0(i) = alphah(V0(i)) / (alphah(V0(i)) + betah(V0(i))); % h
    n0(i) = alphan(V0(i)) / (alphan(V0(i)) + betan(V0(i))); % n
end


NaFraction = m0.^3.*h0;
KFraction = n0.^4;

% Initialize Output
YCoupled.NaFraction(1,:) = NaFraction;
YCoupled.KFraction(1,:) = KFraction;
YCoupled.V(1,:) = V0;
YCoupled.t(1) = t0;
YCoupled.m(1,:) = m0;
YCoupled.h(1,:) = h0;
YCoupled.n(1,:) = n0;

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

% noise strength
n_s = noise_intensity;

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
    % turn scalar values into arrays to account for multiple neurons
    mNoiseVec = randn(nt1,noOfNeurons)*n_s;
    % Imposing bounds on argument of sqrt functions, not directly altering dynamics of the subunits
    mNoise = @(V,m,i,k) sqrt((alpham(V)*(1-m) + betam(V)*m)/NNa) * mNoiseVec(i-1,k);
    hNoiseVec = randn(nt1,noOfNeurons)*n_s;
    hNoise = @(V,h,i,k) sqrt((alphah(V)*(1-h) + betah(V)*h)/NNa) * hNoiseVec(i-1,k);
    
    nNoiseVec = randn(nt1,noOfNeurons)*n_s;
    nNoise = @(V,n,i,k) sqrt((alphan(V)*(1-n) + betan(V)*n)/NK)  * nNoiseVec(i-1,k);
else
    mNoise = @(V,m,i,k) 0;
    hNoise = @(V,h,i,k) 0;
    nNoise = @(V,n,i,k) 0;
end


% Conductance Noise (Linaro et al Voltage Clamp)
if strfind(NoiseModel,'VClamp')
    ConductanceNoise = 1;
    NaWeiner = randn(nt1,7, noOfNeurons)*n_s;
    KWeiner = randn(nt1,4, noOfNeurons)*n_s;
    NaNoise = zeros(noOfNeurons, 7);  % Initialize
    KNoise = zeros(noOfNeurons, 4);  % Initialize
    
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
    NaHat =zeros(8, noOfNeurons);
    KHat = zeros(5,noOfNeurons);
    NaNoise = randn(8,nt1,noOfNeurons)*n_s;
    KNoise = randn(5,nt1,noOfNeurons)*n_s;
    
    
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
        for n = 1:noOfNeurons
        MCNa(1,1,n) = floor(NNa*(1-m0(n))^3*(1-h0(n)));
        MCNa(2,1,n) = floor(NNa*3*(1-m0(n))^2*m0(n)*(1-h0(n)));
        MCNa(3,1,n) = floor(NNa*3*(1-m0(n))^1*m0(n)^2*(1-h0(n)));
        MCNa(4,1,n) = floor(NNa*(1-m0(n))*m0(n)^3*(1-h0(n)));
        MCNa(1,2,n) = floor(NNa*(1-m0(n))^3*(h0(n)));
        MCNa(2,2,n) = floor(NNa*3*(1-m0(n))^2*m0(n)*(h0(n)));
        MCNa(3,2,n) = floor(NNa*3*(1-m0(n))^1*m0(n)^2*(h0(n)));
        MCNa(4,2,n) = NNa - sum(sum(MCNa(:,:,n)));
        MCK(1:4,n) = floor(NK*[(1-n0(n))^4 4*n0(n)*(1-n0(n))^3 6*n0(n)^2*(1-n0(n))^2 4*n0(n)^3*(1-n0(n))^1 ]);
        MCK(5,n) = NK-sum(sum(MCK(:,n)));
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% HERE IS THE SOLVER            %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% USING EULER FOR ODEs,         %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% EULER-MARUYAMA FOR SDEs, and  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GILLESPIE FOR MARKOV CHAIN    %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=2:nt
    
    % Input Current and coupling term
    I_coupled =  K.*Icoupled(V0);
    I = Ifunc(t(i-1)) +I_coupled;
    
    %Update subunits
    % Noise terms are non-zero for Subunit Noise model
    
    %for each neuron calculate subunit kinetics
    for k = 1:noOfNeurons
        m(k) = m0(k) + dt*(alpham(V0(k))*(1-m0(k)) - betam(V0(k))*m0(k)) + mNoise(V0(k),m0(k),i,k)*sqrt(dt);  % shifted to i-1 in function
        h(k) = h0(k) + dt*(alphah(V0(k))*(1-h0(k)) - betah(V0(k))*h0(k)) + hNoise(V0(k),h0(k),i,k)*sqrt(dt);
        n(k) = n0(k) + dt*(alphan(V0(k))*(1-n0(k)) - betan(V0(k))*n0(k)) + nNoise(V0(k),n0(k),i,k)*sqrt(dt);
        
        % Enforce boundary conditions (only necessary for subunit noise model
        m(k) = max(0,min(1,m(k)));
        h(k) = max(0,min(1,h(k)));
        n(k) = max(0,min(1,n(k)));
        
    end
    
    
    % Update Fluctuations if using conductance noise model
    if ~isempty(strfind(NoiseModel,'VClamp')) || ~isempty(strfind(NoiseModel,'FoxLuSystemSize'))
        switch NoiseModel
            case 'VClamp'  % Voltage Clamp (Linaro et al)
                for j =1:noOfNeurons
                    NaNoise(j,:) = NaNoise(j,:) + dt*(-NaNoise(j,:) ./ TauNa(V0(j))) + sqrt(dt)*(SigmaNa(V0(j)).*NaWeiner(i-1,:,j));
                     KNoise(j,:) = KNoise(j,:) + dt*(-KNoise(j,:) ./ TauK(V0(j))) + sqrt(dt)*(SigmaK(V0(j)).*KWeiner(i-1,:,j));
                end
                NaFluctuation = sum(NaNoise,2);
                KFluctuation = sum(KNoise,2);
            case 'FoxLuSystemSize'  % System Size (Fox and Lu)
              
                NaFluctuation = zeros(1,noOfNeurons);
                KFluctuation = zeros(1,noOfNeurons);
                for j =1:noOfNeurons
                    NaBar = [(1-m0(j))^3*(1-h0(j)) , 3*(1-m0(j))^2*m0(j)*(1-h0(j)) , ...
                             3*(1-m0(j))*m0(j)^2*(1-h0(j)) , m0(j)^3*(1-h0(j)) , ...
                             (1-m0(j))^3*h0(j) , 3*(1-m0(j))^2*m0(j)*h0(j) , ...
                            3*(1-m0(j))*m0(j)^2*h0(j) , m0(j)^3*h0(j)];
                  
                    
                    KBar  = [(1-n0(j))^4 , 4*n0(j)*(1-n0(j))^3 , ...
                        6*n0(j)^2*(1-n0(j))^2  , 4*n0(j)^3*(1-n0(j))  , n0(j)^4];
                    
                    
                    NaHat(:,j) = NaHat(:,j) + dt*ANa(V0(j))*NaHat(:,j) ...
                    + sqrt(dt)*SNa(V0(j),NaBar,NNa)*NaNoise(:,i-1,j);
                
                    KHat(:,j) =  KHat(:,j)  + dt*AK(V0(j))*KHat(:,j) ...
                    + sqrt(dt)*SK(V0(j),KBar)*KNoise(:,i-1,j);
                    
                    NaFluctuation(1,j) = NaHat(end,j) ;
                    KFluctuation(1,j)  = KHat(end,j) ;
                end
                
                
        end
    else
        NaFluctuation = 0;
        KFluctuation = 0;
    end
   if strfind(NoiseModel,'MarkovChain')
       for iterator=1:noOfNeurons
        [MCNa(:,:,iterator), MCK(:,iterator)]= MarkovChainFraction(V0(iterator), MCNa(:,:,iterator), MCK(:,iterator), t0,dt);
        NaFraction(iterator) = MCNa(4,2,iterator) / NNa;
        KFraction(iterator) = MCK(5,iterator) / NK;
       end
    elseif strfind(NoiseModel,'FoxLuSystemSize') 
        % Note: Impose bounds on fractions to avoid <0 or >1 in dV/dt equation, this doesn't directly alter the dynamics of the subunits or channels
        for k = 1:noOfNeurons
          
            NaFraction(k) = max(0, min(1, m0(k)^3*h0(k) + NaFluctuation(k)));  % Fluctuations are non-zero for Conductance Noise Models
            KFraction(k) = max(0, min(1, n0(k)^4 + KFluctuation(k)));
            
            
        end
    elseif strfind(NoiseModel,'VClamp')
        % Note: Impose bounds on fractions to avoid <0 or >1 in dV/dt equation, this doesn't directly alter the dynamics of the subunits or channels
        for k = 1:noOfNeurons
            NaFraction(k) = max(0, min(1, m0(k)^3*h0(k) + NaFluctuation(k)));  % Fluctuations are non-zero for Conductance Noise Models
            KFraction(k) = max(0, min(1, n0(k)^4 + KFluctuation(k)));
        end
    else
         for k = 1:noOfNeurons
          
            NaFraction(k) = max(0, min(1, m0(k)^3*h0(k) + NaFluctuation));  % Fluctuations are non-zero for Conductance Noise Models
            KFraction(k) = max(0, min(1, n0(k)^4 + KFluctuation));
            
            function YCoupled = noised_coupling( noOfNeurons, K, noise_intensity, t, Ifunc,SigmaIn, Area, NoiseModel)

        end
    end
    
    % Update Voltage
    for k = 1:noOfNeurons
        Vrhs(k) = (-gNa*(NaFraction(k)).*(V0(k) - ENa)-gK*(KFraction(k)).*(V0(k) - EK) - gL*(V0(k)-EL) + I(k))/C;
        V(k) = V0(k) + dt*Vrhs(k) + sqrt(dt)*VNoise(i-1)/C ;   % VNoise is non-zero for Current Noise Model
        
    end
    
    
    % Save Outputs
 
 
    YCoupled.NaFraction(i,:) = NaFraction;
    YCoupled.KFraction(i,:) = KFraction;
    YCoupled.V(i,:) = V;
    YCoupled.t(i) = t(i);
    YCoupled.m(i,:) = m;
    YCoupled.h(i,:) = h;
    YCoupled.n(i,:) = n;
    YCoupled.I(i,:) = I;
    YCoupled.Icoupled(i,:) = I_coupled;
    
    
    V0 =V;
    m0 = m;
    h0 = h;
    n0 = n;
    
    
end  % End loop over time for SDE solver


end % End Function Definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% END OF SOLVER                 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%coupling function for a row of neurons
function Vcoupled = Icoupled(V0)
Vcoupled =zeros(size(V0));
if size(V0,2) ==1
    Vcoupled = 0;
else
    
    for i = 1:size(V0,2)
        %if first neuron, 
        if i == 1
            Vcoupled(i) = V0(i + 1) - V0(i);
        %if last neuron,
        elseif i == size(V0,2)
            Vcoupled(i) = V0(i - 1) - V0(i);
        else
         % for all neurons in the middle there are two interactions
            Vcoupled(i) = (V0(i - 1) - V0(i)) + (V0(i + 1) - V0(i));
        end
    end
end
end
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