clear

method      = 'KMS';    % Method - either AS or KMS
component   = 1;        % Component of theta to build confidence interval around
DGP         = 20;        % DGP = 1,2,3,4,5,6,7,8.  DGP = 5,6,7 are Games examples.
% DGP=8 is the BCS example.  If DGP = 8 and BCS = 1,
% then also calls computation for BCS confidence interval
additive    = 1;
% Set DGP = 8 and KMS = 0 to run BCS.
alpha       = 0.05;     % Significance level
n           = 1000;     % Sample size
Nmc         = 1001;        % Number of Monte Carlos
sim_lo      = 1;        % BCScode is very slow.  We split the job between many computers.
sim_hi      = 250;      % We run simulations mm = sim_lo ... sim_hi.
% Default is sim_lo = 1 and sim_hi= Nmc
KMSoptions  = KMSoptions_BCS_Simulation();

%% Extract/Save Information to KMSoptions, and set seed
KMSoptions.DGP = DGP;
KMSoptions.additive = additive;
KMSoptions.n = n;
KMSoptions.component = component;
seed = KMSoptions.seed;
B    = KMSoptions.B;                                                        % Bootstraps
stream = RandStream('mlfg6331_64','Seed',seed);
RandStream.setGlobalStream(stream)

%% Parameters
type = 'two-sided';         % Two-sided or one sided test?  Set to 'one-sided-UB' or 'one-sided-LB' or 'two-sided'
kappa = NaN;                 % Default kappa function
phi   = NaN;                % Default GMS function

%% Parameters that depend on DGP
theta_true  = [0.1;  1.1;   0.8; -0.5; -0.1; -0.5;
               0.3; -0.8;  -0.9; -0.2;    0; -0.5;
              -0.1; -0.3;   0.4];                   % True parameter vector

dim_p       = size(theta_true,1);
dX = (size(theta_true,1) -3)/2;                                            % dimension of X.
p = zeros(size(theta_true,1),1);                                       % Projection direction
p(component) = 1;
KMSoptions.S =  0;                                                    % Rho Polytope Constraints
% Param space is theta_1,theta_2 in [0,1], theta_k in
% [0,min(theta_1,theta_2)] for k = 1,2,3.
LB_theta    = [ -2*ones(1,dX*2+2) 0]';
UB_theta    = [  2*ones(1,dX*2)  0,  0,  0.85]';
A_theta     = []; 
b_theta     = [];
% We randomly select theta_0 from the parameter space
theta_0     =0.5*LB_theta + 0.5*UB_theta;   
psuppX = ones(dX,1)/dX;                                                % P(X=x), X is discrete uniform.
KMSoptions.dX = dX;
KMSoptions.psuppX = psuppX;
selp = 0.5*ones(dX,1);                                                    % prob of P(A_1=0,A_2=1) when there is multiplicity.
KMSoptions.selp = selp;                                                % NOTE: THIS IS THE OPPOSITE of DGP6,DGP5.
CVXGEN_name = strcat('csolve_DGP20_dX',num2str(dX));                                           % CVXGEN file name

% %% Compute population identification region
% stream = RandStream('mlfg6331_64','Seed',seed);
% RandStream.setGlobalStream(stream)
% stream.Substream = B + B*10^3 + 2;
% addpath ./MVNorm
% Identification_region = KMS_5_identification_region(theta_true,theta_true,LB_theta,UB_theta,A_theta,b_theta,KMSoptions);

%% Generate data 
% DATA FOR BCS SIMULATION
% Draw random variable
Z = zeros(n,2+dX,Nmc);
beta1 = theta_true(1:dX,1);  beta2 = theta_true(dX+1:2*dX,1);
delta1 = theta_true(2*dX+1); delta2 = theta_true(2*dX+2);
rho = theta_true(2*dX+3);

baseDatas = rand(n,4,Nmc); 
for mm = 1:Nmc
    rng(mm,'twister')
    % Simulate data
    epsilons = mvnrnd([0 0],[1 rho; rho 1],n); % epsilon in the model 
    % X denotes the market type indicator;
    U_X = rand(n,1);
    beta1X = zeros(n,1);
    beta2X = zeros(n,1);
    X = nan(n,dX);
    for j=1:dX
        X(:,j) = (U_X>=(j-1)/dX).*(U_X<j/dX);
        beta1X = beta1X+beta1(j)*X(:,j);
        beta2X = beta2X+beta2(j)*X(:,j);
    end

    % Initialize matrices that will contain both entry decision and market type           
    % Entry decision that indices {A_1=1,A_2=1}
    multiple = rand(n,1); % determines how multiplicity is resolved;          
    select01 = multiple < (X*selp);           

    Y = nan(n,2);
    Y(:,1) = (epsilons(:,1) >= -beta1X-delta1) + ...
        (epsilons(:,1) < -beta1X-delta1).*(epsilons(:,1) >= -beta1X).*((epsilons(:,2) < -beta2X) + ...
             (epsilons(:,2) >= -beta2X).*(epsilons(:,2) < -beta2X-delta2).*(1-select01));
    Y(:,2) = (epsilons(:,2) >= -beta2X-delta2) + ...
        (epsilons(:,2) < -beta2X-delta2).*(epsilons(:,2) >= -beta2X).*((epsilons(:,1) < -beta1X) + ...
             (epsilons(:,1) >= -beta1X).*(epsilons(:,1) < -beta1X-delta1).*select01);

    Z(:,:,mm) = [Y X];               % Save Z to pass to KMS algorithm
end

% Run BCS (if required)
name     = strcat('BCS_LL_KMS_component=',num2str(component),'_samplesize',num2str(n),'_dX',num2str(dX),'_sim_lo',num2str(sim_lo));
filename = strcat('Results/',name,'.mat');

addpath ./BCS
warning('OFF')
disp('MC started')
KMSoptions.theta_true = theta_true;
KMSoptions.alpha = alpha;
KMSoptions.seed = seed;
for mm = 218:sim_hi
    mm
    KMSoptions.W = Z(:,:,mm);
    W = Z(:,:,mm);
    KMSoptions.theta_add = [];
    t1 = tic;
    KMSoptions.BCS_EAM = 1;
    KMSoptions.LL_EAM  = 1;
    [LL_confidence_interval{mm},LL_output{mm}] = KMS_0_Main(W,theta_0,...
        p,[],LB_theta,UB_theta,A_theta,b_theta,alpha,type,method,kappa,phi,CVXGEN_name,KMSoptions);
    t_LL(mm) = toc(t1);   
    
    t2 = tic;
    KMSoptions.BCS_EAM = 0;
    KMSoptions.LL_EAM  = 0;
    [KMS_confidence_interval{mm},KMS_output{mm}] = KMS_0_Main(W,theta_0,...
            p,[],LB_theta,UB_theta,A_theta,b_theta,alpha,type,method,kappa,phi,CVXGEN_name,KMSoptions);
    t_KMS(mm) = toc(t2);
     
    t3 = tic;
    KMSoptions.BCS_EAM = 1;
    KMSoptions.LL_EAM  = 0;
    [BCS_confidence_interval{mm},BCS_output{mm}] = KMS_0_Main(W,theta_0,...
            p,[],LB_theta,UB_theta,A_theta,b_theta,alpha,type,method,kappa,phi,CVXGEN_name,KMSoptions);
    t_BCS(mm) = toc(t3);
    
    save(filename)
end
% %% Cell to vector
BCS_CI = nan(Nmc,2);
KMS_CI = nan(Nmc,2);
LL_CI = nan(Nmc,2);
for mm =  sim_lo:sim_hi
    BCS_CI(mm,:) = BCS_confidence_interval{mm};
    KMS_CI(mm,:) = KMS_confidence_interval{mm};
    LL_CI(mm,:)  = LL_confidence_interval{mm};
end

save(filename)
    


