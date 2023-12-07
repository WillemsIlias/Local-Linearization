clear

method      = 'KMS';    % Method - either AS or KMS
DGP         = 20;        % DGP = 1,2,3,4,5,6,7,8.  DGP = 5,6,7 are Games examples.
% DGP=8 is the BCS example.  If DGP = 8 and BCS = 1,
% then also calls computation for BCS confidence interval
% Set DGP = 8 and KMS = 0 to run BCS.
alpha       = 0.05;     % Significance level
B = 600;
KMSoptions  = KMSoptions_BCS_Simulation();

%% Extract/Save Information to KMSoptions, and set seed
KMSoptions.DGP = DGP;
seed = KMSoptions.seed;                                                       % Bootstraps
stream = RandStream('mlfg6331_64','Seed',seed);
RandStream.setGlobalStream(stream)

%% Parameters
type = 'two-sided';         % Two-sided or one sided test?  Set to 'one-sided-UB' or 'one-sided-LB' or 'two-sided'
kappa = NaN;                 % Default kappa function
phi   = NaN;                % Default GMS function

%% Parameters that depend on DGP
theta_true_set{1} = [  0.1;  1.1;
                       0.3; -0.8; 
                      -0.1; -0.3;   0.4]; 
theta_true_set{2} = [  0.1;  1.1;   -0.1; -0.5;
                       0.3; -0.8;      0; -0.5;
                      -0.1; -0.3;   0.4];  
theta_true_set{3} = [  0.1;  1.1;   0.8; -0.5; -0.1; -0.5;
                       0.3; -0.8;  -0.9; -0.2;    0; -0.5;
                      -0.1; -0.3;   0.4];   
for i = 1:3
    theta_true  = theta_true_set{i};                   % True parameter vector

    dim_p       = size(theta_true,1);
    dX = (size(theta_true,1) -3)/2;                                            % dimension of X.
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

    %% Compute population identification region
    stream = RandStream('mlfg6331_64','Seed',seed);
    RandStream.setGlobalStream(stream)
    stream.Substream = B + B*10^3 + 2;
    addpath ./MVNorm
    Identification_region{i} = KMS_5_identification_region(theta_true,theta_true,LB_theta,UB_theta,A_theta,b_theta,KMSoptions);
end
save('identifiedset')

