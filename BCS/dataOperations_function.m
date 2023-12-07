%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Operations on the data for inference
%   For the sample problems, we only need to use mbar_std
%   For the MR problem, we only need mData and xi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mbar_std,mData,xi] = dataOperations_function(theta,data,kappa,DGP)
% theta denotes the parameter of interest.
% data denotes the dataset.
% kappa denotes the tuning parameter.
if DGP ==8
    % sample size;
    n = size(data,1);

    % dimension of X
    dimSX = size(data,2)/2;

    % computes auxiliary information from data;
    dataP11 = data(:,1:dimSX); % represents P(A_1=1,A_2=1)
    dataP10 = data(:,dimSX+1:2*dimSX); % represents P(A_1=1,A_2=0)

    % creates auxiliary parameters
    auxParameter = [0,theta(3:end)]; % defines beta_q for q=1,...,d_X
    auxParameter1 = theta(1) - auxParameter; % defines theta_1 - beta_q for q=1,...,d_X
    auxParameter2 = theta(2) - auxParameter; % defines theta_2 - beta_q for q=1,...,d_X

    % computes observations whose expectation should be Eq. (5.1)
    mData_1q  = dataP11 - repmat( (1-auxParameter1).*(1-auxParameter2) ,n,1); % Equalities in m_{1,q}
    mData_2q = dataP10 - repmat( auxParameter2.*(1-auxParameter1),n,1); % Inequalities in m_{2,q}
    mData_3q = repmat( auxParameter2 ,n,1) - dataP10; % Inequalities in m_{3,q}
    mData = [mData_2q, mData_3q, mData_1q]; % defines moment (in)equalities data. Note: inequalities should appear first

    % compute studentized sample averages of mData 
    epsilon = 0.000001; % introduces this parameter to avoid division by zero
    mbar_std  = mean(mData)./(std(mData) + epsilon);

    % Additional parameter needed in DR, PR, and MR inference
    xi = (1/kappa)*sqrt(n)*mbar_std; % Slackness measure in GMS
else
    Y = data(:,1:2);
    XX = data(:,3:end);
    dimSX = size(XX,2);
    psuppX = ones(dimSX,1)/dimSX; % P(X=x), X is discrete uniform.
    
    % Extract parameters (easier to read)
    beta1  = theta(1:dimSX)';
    beta2  = theta(dimSX+1:dimSX*2)';
    delta1 = theta(dimSX*2+1);
    delta2 = theta(dimSX*2+2);
    rho    = theta(dimSX*2+3);

    n = size(data,1); % sample size;

    data11 = nan(n,dimSX); data10 = nan(n,dimSX); data00 = nan(n,dimSX);
    g_eq = nan(n,2*dimSX); g_ineq = nan(n,2*dimSX);
    for i = 1:dimSX

        data11(:,i) = (XX(:,i)==1).*(Y(:,1)+Y(:,2)==2);
        data00(:,i) = (XX(:,i)==1).*(Y(:,1)+Y(:,2)==0);
        data10(:,i) = (XX(:,i)==1).*(Y(:,1)==1).*(Y(:,2)==0);            

        % Moment inequalities for entry game (1,0) (See Pg 34, eq 5.3)
        g_ineq(:,(i-1)*2 + 1) = mexBVNcdf([beta1(i),-beta2(i)-delta2],[0 0],[1,-rho;-rho,1])*psuppX(i)-data10(:,i);

        % Moment inequalities for entry game (1,0) (See Pg 34, eq 5.4)
        g_ineq(:,(i-1)*2 + 2) = data10(:,i)...
            -mexBVNcdf([(beta1(i)+delta1),-(beta2(i)+delta2)],[0 0],[1,-rho;-rho,1])*psuppX(i)...
            -mexBVNcdf([-beta1(i)-delta1,-beta2(i)],[0 0],[1,rho;rho,1])*psuppX(i) ...
            +mexBVNcdf([-beta1(i),-beta2(i)],[0 0],[1,rho;rho,1])*psuppX(i);

        % Moment equalities for entry game (0,0) (See Pg 34, eq 5.1)
        g_eq(:,(i-1)*2 + 1) = data00(:,i)-mexBVNcdf([-beta1(i),-beta2(i)],[0 0],[1,rho;rho,1])*psuppX(i);

        % Moment equalities for entry game (1,1) (See Pg 34, eq 5.2)
        g_eq(:,(i-1)*2 + 2) = data11(:,i)-mexBVNcdf([(beta1(i)+delta1),(beta2(i)+delta2)],[0 0],[1,rho;rho,1])*psuppX(i);
    end

    mData = [g_ineq, g_eq]; % defines moment (in)equalities data. Note: inequalities should appear first

    % compute studentized sample averages of mData 
    epsilon = 10^-4; % introduces this parameter to avoid division by zero

    mbar_std  = mean(mData)./max(std(mData,1),epsilon);

    % Additional parameter needed in DR, PR, and MR inference
    xi = (1/kappa)*sqrt(n)*mbar_std; % Slackness measure in GMS
end
end