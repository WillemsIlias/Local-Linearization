function [g_ineq,g_eq] = moments_theta(theta,J1,J2,KMSoptions)
%% USER-SPECIFIED FUNCTION: Moment function that depends on parameter only
% The moment functions are in the form
%
%       E_P[m(W_i,theta)] = E_P[f(W_i)] + g(theta)
%
% where
%
%       E_P[m_j(W_i,theta)] <= 0 for j = 1,...,J1
%       E_P[m_j(W_i,theta)] = 0  for j = J1+1,...,J1+J2
%
% This function computes the model-implied component moment (in)equalities,
% g_j(theta).
%
% The user needs to specify this function.  An example is given below.
%
% INPUT:
%   theta         d-by-1 vector of parameters
%
%   J1            Integer number of moment inequalities
%
%   J2.           Integer number of moment equalities
%
%   KMSoptions.   This is a structure of additional inputs.  The user can
%                 add parameters to KMSoptions, say KMSoptions.params,
%                 and call KMSoptions.params in the user-specified
%                 functions.
%                 For example, in the 2x2 entry game, we include the
%                 support for the covariates and the probability that
%                 a particular support point occurs.
%
% OUTPUT:
%   g_ineq    J1-by-1 vector of moment inequalities g_j(theta) for j=1,...,J1
%
%   g_eq      2*J2-by-1 vector of moment inequalities g_j(theta) for j=1,...,J2
%             Note that we re-write the moment equalities as two moment
%             inequalities.  Thus, we have
%             f_eq = [g(theta) ; - g(theta)], where g(theta) is the
%             vector of moment equalities.
%
% NOTE: User does not need to specify J1 nor J2 here, since it is already
% computed in moments_w(W,KMSoptions).
%
% Below is a list of examples of moment functions.

% Select DGP
DGP = KMSoptions.DGP;

if  DGP == 9
    %% Example: 2-by-2 Entry Game in Kline & Tamer
    % Parameter vector is
    % theta =
    % (beta^1_0, % Constant coefficient for LCC
    % beta^1_1, % Coefficient on Size for LCC
    % beta^1_2, % Coefficient on Presence for LCC
    % beta^2_0, % Constant coefficient for OA
    % beta^2_1, % Coefficient on Size for OA
    % beta^2_2, % Coefficient on Presence for OA
    % delta^1_0, % Constant competition effect for LCC
    % delta^2_0, % Constant competition effect for OA
    % rho) % Correlation in Bivariate Normal Distribution

    % PRELIMINARY
    suppX = KMSoptions.suppX;                   % Support for covariates [x1,x2,x3];
    psuppX = KMSoptions.psuppX;                 % Probability of support point occuring;
    dim_suppX = size(suppX,1);                  % Number of support points                
    g_ineq = zeros(J1,1);                     % Preset moment inequalities and
    g_eq = zeros(J2,1);                       % moment equalities

    % Extract parameters (easier to read)
    beta1  = theta(1:3);
    beta2  = theta(4:6);
    delta1 = theta(7);
    delta2 = theta(8);
    rho    = theta(9);
    
    % MOMENT COMPUTATION
    % For each point in the support, we get the theoretical frequencies of
    % (Y1=y1,Y2=y2,X1=x1,X2=x2,X3=x3).
    % Potential outcomes are:
    % (y1,y2) = (0,0)  (both firms do not enter)
    % (y1,y2) = (1,1)  (both firms enter)
    % (y1,y2) = (1,0)  (only firm 1 enters)
    % (y1,y2) = (0,1)  (Only firm 2 enters)
    for ii = 1:dim_suppX
        % Pick support point (include constant + market charactheristic)
        x1 = [1, suppX(ii,3), suppX(ii,1)];
        x2 = [1, suppX(ii,3), suppX(ii,2)];
       
        % Probability of this support point occuring
        pX= psuppX(ii);

        % Moment inequalities for entry game (See Pg 34, eq 5.3)
        g_ineq((ii-1)*2 + 1,1) = mexBVNcdf([(-x1*beta1-delta1),x2*beta2],[0;0],[1,-rho;-rho,1])*pX;
        
        % Moment inequalities for entry game (See Pg 34, eq 5.4)
        g_ineq((ii-1)*2 + 2,1) = ...
        -mexBVNcdf([(-x1*beta1-delta1),(x2*beta2+delta2)],[0;0],[1,-rho;-rho,1])*pX...
        -mexBVNcdf([-x1*beta1,(-x2*beta2-delta2)],[0;0],[1,rho;rho,1])*pX ...
        +mexBVNcdf([-x1*beta1,-x2*beta2],[0;0],[1,rho;rho,1])*pX;
    
    
        % Moment equalities for entry game (See Pg 34, eq 5.1)
        g_eq((ii-1)*2 + 1,1) = mexBVNcdf([-x1*beta1,-x2*beta2],[0;0],[1,rho;rho,1])*pX;
    
        % Moment equalities for entry game (See Pg 34, eq 5.2)
        g_eq((ii-1)*2 + 2,1) = mexBVNcdf([(x1*beta1+delta1),(x2*beta2+delta2)],[0;0],[1,rho;rho,1])*pX;
    end
    % Concatenate the positive and negative of g_eq to transform into moment
    % inequalities
    g_eq = [g_eq ;-g_eq];

    % The way that the moment functions are written, m = f + g, g enters
    % additively.  The moment functions in this example have the incorrect
    % sign.
    g_ineq  = -g_ineq;
    g_eq    = -g_eq;

elseif  DGP == 20
    %% Example: 2-by-2 Entry Game in Kline & Tamer

    % PRELIMINARY
    suppX = KMSoptions.suppX;                   % Support for covariates [x1,x2,x3];
    psuppX = KMSoptions.psuppX;                 % Probability of support point occuring;
    dim_suppX = size(suppX,1);                  % Number of support points                
    g_ineq = zeros(J1,1);                     % Preset moment inequalities and
    g_eq = zeros(J2,1);                       % moment equalities
    dX = KMSoptions.dX;
    
    % Extract parameters (easier to read)
    beta1  = theta(1:dX)';
    beta2  = theta(dX+1:dX*2)';
    delta1 = theta(dX*2+1);
    delta2 = theta(dX*2+2);
    rho    = theta(dX*2+3);
    
    % MOMENT COMPUTATION
    for ii = 1:dX
       
        % Probability of this support point occuring
        pX= psuppX(ii);

        % Moment inequalities for entry game (See Pg 34, eq 5.3)
        g_ineq((ii-1)*2 + 1,1) = -mexBVNcdf([beta1(ii),-beta2(ii)-delta2],[0 0],[1,-rho;-rho,1])*pX;
%         g_ineq((ii-1)*2 + 1,1) = mexBVNcdf([(-x1*beta1-delta1),x2*beta2],[0;0],[1,-rho;-rho,1])*pX;
        
        % Moment inequalities for entry game (See Pg 34, eq 5.4)
        g_ineq((ii-1)*2 + 2,1) = ...
            mexBVNcdf([(beta1(ii)+delta1),-(beta2(ii)+delta2)],[0 0],[1,-rho;-rho,1])*pX...
            +mexBVNcdf([-beta1(ii)-delta1,-beta2(ii)],[0 0],[1,rho;rho,1])*pX ...
            -mexBVNcdf([-beta1(ii),-beta2(ii)],[0 0],[1,rho;rho,1])*pX;
    
    
        % Moment equalities for entry game (See Pg 34, eq 5.1)
        g_eq((ii-1)*2 + 1,1) = -mexBVNcdf([-beta1(ii),-beta2(ii)],[0 0],[1,rho;rho,1])*pX;
    
        % Moment equalities for entry game (See Pg 34, eq 5.2)
        g_eq((ii-1)*2 + 2,1) = -mexBVNcdf([(beta1(ii)+delta1),(beta2(ii)+delta2)],[0 0],[1,rho;rho,1])*pX;
    end
    % Concatenate the positive and negative of g_eq to transform into moment
    % inequalities
    g_eq = [g_eq ;-g_eq];
    
end


    


end
