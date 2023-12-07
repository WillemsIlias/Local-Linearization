function [f_ineq,f_eq,f_ineq_keep,f_eq_keep,paired_mom,J1,J2,J3,mData_ineq,mData_eq] = moments_w(W,KMSoptions)
%% USER-SPECIFIED FUNCTION: Moment function that depends on data only
% The moment (in)equalities are in the form
%
%       E_P[m(W_i,theta)] = E_P[f(W_i)] + g(theta)
%
% where
%
%       E_P[m_j(W_i,theta)] <= 0 for j = 1,...,J1
%       E_P[m_j(W_i,theta)] = 0  for j = J1+1,...,J1+J2
%
% This function computes the estimator for  E_P[f(W_i)], which is:
% f_j = (1/n)sum_i f_j(W_i).
%
% The user needs to specify this function.  Examples are given below.
%
% INPUT:
%   W.              Data vector W = [w_1, w2, ... , w_K], where w_k is n-by-1
%                   and n is the sample size.
%
%   KMSoptions.     This is a structure of additional inputs.  The user can
%                   add parameters to KMSoptions, say KMSoptions.params,
%                   and call KMSoptions.params in the user-specified
%                   functions.
%                   For example, in the 2x2 entry game, we include the
%                   support for the covariates and the probability that
%                   a particular support point occurs.
%
% OUTPUT:
%   f_ineq      J1-by-1 vector of moment inequalities f_j(W) for j=1,...,J1
%
%   f_eq        2*J2-by-1 vector of moment inequalities f_j(W) for j=1,...,J2
%               Note that we re-write the moment equalities as two moment
%               inequalities.  Thus, we have f_eq = [f(W) ; - f(W)], where
%               f(W) is the vector of moment equalities.
%
%   f_ineq_keep,f_eq_keep      
%               J1-by-1 and 2*J2-by-1 vector of moments to keep.  
%               If empirical moment j is too close to the boundary of the 
%               support of f_j, then we drop this moment in all future 
%               calculations.  The jth entry of f_keep is equal to 1 if we
%               keep f_j, otherwise it is equal to 0.  
%               As a concrete example, suppose that f_j(W_i) is a Bernoulli
%               random variable, so that E_P[f_j(W_i)] bounded between [0,1].
%               If f_j close to 0 or 1, we set f_keep(j,1) = 0.  
%               NB: If this boundary issue is not a problem in your
%               application, set f_keep = ones(J,1).
%
%   paired_mom  J3-by-1 vector of paired moment inequalities.
%               Let k=1,...,J3 by the kth pair of moment inequalities that
%               are highly correlated.  Say that these moment inequalities
%               are indexed by j1,j2 in {1,...,J1}.  Then
%                      paired_mom(j1,1) = paired_mom(j2,1) = k.
%
%   J1          Integer number of moment inequalities
%
%   J2.         Integer number of moment equalities
%
%   J3.         Integer number of pairs of paired moment inequalities.
%               Nb: require that 2*J3 <= J1.
%
% Below is a list of examples of moment functions.

% Threshhold for f_j being to boundary threshhold.
f_keep_threshold  = KMSoptions.f_keep_threshold ;

% Select DGP
DGP = KMSoptions.DGP;

if DGP == 8   
    %% Example: BCS Simulation
    dX = KMSoptions.dX ;                         % Number of covariates
 
    % Number of moment (in)equalities and paried moment inequalities:
    J1 = 2*dX;
    J2 = dX;
    J3 = dX;  

    % NOTE 1: In the 2-by-2 entry game, we specify the number of moment
    % inequalities and equalities to be 2*dimention of the support of the
    % background characteristics.  This is because there are two moment
    % inequalitiesand two moment equalities associated with each point in 
    % the support.

    % NOTE 2: Each moment inequality is paired.  So the number of pairs of
    % moment inequalities is J3 = J1/2.

    f_ineq = zeros(J1,1);                                                   % Preset moment inequalities and 
    f_eq = zeros(J2,1);                                                     % moment equalities
    paired_mom = zeros(J1,1);                                               % paired moment inequalities

    % Extract data:
    dataP11 = W(:,1:dX);            % Both firms enter
    dataP10 = W(:,dX+1:2*(dX));     % Firm 1 enters, firm 2 does not enter
    
    % MOMENT COMPUTATION
    % The moments f_j(Y,X) is equal to the frequency of observing either 
    % Y = (1,1) and X = x_q 
    for ii = 1:dX
        f_ineq((ii-1)*2 + 1,1) = -mean(dataP10(:,ii));                     % Moments m_{2q} in Eq 5.1 BCS
        f_ineq((ii-1)*2 + 2,1) = mean(dataP10(:,ii));                      % Moments m_{3q} in Eq 5.1 BCS
        f_eq(ii,1)             = mean(dataP11(:,ii));                      % Moments m_{1q} in Eq 5.1 BCS
    end
    
    % Concatenate the positive and negative of f_eq to transform into moment
    % inequalities
    f_eq = [f_eq;-f_eq]; 
    
    % All moment inequalities in this example are paired.  
    for jj = 1:J3
        paired_mom((jj-1)*2 + 1,1) = jj;
        paired_mom((jj-1)*2 + 2,1) = jj;
    end
    
    % NB: Upper and lower bound is 0 and 1 for all moment inequalities and
    % equalities
    UB = 1;
    LB = 0;
    ind = find( abs(f_ineq) > LB + f_keep_threshold & abs(f_ineq) < UB - f_keep_threshold);
    f_ineq_keep(ind,1) = 1;
    ind = find( abs(f_eq) > LB + f_keep_threshold & abs(f_eq) < UB - f_keep_threshold);
    f_eq_keep(ind,1) = 1;

elseif DGP == 20
    %% Example: BCS Simulation
    dX = KMSoptions.dX ;                         % Number of covariates
 
    % Number of moment (in)equalities and paried moment inequalities:
    J1 = 2*dX;
    J2 = 2*dX;
    J3 = dX;  

    % NOTE 1: In the 2-by-2 entry game, we specify the number of moment
    % inequalities and equalities to be 2*dimention of the support of the
    % background characteristics.  This is because there are two moment
    % inequalitiesand two moment equalities associated with each point in 
    % the support.

    % NOTE 2: Each moment inequality is paired.  So the number of pairs of
    % moment inequalities is J3 = J1/2.

    % Preset output variables
    f_ineq      = zeros(J1,1);                       % Preset moment inequalities and
    f_eq        = zeros(J2,1);                         % moment equalities
    paired_mom  = zeros(J1,1);                   % paired moment inequalities

    % Extract data (easier to read)
    Y1 = W(:,1);
    Y2 = W(:,2);
    X = W(:,3:end);
    n = size(Y1,1);
    % MOMENT COMPUTATION
    % The moments f_j(Y,X) is equal to the frequency of observing (Y,X):
    % (Y1=y1,Y2=y2,X1=x1,X2=x2,X3=x3).
    % Potential outcomes are:
    % (y1,y2) = (0,0)  (both firms do not enter)
    % (y1,y2) = (1,1)  (both firms enter)
    % (y1,y2) = (1,0)  (only firm 1 enters)
    % (y1,y2) = (0,1)  (Only firm 2 enters)
    for ii = 1:dX

        mData_ineq(:,(ii-1)*2 + 1) =  (Y1 == 1).*(Y2 == 0).*(X(:,ii) == 1);
        mData_ineq(:,(ii-1)*2 + 2) = -(Y1 == 1).*(Y2 == 0).*(X(:,ii) == 1);
        
        mData_eq(:,(ii-1)*2 + 1) =  (Y1==0).*(Y2 == 0).*(X(:,ii) == 1);
        mData_eq(:,(ii-1)*2 + 2) =  (Y1==1).*(Y2 == 1).*(X(:,ii) == 1);
       
    end

    % Moment inequalities for entry game (See Pg 34, eq 5.3)
    f_ineq = (sum(mData_ineq)/n)';

    % Moment equalities for entry game (See Pg 34, eq 5.1)
    f_eq   = (sum(mData_eq)/n)';
    % Concatenate the positive and negative of f_eq to transform into moment
    % inequalities
    f_eq = [f_eq ; -f_eq];
    mData_eq = [mData_eq,  -mData_eq];
    % All moment inequalities are paired in this example.  For each pair of
    % moment inequalities, set the corresponding element in paired_mom to a
    % unique indicator jj = 1,...,J3
    for jj = 1:J3
        paired_mom((jj-1)*2 + 1,1) = jj;
        paired_mom((jj-1)*2 + 2,1) = jj;
    end

    % Select moments to keep.  Note that since f_ineq and f_eq are bounded
    % between [0,1]. If they are close to the boundary, then we drop 
    % them using the 'f_ineq_keep' and 'f_eq_keep' vectors.
    f_ineq_keep = zeros(J1,1);
    f_eq_keep = zeros(2*J2,1);

    % NB: Upper and lower bound is 0 and 1 for all moment inequalities and
    % equalities
    UB = 1;
    LB = 0;
    ind = find( abs(f_ineq) > LB + f_keep_threshold & abs(f_ineq) < UB - f_keep_threshold);
    f_ineq_keep(ind,1) = 1;
    ind = find( abs(f_eq) > LB + f_keep_threshold & abs(f_eq) < UB - f_keep_threshold);
    f_eq_keep(ind,1) = 1;    
end

end
