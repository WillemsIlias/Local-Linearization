function [c_val,minQn,theta_prof] = BCS_32_Critval(lambda,component,KMSoptions,f_ineq,f_eq,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq)

DGP = KMSoptions.DGP;
B           = KMSoptions.B;
seed        = KMSoptions.seed;
theta_true  = KMSoptions.theta_true;
W           = KMSoptions.W;
alpha       = KMSoptions.alpha;
n = size(W,1);
coordinate = component;
LL_EAM = KMSoptions.LL_EAM;
lambdaLL = sqrt(n)*log(n);
if DGP == 20
    %%   
    % Primitives of DGP    
    theta0 = reshape(theta_true,max(size(theta_true)),1);

    dimSX = (size(theta0,1) -3)/2; % dimension of X.

    % Set parameters
    R = B;       % # of bootstrap draws
    
    % Tuning parameters for inference;
    kappa = sqrt(log(n));  % GMS Thresholding parameter, preferred choice Eq. 4.4 in Andrews and Soares (2010).
    tolerance = 0.0001; % define what we mean by "close enough".
    
    
    p = 2*dimSX; % total inequalities;
    k = 4*dimSX; % total (in)equalties;
    
    options = optimset('Display','off','Algorithm','active-set'); % (this is faster).

    % collect entry data in a single matrix;
    data =  W; %data = [dataP11,dataP10];

    % BCS's starting values
    theta_H0 = lambda;
    % determines the value of the other thetas;
    theta_other0 = theta0; theta_other0(coordinate) = [];  
        
    % Compute test statistic for all tests, profile test statistic.
    % Choose starting values for minimization;
    starting_values = NaN(2,size(theta_other0,1)); % initialize;

    % first starting value: true parameter value (unfeasible in practice)
    starting_values(1,:) = theta_other0';

    % second starting value: randomly chosen (feasible in practice)
    lbi = KMSoptions.LB_theta; lbi(coordinate) = [];
    ubi = KMSoptions.UB_theta; ubi(coordinate) = [];
    starting_values(2,:) = (lbi+ubi)/2;  

    % Pick a large value as initial function value for the minimization
    min_value = 10^10;
    min_outcomes = []; % This matrix will collect results of minimization;

    % solve numerical minimization for all starting values;
    for s=1:size(starting_values,1)
        % NOTE: Qn_function defines the S-function. It assumes the coefficient vector to be 2-dimensional. 'x' represents the subvector of 'free'
        % coefficients. 'theta_H0' represents the value of the element of the coefficient vector that is fixed (also denoted 'lambda' above) That
        % is, the one for which we are testing whether there exists a coefficient vector with 'lambda' as its 'coordinate'-th component. 
        [theta_aux,Qn_aux,bandera] =  fmincon(@(x) Qn_function(x,theta_H0,coordinate,data,p,DGP),...
            starting_values(s,:),[],[],[],[],lbi,ubi,[],options);
        % check whether minimization is successful and reduced value
        if Qn_aux  < min_value && bandera >= 1
            Qn_minimizer = theta_aux;
            min_value = Qn_aux;
        end
        % if minimization is successful, collect minimizer and its value;
        if bandera >= 1
            min_outcomes = [min_outcomes;[min_value, Qn_minimizer]]; %#ok<AGROW>
        end
    end

     % at the end of the minimizations, we should have a minimizer 'theta_prof', and a corresponding minimal value
     % of the test statistic 'minQn'.
    if ~isempty(Qn_minimizer)
        minQn = min_value;
        theta_prof = [Qn_minimizer(1:coordinate-1)';lambda;Qn_minimizer(coordinate:end)'];
    else
        theta_prof = [];
        minQn = [];
    end
    
    % Collect minimizers to estimate the indetified set (used for DR)
    if size(min_outcomes,1)>1
        min_outcomes = uniquetol2(min_outcomes,tolerance,'rows'); % set of minimizers;
    end
    In_identified_set_hat = min_outcomes(min_outcomes(:,1)<=min(min_outcomes(:,1)) + tolerance , 2:end); % estimator of id set

    if LL_EAM == 1
        CVXGEN_name = strcat('CVXGEN_dX',num2str(dimSX),'_LL');
        call_cvx = strcat('@(params,settings)',CVXGEN_name,'(params,settings)');
        call_cvx = str2func(call_cvx);
    
        ub = KMSoptions.UB_theta;
        lb = KMSoptions.LB_theta;
        epsilon = sqrt(log(log(n))/n)/2*(ub-lb);
%         options2 = optimoptions('quadprog','Display','off');
        
        settings.verbose = 0;  % disable output of solver progress.
        settings.eps = 1e-9;  % reduce the required objective tolerance, from 1e-6.
            
        for kk = 1:size(In_identified_set_hat,1)

            % Recreate the full parameter vector theta_hat based on the values of the non-restricted elements for which the
            % S-function is minimized.
            t_hat_other = In_identified_set_hat(kk,:);
            t_hat = [t_hat_other(1:coordinate-1),theta_H0,t_hat_other(coordinate:end)];          

            % Compute the (standardized) moments
            [g_ineq, g_eq] = moments_theta(t_hat',p,k-p,KMSoptions);
            std_hat = [f_stdev_ineq; f_stdev_eq(1:k-p)]; % column vector
            mbar_std = [ -(f_ineq + g_ineq) ; (f_eq(1:p) + g_eq(1:p))]'./std_hat'; % KMS have opposite sign for inequalities 

            % Compute the GMS function
            phi = zeros(1,size(std_hat,1));
            phi(1:p) = max(sqrt(n)*mbar_std(1:p)/kappa,0);   

            % Not sure, but this likely computes the matrix G_n as in the algorithm in Bei, 2024. Note again that because m(W, theta) =
            % f(W) + g(theta), std_hat does not depend on theta and can hence be taken out of the derivative.
            Jacobian = diag(std_hat)^(-1)*getgradient(t_hat);
            Jacobian(:, coordinate) = [];        
            Jbar = [Jacobian  [eye(p); zeros(k-p,p)]];
            H = Jbar'*Jbar + lambdaLL/n*[eye(size(theta0,1)-1) zeros(size(theta0,1)-1,p); zeros(p,size(theta0,1)-1+p)];

            lb_boot = sqrt(n)*( lb - t_hat' + epsilon).*(lb - t_hat' + epsilon < 0);  
            lb_boot(coordinate) = []; lb_boot = [lb_boot' -10^10*ones(1,p)];
            
            ub_boot = sqrt(n)*( ub - t_hat' - epsilon).*(ub - t_hat' - epsilon > 0);  
            ub_boot(coordinate) = []; ub_boot = [ub_boot' zeros(1,p)];
            
            Qn_boot = nan(size(In_identified_set_hat,1),B);
            
            params.lb = lb_boot';
            params.ub = ub_boot';
            
            for r=1:B           
                Zbar = [-G_ineq(:,r)',G_eq(1:k-p,r)'] + phi;                    
                f = Jbar'*Zbar';  
%                 [x, fval] = quadprog(H,f,[],[],[],[],lb_boot,ub_boot,zeros(size(theta0,1)-1+p,1),options2);
%                 delta = x(1:end-p);

                params.H = H;
                params.f = f;                
                [vars, status] = call_cvx(params,settings);
                x = vars.x;
                delta = x(1:end-p);
                fval = status.optval;
                Qn_boot(kk,r) = fval*2+ Zbar*Zbar'-lambdaLL/n*(delta'*delta); 
            end       
        end
        c_val = quantile(min(Qn_boot,[],1), 1-alpha); 
    else
        if KMSoptions.additive == 1

            % CV calculation
            parfor r = 1:R
                rng(r,'twister');
                %- Step 1: Discard resampling method
                % Pick a large value as initial function value for the minimization;
                min_value_DR = 10^10;

                minQn_DR_aux = [];
                % Minimize but restricted to points of the estimated indentified set;
                for IDsetHat_index = 1:size(In_identified_set_hat,1)
                    minQn_DR_aux = min(min_value_DR,Qn_MR_function_additive(In_identified_set_hat(IDsetHat_index,:),theta_H0,...
                        coordinate,f_ineq,f_eq,f_stdev_ineq,f_stdev_eq,G_ineq(:,r),G_eq(:,r),KMSoptions,kappa,p,k,1));    
                end

                % compute simulated DR criterion function
                minQn_DR(1,r) = minQn_DR_aux;

                %- Step 2: Penalized resampling method;
                % Starting values for PR minimization;
                starting_values_r =  uniquetol2([starting_values;Qn_minimizer],tolerance,'rows');

                % Pick a large value as initial function value for the minimization;
                min_value_PR = 10^10;
                % check whether minimization is successful and reduced value
                for s=1:size(starting_values_r,1)
                    [theta_aux,Qn_aux,bandera] =  fmincon(@(x) Qn_MR_function_additive(x,theta_H0,coordinate,f_ineq,f_eq,...
                        f_stdev_ineq,f_stdev_eq,G_ineq(:,r),G_eq(:,r),KMSoptions,kappa,p,k,2),...
                        starting_values_r(s,:),[],[],[],[],lbi,ubi,[],options);

                    if Qn_aux  < min_value_PR && bandera >= 1
                        minimizer = theta_aux;
                        min_value_PR = Qn_aux;
                    end
                end

                % compute simulated PR criterion function
                minQn_PR(1,r) = min_value_PR ;

                %- Step 3: combine PR and DR to get MR
                minQn_MR(1,r) = min(minQn_DR(1,r),minQn_PR(1,r));
            end

        else
            % Initialization of simulator %
            stream=RandStream('mlfg6331_64','Seed',seed);
            RandStream.setGlobalStream(stream);
            W2_AA = randn(R,n); % Random draws for RC, RS, and MR tests
            % CV calculation
            parfor r = 1:R
                rng(r,'twister');
                %- Step 1: Discard resampling method
                % Pick a large value as initial function value for the minimization;
                min_value_DR = 10^10;

                minQn_DR_aux = [];
                % Minimize but restricted to points of the estimated indentified set;
                for IDsetHat_index = 1:size(In_identified_set_hat,1)
                    minQn_DR_aux = min(min_value_DR,Qn_MR_function(In_identified_set_hat(IDsetHat_index,:),...
                        theta_H0,coordinate,data,kappa,p,k,W2_AA(r,:),1,DGP));
                end

                % compute simulated DR criterion function
                minQn_DR(1,r) = minQn_DR_aux;

                %- Step 2: Penalized resampling method;
                % Starting values for PR minimization;
                starting_values_r =  uniquetol2([starting_values;Qn_minimizer],tolerance,'rows');

                % Pick a large value as initial function value for the minimization;
                min_value_PR = 10^10;
                % check whether minimization is successful and reduced value
                for s=1:size(starting_values_r,1)
                    [theta_aux,Qn_aux,bandera] =  fmincon(@(x) Qn_MR_function(x,theta_H0,coordinate,data,kappa,...
                        p,k,W2_AA(r,:),2,DGP),starting_values_r(s,:),[],[],[],[],lbi,ubi,[],options);
                    if Qn_aux  < min_value_PR && bandera >= 1
                        minimizer = theta_aux;
                        min_value_PR = Qn_aux;
                    end
                end

                % compute simulated PR criterion function
                minQn_PR(1,r) = min_value_PR ;

                %- Step 3: combine PR and DR to get MR
                minQn_MR(1,r) = min(minQn_DR(1,r),minQn_PR(1,r));
            end
        end
        c_val = quantile(minQn_MR, 1-alpha);  
    end
    
end
end
