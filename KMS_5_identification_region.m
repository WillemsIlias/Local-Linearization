function [Identification_region] =  KMS_5_identification_region(theta_true,theta_0,LB_theta,UB_theta,A_theta,b_theta,KMSoptions)
%% Code Description:
% This function computes the true identification region for the games
% and the BCS example.

%% Extract information
DGP = KMSoptions.DGP;
dim_p = size(theta_true,1);
KMSoptions.dim_p = dim_p;
Identification_region = zeros(dim_p,2);
options_fmincon = KMSoptions.options_fmincon;
options_multistart  = KMSoptions.options_multistart;

%% Population moments
if DGP == 20 
    % True population moment:
    % PRELIMINARY
    dX = KMSoptions.dX;
    psuppX = KMSoptions.psuppX;
    selp = KMSoptions.selp;
    J1 = 2*dX;
    J2 = 2*dX;
    J = J1 + 2*J2;
    KMSoptions.J1 = J1;
    KMSoptions.J2 = J2;
    KMSoptions.J = J;

    % Preset output variables
    f_ineq_pop = zeros(J1,1);                        % Preset moment inequalities and
    f_eq_pop = zeros(J2,1);                         % moment equalities

    % Extract parameters (easier to read)
    beta1 = theta_true(1:dX,1);  
    beta2 = theta_true(dX+1:2*dX,1);
    delta1 = theta_true(2*dX+1); 
    delta2 = theta_true(2*dX+2);
    rho = theta_true(2*dX+3);

    % MOMENT COMPUTATION
    for ii = 1:dX
        p_ub10 = mexBVNcdf([beta1(ii),-beta2(ii)-delta2],[0 0],[1,-rho;-rho,1]);
        p_lb10 = mexBVNcdf([(beta1(ii)+delta1),-(beta2(ii)+delta2)],[0 0],[1,-rho;-rho,1])...
               +mexBVNcdf([-beta1(ii)-delta1,-beta2(ii)],[0 0],[1,rho;rho,1]) ...
               -mexBVNcdf([-beta1(ii),-beta2(ii)],[0 0],[1,rho;rho,1]);
        p_m10 = p_ub10 - p_lb10;
        p_10 = p_lb10 + p_m10*(1-selp(ii));
        
        
        f_ineq_pop((ii-1)*2 + 1,1) =  p_10*psuppX(ii);                                   % Moments m_{2q} in Eq 5.1 BCS
        f_ineq_pop((ii-1)*2 + 2,1) = -p_10*psuppX(ii);

        f_eq_pop((ii-1)*2 + 1,1)   = mexBVNcdf([-beta1(ii),-beta2(ii)],[0 0],[1,rho;rho,1])*psuppX(ii);
        f_eq_pop((ii-1)*2 + 2,1)   = mexBVNcdf([(beta1(ii)+delta1),(beta2(ii)+delta2)],[0 0],[1,rho;rho,1])*psuppX(ii);
    end

    % Concatenate the positive and negative of f_eq to transform into moment
    % inequalities
    f_eq_pop = [f_eq_pop ; -f_eq_pop];
end


%% Identification region
% We solve min/max p'theta in all basis vector directions
options_fmincon.TolCon = 1e-10;
options_fmincon.TolFun = 1e-10;
options_fmincon.TolX = 1e-10;

options_multistart.Display = 'off';
mult_num = 50;

for ii = 1:dim_p
   fprintf('Starting Identification Region Component=%d \n',ii)

   % maximize direction
   p = zeros(dim_p,1);
   p(ii) = 1;
   
   % Objective function and constraint:
   obj_true =  @(theta)KMS_51_IRobjective(theta,p);
   constraint_true =  @(theta)KMS_52_IRconstraint(theta,f_ineq_pop,f_eq_pop,KMSoptions);
    
   % Solve using fmincon from each initial point theta_0_fminimax. 
   %  NB: we select the scale option in fmincon to avoid scaling issues.
  % [x,fval,exitflag] = fmincon(obj_true,theta_true,A_theta,b_theta,[],[],...
  %                  LB_theta,UB_theta,constraint_true,options_fmincon);

   problem = createOptimProblem('fmincon','x0',theta_0,...
        'objective',obj_true,'Aineq', A_theta, 'bineq', b_theta,'lb',LB_theta,'ub',UB_theta,...
        'nonlcon',constraint_true,'options',options_fmincon);
       
    [~,fval,exitflag] = run(options_multistart,problem,mult_num);
    if exitflag > 0
        Identification_region(ii,2) = -fval;
    else
        Identification_region(ii,2) = nan;
    end

   % minimize direction
   p = zeros(dim_p,1);
   p(ii) = -1;
   
   % Objective function and constraint:
   obj_true =  @(theta)KMS_51_IRobjective(theta,p);
   constraint_true =  @(theta)KMS_52_IRconstraint(theta,f_ineq_pop,f_eq_pop,KMSoptions);
   
   problem = createOptimProblem('fmincon','x0',theta_0,...
        'objective',obj_true,'Aineq', A_theta, 'bineq', b_theta,'lb',LB_theta,'ub',UB_theta,...
        'nonlcon',constraint_true,'options',options_fmincon);

    [~,fval,exitflag] = run(options_multistart,problem,mult_num);
    if exitflag > 0
        Identification_region(ii,1) = fval;
    else
        Identification_region(ii,1) = nan;
    end
end



end