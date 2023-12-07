%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Defines the sample criterion function for DR or PR inference
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = Qn_MR_function_additive(theta_to_min,theta_H0,coordinate,f_ineq,f_eq,f_stdev_ineq,f_stdev_eq,G_ineq,G_eq,KMSoptions,kappa,p,k,MR_type)
% theta_to_min denotes the subvector of theta that should be minimized.
% theta_H0 denotes the subvector of theta that is fixed according to H0.
% coodinate denotes the coordinate of interest (1 or 2).
% data denotes the dataset.
% kappa denotes the tuning parameter using by DR or PR.
% p denotes the number of moment inequalities (which should appear first).
% k denotes the number of moment (in)equalities.
% W2_AA is a vector of random variables used to implement the (multiplier) bootstrap.
% MR_type indicates the type of resampling, i.e., DR or PR.

theta = [theta_to_min(1:coordinate-1),theta_H0,theta_to_min(coordinate:end)];
J1 = KMSoptions.J1;
J2 = KMSoptions.J2;
n = KMSoptions.n;
[g_ineq,g_eq] = moments_theta(theta,J1,J2,KMSoptions);
xi = [-f_ineq' - g_ineq', f_eq(1:J2)' + g_eq(1:J2)']./[f_stdev_ineq',f_stdev_eq(1:J2)']*sqrt(n)/kappa;

if MR_type == 1 % DR method;
    value = S_function([-G_ineq' G_eq(1:J2)'] + repmat(phi_function(xi,p,k-p),size(G_ineq',1),1),p);
else % PR method;
    value = S_function([-G_ineq' G_eq(1:J2)'] + repmat(xi,size(G_ineq',1),1),p);
end
