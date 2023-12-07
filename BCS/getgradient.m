function G = getgradient(theta)
%   GRADIENT COMPUTATION
%   Detailed explanation goes here

% Extract parameters (easier to read)
theta = reshape(theta,max(size(theta)),1);
dimSX = (size(theta,1) -3)/2; % dimension of X.
beta1 = theta(1:dimSX,1);  beta2 = theta(dimSX+1:2*dimSX,1);
delta1 = theta(2*dimSX+1); delta2 = theta(2*dimSX+2);
rho = theta(2*dimSX+3);

mu   = [0,0];
sigma_rho = [1,rho;rho,1];
r= -rho;
sigma_r = [1,r; r,1];
psuppX = ones(dimSX,1)/dimSX; % P(X=x), X is discrete uniform.
for ii = 1:dimSX
    % Pick support point (include constant)
    x = zeros(1,dimSX);
    x(ii) = 1;

    %%%% Gradient WRT beta1,beta2,delta1,delta2 %%%%%%%%%%%%%%%%%%%%%%%
    % Gradients for moment inequalities for entry game (See Pg 34, eq
    % 5.3)(0,1)
    % NB: g3 has r correlation. Derivatives wrt beta should not change
    % its form much
%     Dg3b1c = x1.*normpdf(-x1*beta1-delta1) .*normcdf((x2*beta2-r*(-x1*beta1-delta1))/sqrt(1-r^2))*pX;
%     Dg3b2c = -x2.*normpdf(x2*beta2).*normcdf((-x1*beta1-delta1-r*(x2*beta2))/sqrt(1-r^2))*pX;
%     Dg3d1c = normpdf(-x1*beta1-delta1) .*normcdf((x2*beta2-r*(-x1*beta1-delta1))/sqrt(1-r^2))*pX;
%     Dg3d2c = zeros(1,1);
    
    Dg3b1c = x.*normpdf(x*beta1).*normcdf((-(x*beta2+delta2)-r*(x*beta1))/sqrt(1-r^2))*psuppX(ii);
    Dg3b2c = -x.*normpdf(-x*beta2-delta2).*normcdf((x*beta1-r*(-x*beta2-delta2))/sqrt(1-r^2))*psuppX(ii);
    Dg3d1c = zeros(1,1);
    Dg3d2c = -normpdf(-x*beta2-delta2).*normcdf((x*beta1-r*(-x*beta2-delta2))/sqrt(1-r^2))*psuppX(ii);
      
    Dg4b1c = -x.*normpdf(x*beta1+delta1).*normcdf((-(x*beta2+delta2)-r*(x*beta1+delta1))/sqrt(1-r^2))*psuppX(ii) ...
        +(x.*normpdf(-x*beta1-delta1).*normcdf((-x*beta2-rho*(-x*beta1-delta1))/sqrt(1-rho^2)))*psuppX(ii)...
        -(x.*normpdf(-x*beta1).*normcdf((-x*beta2-rho*(-x*beta1))/sqrt(1-rho^2)))*psuppX(ii);
    
    Dg4b2c = x.*normpdf(x*beta2+delta2).*normcdf((x*beta1+delta1-r*(-x*beta2-delta2))/sqrt(1-r^2))*psuppX(ii) ...
        +(x.*normpdf(-x*beta2).*normcdf((-x*beta1-delta1-rho*(-x*beta2))/sqrt(1-rho^2)))*psuppX(ii)...
        - x.*normpdf(-x*beta2).*normcdf((-x*beta1-rho*(-x*beta2))/sqrt(1-rho^2))*psuppX(ii);
    
    Dg4d1c = -normpdf(x*beta1+delta1).*normcdf((-(x*beta2+delta2)-r*(x*beta1+delta1))/sqrt(1-r^2))*psuppX(ii) ...
        +(normpdf(-x*beta1-delta1).*normcdf((-x*beta2-rho*(-x*beta1-delta1))/sqrt(1-rho^2)))*psuppX(ii);
    Dg4d2c = normpdf(x*beta2+delta2).*normcdf((x*beta1+delta1-r*(-x*beta2-delta2))/sqrt(1-r^2))*psuppX(ii);

    % Gradients for moment equalities for entry game (See Pg 34, eq 5.1) (0,0)
    Dg1b1c = x.*normpdf(-x*beta1).*normcdf((-x*beta2-rho*(-x*beta1))/sqrt(1-rho^2))*psuppX(ii);
    Dg1b2c = x.*normpdf(-x*beta2).*normcdf((-x*beta1-rho*(-x*beta2))/sqrt(1-rho^2))*psuppX(ii);
    Dg1d1c = zeros(1,1);
    Dg1d2c = zeros(1,1);

    % Gradients for moment equalities for entry game (See Pg 34, eq 5.2) (1,1)
    Dg2b1c = -x.*normpdf(x*beta1+delta1).*(normcdf((x*beta2+delta2-rho*(x*beta1+delta1))/sqrt(1-rho^2)))*psuppX(ii);
    Dg2b2c = -x.*normpdf(x*beta2+delta2).*(normcdf((x*beta1+delta1-rho*(x*beta2+delta2))/sqrt(1-rho^2)))*psuppX(ii);
    Dg2d1c = -normpdf(x*beta1+delta1).*(normcdf((x*beta2+delta2-rho*(x*beta1+delta1))/sqrt(1-rho^2)))*psuppX(ii);
    Dg2d2c = -normpdf(x*beta2+delta2).*(normcdf((x*beta1+delta1-rho*(x*beta2+delta2))/sqrt(1-rho^2)))*psuppX(ii);

    %%%% Gradient WRT rho %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gradients for moments inequalities w.r.t. rho
    sgn_r   = -1;
    upper3 = [-x*beta2-delta2,x*beta1];
    Dg3rho = sgn_r*get_drho(r,mu,sigma_r,upper3)*psuppX(ii);

    % Gradients for moments inequalities w.r.t. rho
    upper4_1 = [-x*beta2-delta2,x*beta1+delta1];
    upper4_2 = [-x*beta2,-x*beta1-delta1];
    upper4_3 = [-x*beta2,-x*beta1];
    Dg4rho = -sgn_r*get_drho(r,mu,sigma_r,upper4_1)*psuppX(ii)...
        -  get_drho(rho,mu,sigma_rho,upper4_2)*psuppX(ii)...
        +  get_drho(rho,mu,sigma_rho,upper4_3)*psuppX(ii);

    upper1 = [-x*beta1,-x*beta2];
    Dg1rho = -get_drho(rho,mu,sigma_rho,upper1)*psuppX(ii);

    upper2 = [x*beta1+delta1,x*beta2+delta2];
    Dg2rho = -get_drho(rho,mu,sigma_rho,upper2)*psuppX(ii);

    %%%% Compile Moments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dg_ineq((ii-1)*2 + 1,:) = [Dg3b1c Dg3b2c Dg3d1c Dg3d2c Dg3rho];
    Dg_ineq((ii-1)*2 + 2,:) = [Dg4b1c Dg4b2c Dg4d1c Dg4d2c Dg4rho];
    Dg_eq((ii-1)*2 + 1,:)   = [Dg1b1c Dg1b2c Dg1d1c Dg1d2c Dg1rho];
    Dg_eq((ii-1)*2 + 2,:)   = [Dg2b1c Dg2b2c Dg2d1c Dg2d2c Dg2rho];

end

G = [Dg_ineq; Dg_eq];
end

function drho = get_drho(rho,mu,sigma_rho,upper)
c1 = mexBVNcdf(upper,mu,sigma_rho)./((1-rho.^2).^2);
c2 = rho - rho.^3;
c3 = (1+rho.^2).*get_EXiXj(1,2,mu,sigma_rho,upper);
c4 = -rho.*(get_EXiXj(1,1,mu,sigma_rho,upper)+get_EXiXj(2,2,mu,sigma_rho,upper));
drho = c1.*(c2+c3+c4);
end

function EXiXj = get_EXiXj(i,j,mu,sigma_rho,upper)
if i==1 && j==2
    EXiXj = sigma_rho(1,2) ...
        + sigma_rho(2,1).*(-upper(:,1).*get_Fk(upper(:,1),1,mu,sigma_rho,upper))...
        + sigma_rho(1,2).*(-upper(:,2).*get_Fk(upper(:,2),2,mu,sigma_rho,upper))...
        + (sigma_rho(1,1).*sigma_rho(2,2)-sigma_rho(1,2).*sigma_rho(2,1))...
        .*get_Fqr(upper,mu,sigma_rho,upper);
elseif i==1 && j==1
    EXiXj= sigma_rho(1,1)...
        + sigma_rho(1,1).*(-upper(:,1).*get_Fk(upper(:,1),1,mu,sigma_rho,upper))...
        + sigma_rho(1,2).^2./sigma_rho(2,2)...
        .*(-upper(:,2).*get_Fk(upper(:,2),2,mu,sigma_rho,upper))...
        + (sigma_rho(1,2).*sigma_rho(1,1)-sigma_rho(1,2).^2.*sigma_rho(2,1)./sigma_rho(2,2))...
        .*get_Fqr(upper,mu,sigma_rho,upper);
elseif i==2 && j==2
    EXiXj= sigma_rho(2,2)...
        + sigma_rho(2,2).*(-upper(:,2).*get_Fk(upper(:,2),2,mu,sigma_rho,upper))...
        + sigma_rho(2,1).^2./sigma_rho(1,1)...
        .*(-upper(:,1).*get_Fk(upper(:,1),1,mu,sigma_rho,upper))...
        + (sigma_rho(2,1).*sigma_rho(2,2)-sigma_rho(2,1).^2.*sigma_rho(1,2)./sigma_rho(1,1))...
        .*get_Fqr(upper,mu,sigma_rho,upper);
end
end

function Fqr = get_Fqr(X,mu,sigma_rho,upper)
% X: n-by-2
% upper: n-by-2
% mean: 1-by-2
% sigma: 2-by-2
alpha = mexBVNcdf(upper,mu,sigma_rho);
Fqr = mvnpdf(X,mu,sigma_rho)./alpha;
end

function Fk = get_Fk(xn,i,mu,sigma_rho,upper)
n = length(xn);
C = sigma_rho;
A = inv(sigma_rho);
if i==1
    j=2;
elseif i==2
    j=1;
end
A_1 = A(j,j);
A_1_inv = inv(A_1);
C_1 = C(j,j);
c_nn = C(i,i);
c = C(j,i);
mu_1 = mu(j);
mu_n = mu(i);
f_xn = zeros(n,1);
p = mexBVNcdf(upper,mu,sigma_rho);
for l=1:n
    m = mu_1 + (xn(l) - mu_n) .* c/c_nn;
    f_xn(l) = exp(-0.5 .* (xn(l) - mu_n).^2./c_nn) .* normcdf(upper(l,j),m,sqrt(A_1_inv));
end
Fk = 1./p .* 1./sqrt(2 * pi * c_nn) .* f_xn;
end
