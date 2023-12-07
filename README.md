# Local-Linearization
Code for Bei (JoE, 2024)

I combine my test and the E-A-M algorithm in Kaido, Molinari, and Stoye (2019) to calculate the confidence interval. The command is at BCS_LL_KMS_Simulation_1k_dX2.m line 124-125:

[LL_confidence_interval,LL_output] = KMS_0_Main(W,theta_0,...
        p,[], LB_theta, UB_theta, A_theta, b_theta, alpha, type, method, kappa, phi, CVXGEN_name, KMSoptions)
        
The test statistic and critical value at a single point are given in /BCS/BCS_32_Critval.m (with LL_EAM = 1 at line 87). 
 
I use CVXGEN for quadratic programming, and you can download their code here:
https://cvxgen.com/docs/index.html
with the problem setup in 00CVXGEN.txt, where dim_p is the dimension of parameter theta, J1 is the number of moment inequalities and J2 is the number of moment equalities.
