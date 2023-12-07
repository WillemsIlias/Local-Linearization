% csolve  Solves a custom quadratic program very rapidly.
%
% [vars, status] = csolve(params, settings)
%
% solves the convex optimization problem
%
%   minimize((1/2)*quad_form(x, H) + f'*x)
%   subject to
%     x <= ub
%     x(1) >= lb(1)
%     x(2) >= lb(2)
%     x(3) >= lb(3)
%     x(4) >= lb(4)
%     x(5) >= lb(5)
%     x(6) >= lb(6)
%     x(7) >= lb(7)
%     x(8) >= lb(8)
%     x(9) >= lb(9)
%     x(10) >= lb(10)
%     x(11) >= lb(11)
%     x(12) >= lb(12)
%     x(13) >= lb(13)
%     x(14) >= lb(14)
%     x(15) >= lb(15)
%     x(16) >= lb(16)
%     x(17) >= lb(17)
%     x(18) >= lb(18)
%
% with variables
%        x  34 x 1
%
% and parameters
%        H  34 x 34   PSD
%        f  34 x 1
%       lb  34 x 1
%       ub  34 x 1
%
% Note:
%   - Check status.converged, which will be 1 if optimization succeeded.
%   - You don't have to specify settings if you don't want to.
%   - To hide output, use settings.verbose = 0.
%   - To change iterations, use settings.max_iters = 20.
%   - You may wish to compare with cvxsolve to check the solver is correct.
%
% Specify params.H, ..., params.ub, then run
%   [vars, status] = csolve(params, settings)
% Produced by CVXGEN, 2022-06-01 01:00:28 -0400.
% CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com.
% The code in this file is Copyright (C) 2006-2017 Jacob Mattingley.
% CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
% applications without prior written permission from Jacob Mattingley.

% Filename: csolve.m.
% Description: Help file for the Matlab solver interface.
