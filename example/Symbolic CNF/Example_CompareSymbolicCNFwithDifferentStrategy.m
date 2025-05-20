% This m file provide a example to demonstrate how to use direct parametrisation method for 
% invariant manifold (DPIM) with BP (appling Bloch's assumption onto the physic model) and 
% PN (imposing the periodic assumption to the normal coordinates of CNF) strategy to obtain 
% the symbolic CNFs. 
% 
% The concerned dynamic model is a 2 DOFs locally resonant metamaterial lattice, with expression as
%
% \mathbf{M}\ddot{\mathbf{x}}_n + \sum_{i=-1}^{1}{\mathbf{C}_i\dot{\mathbf{x}}_{n+i}}
% + \sum_{i=-1}^{1}{\mathbf{K}_i\mathbf{x}_{n+i}} +\mathbf{F}_n\left(\mathbf{x}_{n}\right)
% +\mathbf{F}_{nb,-1}\left(\mathbf{x}_{n}-\mathbf{x}_{n-1}\right) 
% +\mathbf{F}_{nb,+1}\left(\mathbf{x}_{n}-\mathbf{x}_{n+1}\right) = 0
%
% where
% \mathbf{M} = \begin{bmatrix}m_h+m_a & m_a\\ m_a & m_a \end{bmatrix},
% \mathbf{K}_{\pm1} = \begin{bmatrix}-\kappa_1 & 0\\  0 & 0 \end{bmatrix},
% \mathbf{K}_0 = \begin{bmatrix}2\kappa_1 & 0\\ 0 & \kappa_2\end{bmatrix}, 
% \mathbf{C}_{\pm1} = \begin{bmatrix}-c_1 & 0\\0 & 0\end{bmatrix},
% \mathbf{C}_0 = \begin{bmatrix} 2c_1 & 0\\ 0 & c_2\end{bmatrix}
% \mathbf{x}_{n+i} = \begin{bmatrix}u_{n+i}\\v_{n+i}\end{bmatrix}, 
% \mathbf{F}_{nb,\pm1} = \begin{bmatrix}
%     \mp\beta_1 (u_n-u_{n\pm1})^2 + \gamma_1 (u_n-u_{n\pm1})^3\\
%     0
% \end{bmatrix}, \mathbf{F}_{n} = \begin{bmatrix}
%     0\\
%     \beta_2 v_n^2 + \gamma_2 v_n^3
% \end{bmatrix}.

%% clear date
clear
close all
clc

%% Set dynamic's parameters
syms k beta mh ma kappa1 kappa2 gamma2 gamma1 beta1 beta2 c1 c2 real
parameters_model = struct('mh',1,'ma',0.5,'kappa1',1,'beta1',0,'gamma1',0,...
    'k',k,'kappa2',0.4,'beta2',beta2,'gamma2',0,'c1',0,'c2',0);

%% Using the DPIM to calculatic the symbolic CNF-PN and CNF-BP
mode = 1; % concerned dispersion branch, 1 for acoustic branch, 2 for optic branch
ReducedModel_CNF_PN = ObtainDispersion_Analytic_CNF_PN(parameters_model,mode);
ReducedModel_CNF_BP = ObtainDispersion_Analytic_CNF_BP(parameters_model);

%% Output the coefficients of z_1^2*z_2 in the CNF
disp('The expression of \varpi_1 of CNF-PN is as:')
disp(latex(ReducedModel_CNF_PN.f{3}(1,2)))
disp('The expression of \varpi_1 of CNF-BP is as:')
disp(latex(ReducedModel_CNF_BP.f{3}(1,2)))