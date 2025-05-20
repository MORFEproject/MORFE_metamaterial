function Model = Build_cell_dpim(parameter,k)
% This function to build the dynamic model for the following calculation of
% DPIM, the concerned dynamic model is a 2 DOFs locally resonant metamaterial lattice, with expression as
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

%% setting model's parameter
mh = parameter.mh;                              % host oscillator's mass
ma = parameter.ma;                              % attachment's mass
kappa1 = parameter.kappa1;                 % host oscillator's stiffness
kappa2 = parameter.kappa2;                 % attachment's stiffness
beta1 = parameter.beta1;                       % host oscillator's quadratic nonlinear stiffness
beta2 = parameter.beta2;                        % attachment's quadratic nonlinear stiffness
gamma1 = parameter.gamma1;              % host oscillator's cubic nonlinear stiffness
gamma2 = parameter.gamma2;              % attachment's cubic nonlinear stiffness
c1 = parameter.c1;                                  % host oscillator's linear damping
c2 = parameter.c2;                                  % attachment's linear damping

%% Build linear matrix
M = [mh+ma, ma;ma, ma]; K = diag([2*kappa1*(1-cos(k)), kappa2]);
C = diag([2*c1*(1-cos(k)), c2]);
Cn = {diag([2*c1, c2]), diag([-c1,0])};
Kn = {diag([2*kappa1, kappa2]),diag([-kappa1,0])};

%% Build nonlinear force
fnl = cell(3,1);
fnlb = cell(3,1);

if beta2 == 0
    fnl{2} = {};
else
    I2 = [2 2]';
    V2 = [0;beta2];
    fnl{2} = struct('I',I2,'vector',V2);
end

if gamma2 == 0
    fnl{3} = {};
else
    I3 = [2 2 2]';
    V3 = [0;gamma2];
    fnl{3} = struct('I',I3,'vector',V3);
end

if beta1 == 0
    fnlb{2} = {};
else
    I2b = [1 1]';
    V2b = [beta1;0];
    fnlb{2} = struct('I',I2b,'vector',V2b);
end

if gamma1 == 0
    fnlb{3} = {};
else
    I3b = [1 1 1]';
    V3b = [gamma1;0];
    fnlb{3} = struct('I',I3b,'vector',V3b);
end

%% Construct a struct with the model's information
Model = struct('M',M,'K',K,'C',C);
Model.Cn = Cn;
Model.Kn = Kn;
Model.k = k;
Model.n = size(M,2);
Model.nonlinear_force = fnl;
Model.nonlinear_force_b = fnlb;
end