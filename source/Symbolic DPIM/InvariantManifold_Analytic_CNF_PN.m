function ReducedModel = InvariantManifold_Analytic_CNF_PN(Model, mode)
% Making use of the DPIM to calculate the CNF-PN and corresponding
% nonlinear maping or invariant manifold

%% Set up the parameter
n_p = length(Model.Indp);% Truncation order of monomial

M = Model.M;
C0 = Model.Cn{1};
Cb = Model.Cn{2};
K0 = Model.Kn{1};
Kb = Model.Kn{2};
Cx = @(kx) C0 + 2*cos(kx)*Cb;
Kx = @(kx) K0 + 2*cos(kx)*Kb;
k = Model.k;
Lambda = diag(Model.spectrum.L);

Fn = Model.nonlinear_force;
Fn_b = Model.nonlinear_force_b;

f = cell(n_p,1);
Psi = cell(n_p,1);
Upsilon = cell(n_p,1);
Indp = Model.Indp;
n = 2;
N = Model.n;
Y = Model.spectrum.Y;

%% Build the first-order CNF and invariant manifold using the linear dispersion spectrum
f{1} = Model.spectrum.L;                                            % First-order CNF
Psi{1} = Model.spectrum.Y(1:N,:);                             % First-order invariant manifold of displacement variables
Upsilon{1} = Model.spectrum.Y(N+1:end,:);             % First-order invariant manifold of velocity variables

%% Build the p-order CNF and invariant manifold using DPIM framework
for p_i = 2:n_p
I_pi = Indp{p_i};
n_Ipi = size(I_pi,2);
f{p_i} = sym(sparse(n, n_Ipi));
Psi{p_i} = sym(sparse(N, n_Ipi));
Upsilon{p_i} = sym(sparse(N, n_Ipi));

for l_i = 1:n_Ipi
    I = I_pi(:, l_i);
    theta_I = I'*Lambda;

    mu_I = getRInd(Psi, f, Indp, I, 0, k);
    mu_Im1 = getRInd(Psi, f, Indp, I, -1, k);
    mu_Ip1 = getRInd(Psi, f, Indp, I, 1, k);
    nu_I = getRInd(Upsilon, f, Indp, I, 0, k);
    G_I = getFnI(Fn, Psi, Indp, I, 2, 0, k);
    G_Im1 = getFnI(Fn_b, Psi, Indp, I, 2, -1, k);
    G_Ip1 = getFnI(Fn_b, Psi, Indp, I, 2, 1, k);
    H_I = getFnI(Fn, Psi, Indp, I, 3, 0, k);
    H_Im1 = getFnI(Fn_b, Psi, Indp, I, 3, -1, k);
    H_Ip1 = getFnI(Fn_b, Psi, Indp, I, 3, 1, k);
    Xi_I = simplify(- G_I-G_Im1-G_Ip1 - H_I- H_Im1- H_Ip1 - M*nu_I ...
        - (theta_I*M+C0)*mu_I-Cb*mu_Im1-Cb*mu_Ip1);

    na = I(2)-I(1);
    
    switch abs(na)
        case 0 & (mode == 1)
            % acoustic resonance of alpha(p,q) = [n,n]
            V1 = Psi{1}(:,1); V2 = Psi{1}(:,2); Vr = [V1,V2];
            lambda1 = f{1}(1,1); lambda2 = f{1}(2,2);
            A = [Kx(0)+theta_I*Cx(0)+theta_I^2*M,(Cx(k)+(lambda1+theta_I)*M)*V1,(Cx(k)+(lambda2+theta_I)*M)*V2;
                zeros(1,2),1,1;
                1,0,0,0];
            R = [Xi_I;0;0];
            A = A(1:3,2:end);
            R = R(1:3);
            Phf = A\R;
            Psi{p_i}(2:end,l_i) = simplify(Phf(1:N-1),'Steps',30);
            f{p_i}(1:2,l_i) = simplify(Phf(N:end),'Steps',30);
            Upsilon{p_i}(:,l_i) = simplify(theta_I*Phf(1:N-1) + Vr*Phf(N:end) + mu_I,'Steps',30);

        case 1
            % trivial resonance of alpha(p,q) = [n+1,n] or [n,n+1]
            IR = (na == 1) + 1; % na=1 for mode 2，otherwise mode 1
            Vr = Y(1:N,IR);
            A11 = theta_I^2*M + theta_I*Cx(k) + Kx(k);
            A12 = ((theta_I+f{1}(IR,IR))*M + Cx(k))*Vr; 
            A21 = Vr.'*((theta_I + f{1}(IR,IR))*M +Cx(k));
            A22 = Vr.'*M*Vr;
            A = [A11 A12;A21 A22];
            B = [Xi_I;  -Vr.'*M*mu_I];
            Phf = A\B;
            Psi{p_i}(:,l_i) = Phf(1:N);
            f{p_i}(IR,l_i) = Vr.'*Xi_I;
            Upsilon{p_i}(:,l_i) = theta_I*Phf(1:N) + Vr*Vr.'*Xi_I + mu_I;
        otherwise
            % non-resonant case
            A = theta_I^2*M + theta_I*Cx(na*k) + Kx(na*k);
            Phf = A\Xi_I;
            Psi{p_i}(:,l_i) = simplify(Phf,'Steps',30);
            Upsilon{p_i}(:,l_i) = simplify(theta_I*Phf + mu_I,'Steps',30);
    end
    disp(['Completed the calculation of monomial of alpha(p,q) = ', mat2str(I)])
end
end

%% Output the results
ReducedModel.f = f;
ReducedModel.Psi = Psi;
ReducedModel.Upsilon = Upsilon;
ReducedModel.I = Indp;
end