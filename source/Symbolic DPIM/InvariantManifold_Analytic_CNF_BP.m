function ReducedModel = InvariantManifold_Analytic_CNF_BP(Model)
% Making use of the DPIM to calculate the CNF-BP and corresponding
% nonlinear maping or invariant manifold

%% Set up the parameter
n_p = length(Model.Indp);

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
N = Model.n;
n = 2;
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

    mu_I = getRInd(Psi, f, Indp, I, 0);
    nu_I = getRInd(Upsilon, f, Indp, I, 0);
    G_I = getFnI(Fn, Psi, Indp, I, 2, 0);
    G_Ib = ((1-exp(1i*k))^2-(1-exp(-1i*k))^2)*getFnI(Fn_b, Psi, Indp, I, 2, 0);
    H_I = getFnI(Fn, Psi, Indp, I, 3, 0);
    H_Ib = ((1-exp(1i*k))^3+(1-exp(-1i*k))^3)*getFnI(Fn_b, Psi, Indp, I, 3, 0);
    Xi_I = - G_I - H_I -G_Ib -H_Ib- M*nu_I - (theta_I*M + Cx(k))*mu_I;

    na = I(2)-I(1);
    
    switch abs(na)
        case 1
            % trivial resonance of alpha(p,q) = [n+1,n] or [n,n+1]
            IR = (na == 1) + 1; % na=1 for mode 2，otherwise mode 1
            Vr = Y(1:N,IR);
            A11 = theta_I^2*M + theta_I*Cx(k) + Kx(k);
            A12 = ((theta_I+f{1}(IR,IR))*M + Cx(k))*Vr; 
            A21 = Vr.'*((theta_I + f{1}(IR,IR))*M +Cx(k));
            A22 = Vr.'*M*Vr;
            A = [A11 A12;A21 A22];
            B = [Xi_I; - Vr.'*M*mu_I];
            Phf = A\B;
            Psi{p_i}(:,l_i) = Phf(1:N);
            f{p_i}(IR,l_i) = Vr.'*Xi_I;
            Upsilon{p_i}(:,l_i) = theta_I*Phf(1:N) + Vr*Vr.'*Xi_I + mu_I;
        otherwise
            % non-resonant case
            A = theta_I^2*M + theta_I*Cx(k) + Kx(k);
            Phf = A\Xi_I;
            Psi{p_i}(:,l_i) = simplify(Phf,'Steps',30);
            Upsilon{p_i}(:,l_i) = simplify(theta_I*Phf + mu_I,'Steps',30);
    end
    disp(['Completed the calculation of monomial of alpha(p,q) = ', mat2str(I)])
end
end
ReducedModel.f = f;
ReducedModel.Psi = Psi;
ReducedModel.Upsilon = Upsilon;
ReducedModel.I = Indp;
end
