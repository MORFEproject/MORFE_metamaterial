% Compare MSM, HBM, and NF-PN(9) in predicting the wave attenuation of 
% optic dispersion solution at k=0.8pi for complex nonlinear chain
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
clear
close all
clc

%% Set up the parameters
mh = 1;                     % host's mass
ma = 0.5;                   % attachment's mass
kappa1 = 1;               % linear stiffness in host oscillator
kappa2 = 0.4;            % linear stiffness in attachment
beta1 = 1;                  % quadratic nonlinear stiffness in host oscillator
beta2 = -1;                  % quadratic nonlinear stiffness in attachment
gamma1 = 2;             % cubic nonlinear stiffness in host oscillator
gamma2 = 3;             % cubic nonlinear stiffness in attachment
c1 = 5e-3;                       % linear damping in host oscillator
c2 = 2e-2;                       % linear damping in attachment

% The parameter of model
parameters_model = struct('mh',mh,'ma',ma,'kappa1',kappa1,...
    'kappa2',kappa2,'beta1',beta1,'beta2',beta2,'gamma1',gamma1,'gamma2',gamma2,'c1',c1,'c2',c2);

%% Calculate the dispersion solutions with HBM, CNF, MMS
kn = 0.8*pi;
order1 = 3; %DPIM's order
order2 = 9; %DPIM's order
H = 5; % Order of harmonic truncation

% The dispersion solution of acoustic branch
log10a_s = -4.5; log10a_e = -0.4; mode = 1;
[omega{1}, nc{1}, umax{1}, X_HB, manifold{1}] = ObtainDispersion_HBM(parameters_model,H,mode,log10a_s,log10a_e, kn);
manifold{1} = ObtainManifold(manifold{1},[1 2 3 4]);
rho = 10.^(linspace(-4,-0.4,2e2));
[omega{2}, nc{2}, umax{2},manifold{2}, fun_pm] = ObtainDispersion_MMS(parameters_model, mode, rho, kn);
manifold{2} = ObtainManifold(manifold{2},[1 2 3 4]);
rho = 10.^(linspace(-4,-0.4,2e2));
[omega{3}, nc{3}, umax{3}, backbone, manifold{3}] = ObtainDispersion_CNF_PN(parameters_model, mode, rho, kn, order2);
manifold{3} = ObtainManifold(manifold{3},[1 2 3 4]);

figure 
hold on
box on
grid on
ni = 1;
plot(nc{1}{ni},log10(umax{1}{ni}(:,1)),'-k')
plot(nc{2}{ni},log10(umax{2}{ni}(:,1)),'--r')
plot(nc{3}{ni},log10(umax{3}{ni}(:,1)),'--b')
ylim([-3, -1.4])
xlim([0.017 0.019])
ylabel('Amplitude of $\log_{10}{u_n}$','FontName','Times New Roman','FontSize',10,'Interpreter','latex')
xlabel('Damping ratio {\it\zeta_n}','FontName','Times New Roman','FontSize',10)
set(gca,'FontName','Times New Roman','FontSize',10);
set(gcf, 'Position', [600, 600, 200, 170]);

%% Build the dynamic model of annular chain with 100 lattices
n_element = 100;
Model = build_NonlinearAnnularChain(n_element, parameters_model);

%% Choose the initial condition for numerical integration
Amplitude = 10^(-1.4);
[~, ind0] = min(abs(abs(umax{1}{1}(:,1))-Amplitude));
X_HB0 = X_HB{1}(:,ind0);
cell_ndof = 2;

[U0,  V0] = ObtainWavepocket(X_HB0, kn, n_element, cell_ndof, 0);

%% Obtain the wave attenuation process using the numerical integration
X0 = zeros(2*Model.ndof,1);
X0(1:Model.ndof) = U0; X0((Model.ndof+1):2*Model.ndof) = V0;

tspan = linspace(0,200,2e3);
[Tn,Xn] = ode113(@(t,X) Model.Fun(t,X), tspan, X0);

N_time = length(Tn);
Speedfactor = 50;

figure
for ni = 1:N_time/Speedfactor
    plot(1:n_element, Xn(Speedfactor*ni, 1:2:Model.ndof), '-r', ...
        1:n_element, Xn(Speedfactor*ni, 2:2:Model.ndof),'-b');
    xlabel('Number of elements','FontName','Times New Roman');
    ylabel('Displacement','FontName','Times New Roman');
    le = legend('{\itx}_m','{\itx}_n');
    set(le,'FontName','Times New Roman');
    title(['Snapshot at ', num2str(Tn(Speedfactor*ni)),' s'],'FontName','Times New Roman');
    xlim([1 n_element]);
    ylim(1.5*10^(X_HB0(end))*[-1 1]);
    pause(0.01);
end

%% Obtain the wave attenuation process based on the reduced-order model
% HBM's results
Tspan = [0, 200];
[Orbits_HBM, Tn_HBM] =  ObtainWaveAttenuation_HBM(X_HB{1}, omega{1}{1}, nc{1}{1}, X_HB0, Tspan);

% MMS's results
[Orbits_MMS, Tn_MMS] = ObtainWaveAttenuation_MMS(fun_pm, U0, V0, Tspan);

% CNF's results
[Orbits_CNF, Tn_CNF] = ObtainWaveAttenuation_CNF(backbone, U0, V0, Tspan);

%% Compare HBM, MMS, CNF in predicting the wave attenuation
figure
hold on
box on
grid on
plot(Tn,Xn(:,200+2*n_element-1),'-k','LineWidth',1)
plot(Tn_HBM,Orbits_HBM(3,:),'--','Color',[161, 71, 67]/255,'LineWidth',1)
plot(Tn_MMS,Orbits_MMS(3,:),'--b','LineWidth',1)
plot(Tn_CNF,Orbits_CNF(3,:),'--','Color',[1, 0.6, 0],'LineWidth',1)
le = legend('Numerical integration','EPMC-HBM','MMS','CNF-PN 9','FontName','Times New Roman','FontSize',10);
set(le,'box','off')
xlabel('Dimensionless time {\it\tau}','FontName','Times New Roman','FontSize',10)
ylabel('Velocity $\dot{u}_n$','FontName','Times New Roman','FontSize',10,'Interpreter','latex')
set(gca,'FontName','Times New Roman','FontSize',10);
set(gcf, 'Position', [600, 600, 300, 230]);

figure
hold on
box on
plot(Tn,Xn(:,200+2*n_element-1),'-k','LineWidth',1)
plot(Tn_HBM,Orbits_HBM(3,:),'--','Color',[161, 71, 67]/255,'LineWidth',1)
plot(Tn_MMS,Orbits_MMS(3,:),'--b','LineWidth',1)
plot(Tn_CNF,Orbits_CNF(3,:),'--','Color',[1, 0.6, 0],'LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',10);
set(gcf, 'Position', [600, 600, 300, 230]);

% Wave attenuation on the invariant manifold from HBM
figure
hold on
box on
grid on
threshold = 10^(-1.38); int = 2; ni = 1;
[~,ind_s1] = min(abs(max(abs(manifold{1}{ni}.out{1}),[],1)-threshold));
s1 = surf(manifold{1}{ni}.out{2}(:,1:ind_s1),manifold{1}{ni}.out{3}(:,1:ind_s1),...
    manifold{1}{ni}.out{4}(:,1:ind_s1), 'FaceColor',[0.5, 0.5, 0.5],'EdgeColor','none');
plot3(Xn(:,200), Xn(:,399), Xn(:,400),'-','Color',[161, 71, 67]/255,'LineWidth',2)
le = legend('Invariant Manifold','Numerical Orbits','FontName','Times New Roman','FontSize',10);
set(le,'box','off')
s1.FaceAlpha = 0.5;  
s1.LineWidth = 0.1; 
xlabel('$v_n$','FontName','Times New Roman','FontSize',10,'Interpreter','latex')
ylabel('$\dot{u}_n$','FontName','Times New Roman','FontSize',10,'Interpreter','latex')
zlabel('$\dot{v}_n$','FontName','Times New Roman','FontSize',10,'Interpreter','latex')
set(gca,'FontName','Times New Roman','FontSize',10);
set(gcf, 'Position', [600, 600, 300, 230]);
view([15 45])

% Wave attenuation on the invariant manifold from MMS
figure
hold on
box on
grid on
threshold = 10^(-1.3); int = 2; ni = 1;
[~,ind_s1] = min(abs(max(abs(manifold{2}{ni}.out{1}),[],1)-threshold));
s1 = surf(manifold{2}{ni}.out{2}(:,1:int:ind_s1),manifold{2}{ni}.out{3}(:,1:int:ind_s1),...
    manifold{2}{ni}.out{4}(:,1:int:ind_s1), 'FaceColor',[0.5, 0.5, 0.5],'EdgeColor','none');
plot3(Xn(:,200), Xn(:,399), Xn(:,400),'-','Color',[161, 71, 67]/255,'LineWidth',2)
s1.FaceAlpha = 0.5;  
s1.LineWidth = 0.1; 
xlabel('$v_n$','FontName','Times New Roman','FontSize',10,'Interpreter','latex')
ylabel('$\dot{u}_n$','FontName','Times New Roman','FontSize',10,'Interpreter','latex')
zlabel('$\dot{v}_n$','FontName','Times New Roman','FontSize',10,'Interpreter','latex')
set(gca,'FontName','Times New Roman','FontSize',10);
set(gcf, 'Position', [600, 600, 300, 230]);
view([15 45])

% Wave attenuation on the invariant manifold from CNF
figure
hold on
box on
grid on
threshold = 10^(-1.4); int = 2; ni = 1;
[~,ind_s1] = min(abs(max(abs(manifold{3}{ni}.out{1}),[],1)-threshold));
s1 = surf(manifold{3}{ni}.out{2}(:,1:int:ind_s1),manifold{3}{ni}.out{3}(:,1:int:ind_s1),...
    manifold{3}{ni}.out{4}(:,1:int:ind_s1), 'FaceColor',[0.5, 0.5, 0.5],'EdgeColor','none');
plot3(Xn(:,200), Xn(:,399), Xn(:,400),'-','Color',[161, 71, 67]/255,'LineWidth',2)
s1.FaceAlpha = 0.5;  
s1.LineWidth = 0.1; 
xlabel('$v_n$','FontName','Times New Roman','FontSize',10,'Interpreter','latex')
ylabel('$\dot{u}_n$','FontName','Times New Roman','FontSize',10,'Interpreter','latex')
zlabel('$\dot{v}_n$','FontName','Times New Roman','FontSize',10,'Interpreter','latex')
set(gca,'FontName','Times New Roman','FontSize',10);
set(gcf, 'Position', [600, 600, 300, 230]);
view([15 45])