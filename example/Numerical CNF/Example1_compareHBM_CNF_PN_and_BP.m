% Comparison of the effect of periodic conditions applied at physical coordinates or 
% normal coordinates on dispersion solutions (Case: quadratic nonlinearity)
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
beta1 = 0;                  % quadratic nonlinear stiffness in host oscillator
beta2 = 2;                  % quadratic nonlinear stiffness in attachment
gamma1 = 0;             % cubic nonlinear stiffness in host oscillator
gamma2 = 0;             % cubic nonlinear stiffness in attachment
c1 = 0;                       % linear damping in host oscillator
c2 = 0;                       % linear damping in attachment

% The parameter of model
parameters_model = struct('mh',mh,'ma',ma,'kappa1',kappa1,...
    'kappa2',kappa2,'beta1',beta1,'beta2',beta2,'gamma1',gamma1,'gamma2',gamma2,'c1',c1,'c2',c2);

%% Using HBM, DPIM with BP and PN strategies to solve the nonlinear dispersion solutions
% under different wavenumber
N = 20;
kn = linspace(pi/N,pi,N);
order = 9; %CNF's truncation order
H = 5; % Order of harmonic truncation

% The dispersion solution of acoustic branch
log10a_s = -4.5; log10a_e = 0; mode = 1;
[omega1{1}, ~, umax1{1}] = ObtainDispersion_HBM(parameters_model,H,mode,log10a_s,log10a_e, kn);
rho = 10.^(linspace(-4,-0.5,2e2));
[omega1{2}, ~, umax1{2}] = ObtainDispersion_CNF_BP(parameters_model, mode, rho, kn, order);
rho = 10.^(linspace(-4,-0.5,2e2));
[omega1{3}, ~, umax1{3}] = ObtainDispersion_CNF_PN(parameters_model, mode, rho, kn, order);

% The dispersion solution of optic branch
log10a_s = -4.5; log10a_e = 0; mode = 2;
[omega2{1}, ~, umax2{1}] = ObtainDispersion_HBM(parameters_model,H,mode,log10a_s,log10a_e, kn);
rho = 10.^(linspace(-4,-0.5,2e2));
[omega2{2}, ~, umax2{2}] = ObtainDispersion_CNF_BP(parameters_model, mode, rho, kn, order);
rho = 10.^(linspace(-4,-0.5,2e2));
[omega2{3}, ~, umax2{3}] = ObtainDispersion_CNF_PN(parameters_model, mode, rho, kn, order);

%% Plot 3D dispersion manifold in space of k-omega_n-log_10(u_n)
Umax = [5e-3, 1.5e-2, 1e-1];
figure 
hold on
box on
grid on
for ni = 1:N
    plot3(kn(ni)*ones(size(omega1{1}{ni})),omega1{1}{ni},log10(umax1{1}{ni}(:,1)),'-','Color',[161, 71, 67]/255,'LineWidth',1);
    plot3(kn(ni)*ones(size(omega2{1}{ni})),omega2{1}{ni},log10(umax2{1}{ni}(:,1)),'-k','LineWidth',1);

    Results1.om = omega1; Results1.um = umax1; 
    Results2.om = omega2; Results2.um = umax2;
    umax =Umax(1);
    Flag = 0;
    Dispersion1 = getDispersion(Results1, kn, umax, Flag);
    Dispersion2 = getDispersion(Results2, kn, umax, Flag);
    
    for nj = 1:length(kn)
        plot3(Dispersion1{1}{nj}.k, Dispersion1{1}{nj}.om, log10(umax)*ones(size(Dispersion1{1}{nj}.k)),'o','Color',...
            [161, 71, 67]/255,'MarkerFaceColor',[161, 71, 67]/255,'MarkerSize',2)
        plot3(Dispersion2{1}{nj}.k, Dispersion2{1}{nj}.om, log10(umax)*ones(size(Dispersion2{1}{nj}.k)),'ok','MarkerFaceColor','k','MarkerSize',2)
    end

    umax =Umax(2);
    Flag = 0;
    Dispersion1 = getDispersion(Results1, kn, umax, Flag);
    Dispersion2 = getDispersion(Results2, kn, umax, Flag);
    
    for nj = 1:length(kn)
        plot3(Dispersion1{1}{nj}.k, Dispersion1{1}{nj}.om, log10(umax)*ones(size(Dispersion1{1}{nj}.k)),'o','Color',...
            [161, 71, 67]/255,'MarkerFaceColor',[161, 71, 67]/255,'MarkerSize',2)
        plot3(Dispersion2{1}{nj}.k, Dispersion2{1}{nj}.om, log10(umax)*ones(size(Dispersion2{1}{nj}.k)),'ok','MarkerFaceColor','k','MarkerSize',2)
    end

    umax =Umax(3);
    Flag = 0;
    Dispersion1 = getDispersion(Results1, kn, umax, Flag);
    Dispersion2 = getDispersion(Results2, kn, umax, Flag);
    
    for nj = 1:length(kn)
        plot3(Dispersion1{1}{nj}.k, Dispersion1{1}{nj}.om, log10(umax)*ones(size(Dispersion1{1}{nj}.k)),'o','Color',...
            [161, 71, 67]/255,'MarkerFaceColor',[161, 71, 67]/255,'MarkerSize',2)
        plot3(Dispersion2{1}{nj}.k, Dispersion2{1}{nj}.om, log10(umax)*ones(size(Dispersion2{1}{nj}.k)),'ok','MarkerFaceColor','k','MarkerSize',2)
    end
end
y = [-5 5]; x = [-2*pi 2*pi];
[X, Y] = meshgrid(x, y);
Z1 = log10(5e-3)*ones(size(X));
Z2 = log10(1.5e-2)*ones(size(X));
Z3 = log10(10^(-1))*ones(size(X));
surf(X, Y, Z1, 'FaceColor', 'r', 'FaceAlpha', 0.2)
surf(X, Y, Z2, 'FaceColor', 'b', 'FaceAlpha', 0.2)
surf(X, Y, Z3, 'FaceColor', 'c', 'FaceAlpha', 0.2)
le = legend('ABs','OBs','FontName','Times New Roman','FontSize',10);
set(le,'box','off')
view([15 45])
zlim([-2.5, -0.5])
ylim([0 2.5])
xlim([0 pi])
zlabel('Amplitude of $\log_{10}{u_n}$','FontName','Times New Roman','FontSize',10,'Interpreter','latex')
xlabel('Wave number {\itk}','FontName','Times New Roman','FontSize',10)
ylabel('Nonlinear Frequency {\it\omega_{nl}}','FontName','Times New Roman','FontSize',10)
set(gca,'FontName','Times New Roman','FontSize',10);
set(gcf, 'Position', [600, 600, 300, 230]);

%% Plot dispersion spectrum under umax=5e-3
Results1.om = omega1; Results1.um = umax1; 
Results2.om = omega2; Results2.um = umax2;
umax = 5e-3;
Flag = 1;
Dispersion1 = getDispersion(Results1, kn, umax, Flag);
Dispersion2 = getDispersion(Results2, kn, umax, Flag);

figure
hold on
box on
grid on
for ni = 1:length(kn)
    plot(Dispersion1{1}{ni}.k, Dispersion1{1}{ni}.om,'s','Color',[161, 71, 67]/255,'LineWidth', 1)
    plot(Dispersion1{2}{ni}.k, Dispersion1{2}{ni}.om,'o','Color',[161, 71, 67]/255,'LineWidth', 1)
    plot(Dispersion1{3}{ni}.k, Dispersion1{3}{ni}.om,'x','Color',[161, 71, 67]/255,'LineWidth', 1)

    plot(Dispersion2{1}{ni}.k, Dispersion2{1}{ni}.om,'sk','LineWidth', 1)
    plot(Dispersion2{2}{ni}.k, Dispersion2{2}{ni}.om,'ok','LineWidth', 1)
    plot(Dispersion2{3}{ni}.k, Dispersion2{3}{ni}.om,'xk','LineWidth', 1)
end
xlim([0, pi])
ylim([0 2.5])
le = legend('EPMC-HBM (A)','CNF-BP (A)','CNF-PN (A)',...
    'EPMC-HBM (O)','CNF-BP (O)','CNF-PN (O)','FontName','Times New Roman','FontSize', 10);
set(le,'box','off')
xlabel('Wave number {\itk}','FontName','Times New Roman','FontSize',10)
ylabel('Nonlinear Frequency {\it\omega_{nl}}','FontName','Times New Roman','FontSize',10)
set(gca,'FontName','Times New Roman','FontSize',10);
set(gcf, 'Position', [600, 600, 300, 230]);

%% Plot dispersion spectrum under umax=1.5e-2
Results1.om = omega1; Results1.um = umax1; 
Results2.om = omega2; Results2.um = umax2;
umax = 1.5e-2;
Flag = 1;
Dispersion1 = getDispersion(Results1, kn, umax, Flag);
Dispersion2 = getDispersion(Results2, kn, umax, Flag);

figure
hold on
box on
grid on
for ni = 1:length(kn)
    plot(Dispersion1{1}{ni}.k, Dispersion1{1}{ni}.om,'s','Color',[161, 71, 67]/255,'LineWidth', 1)
    plot(Dispersion1{2}{ni}.k, Dispersion1{2}{ni}.om,'o','Color',[161, 71, 67]/255,'LineWidth', 1)
    plot(Dispersion1{3}{ni}.k, Dispersion1{3}{ni}.om,'x','Color',[161, 71, 67]/255,'LineWidth', 1)

    plot(Dispersion2{1}{ni}.k, Dispersion2{1}{ni}.om,'sk','LineWidth', 1)
    plot(Dispersion2{2}{ni}.k, Dispersion2{2}{ni}.om,'ok','LineWidth', 1)
    plot(Dispersion2{3}{ni}.k, Dispersion2{3}{ni}.om,'xk','LineWidth', 1)
end
xlim([0,pi])
ylim([0 2.5])
xlabel('Wave number {\itk}','FontName','Times New Roman','FontSize',10)
ylabel('Nonlinear Frequency {\it\omega_{nl}}','FontName','Times New Roman','FontSize',10)
set(gca,'FontName','Times New Roman','FontSize',10);
set(gcf, 'Position', [600, 600, 300, 230]);

%% Plot dispersion spectrum under umax=1e-1
Results1.om = omega1; Results1.um = umax1; 
Results2.om = omega2; Results2.um = umax2;
umax = 10^(-1);
Flag = 1;
Dispersion1 = getDispersion(Results1, kn, umax, Flag);
Dispersion2 = getDispersion(Results2, kn, umax, Flag);

figure
hold on
box on
grid on
for ni = 1:length(kn)
    % plot(Dispersion1{1}{ni}.k, Dispersion1{1}{ni}.om,'s','Color',[161, 71, 67]/255,'LineWidth', 1)
    % plot(Dispersion1{2}{ni}.k, Dispersion1{2}{ni}.om,'o','Color',[161, 71, 67]/255,'LineWidth', 1)
    % plot(Dispersion1{3}{ni}.k, Dispersion1{3}{ni}.om,'x','Color',[161, 71, 67]/255,'LineWidth', 1)

    plot(Dispersion2{1}{ni}.k, Dispersion2{1}{ni}.om,'sk','LineWidth', 1)
    plot(Dispersion2{2}{ni}.k, Dispersion2{2}{ni}.om,'ok','LineWidth', 1)
    plot(Dispersion2{3}{ni}.k, Dispersion2{3}{ni}.om,'xk','LineWidth', 1)
end
xlim([0,pi])
ylim([0 2.5])
xlabel('Wave number {\itk}','FontName','Times New Roman','FontSize',10)
ylabel('Nonlinear Frequency {\it\omega_{nl}}','FontName','Times New Roman','FontSize',10)
set(gca,'FontName','Times New Roman','FontSize',10);
set(gcf, 'Position', [600, 600, 300, 230]);