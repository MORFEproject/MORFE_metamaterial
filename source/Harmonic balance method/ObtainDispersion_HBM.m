function [Om, N_c, Umax, X_HB, Manifold] = ObtainDispersion_HBM(parameters, H, mode, log10a_s, log10a_e, kn, varargin)
% Using the EPMC-HBM framework to calculate the nonlinear dispersion
% solutions

%% Set up the parameters
mh = parameters.mh;
ma = parameters.ma;
kappa1 = parameters.kappa1;
kappa2 = parameters.kappa2;

n_kn = length(kn);
Om = cell(n_kn,1);
N_c = cell(n_kn,1);
Umax = cell(n_kn,1);
X_HB = cell(n_kn, 1);
Manifold = cell(n_kn, 1);

%% Calculate the dispersion solution for different wavenumber kn(ni)
for ni = 1:n_kn
    %% Providing the linear dispersion analysis
    M = [mh+ma, ma;ma,ma];
    K = [2*kappa1*(1-cos(kn(ni))),0;0,kappa2];
    [V,D] = eig(K,M);
    omega = sqrt(diag(D));
    V = V./sqrt(diag(V'*M*V));

    %% Calculation setting
    parameters.k = kn(ni);
    N = 2^7; inorm = 1;
    fscl = 1e-2;
    ds = 0.01;
    X0 = zeros(2*(2*H+1)+2,1);
    imode = mode; % Concerned dispersion branch
    X0(3:4) = V(:,imode);
    X0(end-1) = omega(imode);
    
    %% Using the Nlvib toolbox to perform a numerical continuation
    Sopt = struct('Dscale',[1e0*ones(size(X0,1)-2,1);X0(end-1);1e0;1e0], ...
        'dynamicDscale',1,'stepmax',2e3,'ds',ds,'dsmin',ds/10,'dsmax',5*ds);
    X_HB{ni} = solve_and_continue(X0, @(X) HB_residual_ChainSystem(X,parameters,...
        H,N,inorm,fscl),log10a_s,log10a_e,ds,Sopt);
    if ~isempty(varargin)
        Npoint = varargin{1};
    else
        Npoint = 2^6;
    end
    Manifold{ni}.X_HB = X_HB{ni};
    Manifold{ni} = getOrbits(Manifold{ni},2,Npoint,'nma');
    Umax{ni} = zeros(size(X_HB{ni},2),1);
    for nj = 1:length(Umax{ni})
        Umax{ni}(nj,1) = max(abs(Manifold{ni}.Orbits{nj}(:,1)));
        Umax{ni}(nj,2) = max(abs(Manifold{ni}.Orbits{nj}(:,2)));
    end
    Om{ni} = X_HB{ni}(end-2,:);
    N_c{ni} = X_HB{ni}(end-1,:);
end
end

function Manifold = getOrbits(Manifold, ndof, Npoint, type)
if nargin < 4
    type = 'nma';
end

X_HB = Manifold.X_HB;
norbits = size(X_HB, 2);
Orbits = cell(1, norbits);
Tn = Orbits;

for ni = 1:norbits
    switch type
        case 'nma'
            tn = linspace(0, 2*pi/X_HB(end-2, ni), Npoint);
        case 'frf'
            tn = linspace(0, 2*pi/X_HB(end, ni), Npoint);
    end
    [Yn, dYn] = yfunction(X_HB, ni, ndof, tn, zeros(ni,1), type);
    Orbits{ni} = [Yn dYn];
    Tn{ni} = tn;
end
Manifold.Orbits = Orbits;
Manifold.Tn = Tn;
end

function [yn,dyn] = yfunction(X_HB,ind,n,t,Phi_HB,type)
if nargin <= 5
    type = 'nma';
end

switch type
    case 'nma'
        if nargin <= 4
            phi = 0;
        else
            phi = Phi_HB(ind);
        end
        
        w = X_HB(end-2,ind);
        a = 10.^X_HB(end,ind);
        Qh = a*IVec(X_HB(1:end-3,ind),n);
        
        H = (size(Qh,2) - 1)/2;
        
        Hamonic = zeros(2*H+1,size(t,2));
        Hamonic(1,:) = ones(size(t));
        for i = 1:H
            Hamonic(2*i,:) = cos(i*w*(t - phi/w));
            Hamonic(2*i+1,:) = sin(i*w*(t - phi/w));
        end
        yn = (Qh*Hamonic)';
        dyn = (w*Qh*Delta(H)*Hamonic)';
    case 'frf'
        w = X_HB(end,ind);
        Qh = IVec(X_HB(1:end-1,ind),n);

        H = (size(Qh,2) - 1)/2;

        Hamonic = zeros(2*H+1,size(t,2));
        Hamonic(1,:) = ones(size(t));
        for i = 1:H
            Hamonic(2*i,:) = cos(i*w*t);
            Hamonic(2*i+1,:) = sin(i*w*t);
        end
        yn = (Qh*Hamonic)';
        dyn = (w*Qh*Delta(H)*Hamonic)';
end
end