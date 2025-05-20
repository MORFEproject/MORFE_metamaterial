function [omega, nc, umax, Backbone, manifold] = ObtainDispersion_CNF_PN(parameters, mode, rho, kn, order, varargin)
n_kn = length(kn);
omega = cell(n_kn,1);
nc = omega;
umax = omega;
Backbone = umax;
manifold = umax;

Nmax = 2;
p = order;

for ni = 1:n_kn
model_DPIM = Build_cell_dpim(parameters, kn(ni));
imod = mode;           % Concerned dispersion branch

% Providing the linear dispersion spectrum
model_DPIM = spectrum_analysis(model_DPIM,Nmax,'numerical');

% Select master modes and provide the multiindex for n DOF variables up to
% order p
Master = [imod imod+Nmax];
Indp = getIndp(length(Master), p);
model_DPIM.Indp = Indp;
model_DPIM.Master = Master;

% Performing the DPIM to calculate the CNF-PN
tic
Compact_ReducedModel = InvariantManifold_CNF_PN(model_DPIM, imod);
toc

% Calculating the backbone and damping ratio for dispersion solutions
if ~isempty(varargin)
    Npoint = varargin{1};
else
    Npoint = 2^6;
end
outdof = 1:2;
Opts = struct('rho',rho, 'theta', linspace(0,2*pi,Npoint));
Backbone{ni} = getBackbone(Compact_ReducedModel, outdof, Opts);
omega{ni} = Backbone{ni}.omega;
nc{ni} = Backbone{ni}.nc;
umax{ni} = Backbone{ni}.xmax;

ManifoldFun = getManifoldFun(Compact_ReducedModel);
ReducedModelFun = getReducedModelFun(Compact_ReducedModel);
Backbone{ni}.ManifoldFun = ManifoldFun;
Backbone{ni}.ReducedModelFun = ReducedModelFun;
Orbits = cell(length(rho),1);
Tn = cell(length(rho),1);
for nj = 1:length(rho)
z1 = rho(nj)*exp(1i*linspace(0,2*pi,Npoint))/2;
z2 = rho(nj)*exp(-1i*linspace(0,2*pi,Npoint))/2;
z = [z1;z2].';
Xn = ManifoldFun.Psi(z);
Vn = ManifoldFun.Upsilon(z);
Orbits{nj} = real([Xn, Vn]);
Tn{nj} = linspace(0,2*pi, Npoint);
end
manifoldni.Orbits = Orbits;
manifoldni.Tn = Tn;
end
manifold{ni} = manifoldni;
end
