function [Orbits, Tn] = ObtainWaveAttenuation_CNF(backbone, U0, V0, Tspan)
xn0 = U0(end-1:end); vn0 = V0(end-1:end);
manifold_fun = backbone{1}.ManifoldFun;
ReducedModel_fun = backbone{1}.ReducedModelFun;
f_residual = @(var) norm([real(manifold_fun.Psi([var(1)*exp(1i*var(2))/2 var(1)*exp(-1i*var(2))/2]))-xn0;
    + real(manifold_fun.Upsilon([var(1)*exp(1i*var(2))/2 var(1)*exp(-1i*var(2))/2]))-vn0]);
x0 = [1, 0]; % 初始点
opt_vars = fmincon(f_residual, x0,[],[],[],[],[0,0],[1,2*pi]);
an0 = opt_vars(1);
phan0 = opt_vars(2);
z10 = an0*exp(phan0*1i)/2;
z20 = an0*exp(-phan0*1i)/2;
z0 = [z10;z20];
tspan = linspace(Tspan(1), Tspan(end), 2^14);
[Tn,Xn] = ode45(@(t,z) ReducedModel_fun(t,z), tspan, z0);
% da = @(t, a) -backbone{1}.Dampingfun(a)*a;
% tspan = linspace(Tspan(1), Tspan(end), 2^14);
% [Tn, Xn] = ode45(@(t,a) da(t,a), tspan, an0);
% a_fun = @(t) interp1(Tn, Xn, t);
% Pha = @(s, pha0) arrayfun(@(x) integral(@(t) backbone{1}.Omegafun(a_fun(t)), 0, x), s) + pha0;
% phan = Pha(Tn,phan0);
Orbits = zeros(4, length(Xn));
% for ni = 1:length(Xn)
%     yn = real(manifold_fun.Psi([Xn(ni)*exp(1i*phan(ni))/2,Xn(ni)*exp(-1i*phan(ni))/2]));
%     dyn = real(manifold_fun.Upsilon([Xn(ni)*exp(1i*phan(ni))/2,Xn(ni)*exp(-1i*phan(ni))/2]));
%     Orbits(:, ni) = [yn;dyn];
% end
for ni = 1:length(Xn)
    yn = real(manifold_fun.Psi([Xn(ni,1),Xn(ni,2)]));
    dyn = real(manifold_fun.Upsilon([Xn(ni,1),Xn(ni,2)]));
    Orbits(:, ni) = [yn;dyn];
end
end