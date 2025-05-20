function [Orbits, Tn] = ObtainWaveAttenuation_HBM(X_HB, om, nc, X_HB0, Tspan)
X_HB_fun = @(log10a) interp1(X_HB(end,:), X_HB(1:end-1,:)', log10a);
nc_fun = @(log10a) interp1(X_HB(end,:), nc, log10a);
om_fun = @(log10a) interp1(X_HB(end,:), om, log10a);
manifold_fun = @(a, phi, ndof) ConstructManifold_HB(a, phi, ndof, X_HB_fun(log10(a))');
da = @(t, a) -nc_fun(log10(a))*om_fun(log10(a))*a;

tspan = linspace(Tspan(1), Tspan(end), 2^14);
[Tn, Xn] = ode45(@(t,a) da(t,a), tspan, 10^(X_HB0(end)));
a_fun = @(t) interp1(Tn, Xn, t);
Pha = @(s, pha0) arrayfun(@(x) integral(@(t) om_fun(log10(a_fun(t))), 0, x), s) + pha0;
phan = Pha(Tn,0);
Orbits = zeros(4, length(Xn));
for ni = 1:length(Xn)
    [yn, dyn] = manifold_fun(Xn(ni), phan(ni), 2);
    Orbits(:, ni) = [yn.';dyn.'];
end
end