function [Orbits, Tn] = ObtainWaveAttenuation_MMS(fun_pm, U0, V0, Tspan)
xn0 = U0(end-1:end); vn0 = V0(end-1:end);
f_residual = @(var) norm([fun_pm.Psi(var(1),var(2))-xn0;
    fun_pm.Upsilon(fun_pm.om(var(1)),var(1),var(2))-vn0]);
x0 = [1, 0]; % 初始点
opt_vars = fmincon(f_residual, x0,[],[],[],[],[0,0],[1,2*pi]);
an0 = opt_vars(1);
phan0 = opt_vars(2);
da = @(t, a) -fun_pm.om(a)*fun_pm.nc(fun_pm.om(a),a)*a;
tspan = linspace(Tspan(1),Tspan(end),2^14);
[Tn, Xn] = ode45(@(t,a) da(t,a), tspan, an0);
a_fun = @(t) interp1(Tn, Xn, t);
Pha = @(s, pha0) arrayfun(@(x) integral(@(t) fun_pm.om(a_fun(t)), 0, x), s) + pha0;
phan = Pha(Tn,phan0);
Orbits = zeros(4, length(Xn));
for ni = 1:length(Xn)
    yn = fun_pm.Psi(Xn(ni),phan(ni));
    dyn = fun_pm.Upsilon(fun_pm.om(Xn(ni)), Xn(ni), phan(ni));
    Orbits(:, ni) = [yn;dyn];
end
end