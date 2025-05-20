function [omega, nc, umax, manifold, fun] = ObtainDispersion_MMS(parameters, mode, rho, kn, varargin)
n_kn = length(kn);
omega = cell(n_kn,1);
nc = omega;
umax = nc;

if ~isempty(varargin)
    Npoint = varargin{1};
else
    Npoint = 2^6;
end
manifold = cell(n_kn,1);
for ni = 1:n_kn
parameters.k = kn(ni);
[omega{ni}, nc{ni}, umax{ni}, ~, manifoldni, fun] = Omega_pm_chain(parameters,rho,mode, Npoint);
manifold{ni} = manifoldni;
end
end

function [om, n_c, Amax_h,  Amax_a, manifold, Function] = Omega_pm_chain(parameter, rho, Mode, Npoint)
mh = parameter.mh;
ma = parameter.ma;
kappa1 = parameter.kappa1;
kappa2 = parameter.kappa2;
beta1 = parameter.beta1;
beta2 = parameter.beta2;
gamma1 = parameter.gamma1;
gamma2 = parameter.gamma2;
c1 = parameter.c1;
c2 = parameter.c2;
k = parameter.k; % wavenumber

M = @(kx) [mh+ma ma;ma ma];
K = @(kx) [2*kappa1*(1-cos(kx)) 0;0 kappa2];
C = @(kx) [2*c1*(1-cos(kx)) 0;0 c2];

[V,D] = eig(K(k),M(k));
om = sqrt(diag(D));
V = real(V/sqrt(diag(diag(V'*M(k)*V))));
om_acoustic = om(1);v_acoustic = V(:,1);
om_optic = om(2);v_optic = V(:,2);

if Mode == 1
    omega0 = om_acoustic;
    A0 = v_acoustic;
else
    omega0 = om_optic;
    A0 = v_optic;
end

ac0 = [0;beta2*A0(2)^2/2];
ac2 = [beta1*A0(1)^2*((1-exp(1i*k))^2-(1-exp(-1i*k))^2)/4;beta2*A0(2)^2/4];
C0 = zeros(2,1);
Kn = K(0);
C0(2) = -Kn(2,2)\ac0(2);
C2 = -(K(2*k)-(2*omega0)^2*M(2*k))\ac2;
ac3 = [gamma1*A0(1)^3*((1-exp(1i*k))^3+(1-exp(-1i*k))^3)/8;
    gamma2*A0(2)^3/8] + [beta1*A0(1)*C2(1)*((1-exp(2i*k))*(1-exp(1i*k))...
    -(1-exp(-2i*k))*(1-exp(-1i*k))); beta2*A0(2)*C2(2)];
C3 = -(K(3*k)-(3*omega0)^2*M(3*k))\ac3;
omega2 = rho.^2*beta2*A0(2)^2*(C0(2)+real(C2(2)))/omega0 ...
    + 3*rho.^2*gamma2*A0(2)^4/8/omega0 ...
    + rho.^2*beta1*A0(1)^2*real(C2(1)*((1-exp(2i*k))*(1-exp(-1i*k)) ...
    - (1-exp(-2i*k))*(1-exp(1i*k))))/omega0 ...
    + 3*rho.^2*gamma1*A0(1)^4*(1-cos(k))^2/2/omega0;
om = omega0 + omega2;
n_c = A0'*C(k)*A0/2./om + rho.^3*beta1*A0(1)^2*imag(C2(1)*((1-exp(2i*k))*(1-exp(-1i*k)) ...
    - (1-exp(-2i*k))*(1-exp(1i*k))))/om.^2 ...
    + rho.^3*beta2*A0(2)^2*imag(C2(2))/om.^2;
phan = linspace(0,2*pi,Npoint);
uh = rho'*A0(1)*cos(phan) + (rho.^2)'*C0(1)*ones(size(phan)) ...
    +2*(rho.^2)'*real(C2(1))*cos(2*phan)-2*(rho.^2)'*imag(C2(1))*sin(2*phan) ...
    +2*(rho.^3)'*real(C3(1))*cos(3*phan)-2*(rho.^3)'*imag(C3(1))*sin(3*phan);
Amax_h = max(abs(uh),[],2);
ua = rho'*A0(2)*cos(phan) + (rho.^2)'*C0(2)*ones(size(phan)) ...
    +2*(rho.^2)'*real(C2(2))*cos(2*phan)-2*(rho.^2)'*imag(C2(2))*sin(2*phan)...
    +2*(rho.^3)'*real(C3(2))*cos(3*phan)-2*(rho.^3)'*imag(C3(2))*sin(3*phan);
Amax_a = max(abs(ua),[],2);
vh = -om'.*rho'*A0(1)*sin(phan) - 2*om'.*2.*(rho.^2)'*real(C2(1))*sin(2*phan) ...
    - 2*om'.*2.*(rho.^2)'*imag(C2(1))*cos(2*phan) - 3*om'.*2.*(rho.^3)'*real(C3(1))*sin(3*phan) ...
    - 3*om'.*2.*(rho.^3)'*imag(C3(1))*cos(3*phan);
va = -om'.*rho'*A0(2)*sin(phan) - 2*om'.*2.*(rho.^2)'*real(C2(2))*sin(2*phan) ...
    - 2*om'.*2.*(rho.^2)'*imag(C2(2))*cos(2*phan)- 3*om'.*2.*(rho.^3)'*real(C3(2))*sin(3*phan) ...
    - 3*om'.*2.*(rho.^3)'*imag(C3(2))*cos(2*phan);

Orbits = cell(length(rho),1);
Tn = Orbits;
for ni = 1:length(rho)
    Orbits{ni} = [uh(ni,:)',ua(ni,:)',vh(ni,:)',va(ni,:)'];
    Tn{ni} = phan;
end
manifold.Orbits = Orbits; manifold.Tn = Tn;

Psi = @(r, p) [r*A0(1)*cos(p)+(r^2)*C0(1)*ones(size(p))+2*r^2*real(C2(1))*cos(2*p)-2*r^2*imag(C2(1))*sin(2*p)+2*r^3*real(C3(1))*cos(3*p)-2*r^3*imag(C3(1))*sin(3*p);
    r*A0(2)*cos(p)+r^2*C0(2)*ones(size(p))+2*r^2*real(C2(2))*cos(2*p)-2*r^2*imag(C2(2))*sin(2*p)+2*r^3*real(C3(2))*cos(3*p)-2*r^3*imag(C3(2))*sin(3*p)];
Upsilon = @(o, r, p) [-o*r*A0(1)*sin(p) - 2*o*2*(r^2)*real(C2(1))*sin(2*p)- 2*o*2*(r^2)*imag(C2(1))*cos(2*p)-3*o*2*(r^3)*real(C3(1))*sin(3*p)- 3*o*2*(r^3)*imag(C3(1))*cos(3*p);
    -o*r*A0(2)*sin(p) - 2*o*2.*(r^2)*real(C2(2))*sin(2*p)- 2*o*2.*(r^2)*imag(C2(1))*cos(2*p)-3*o*2.*(r^3)*real(C3(2))*sin(3*p)- 3*o*2.*(r^3)*imag(C3(1))*cos(3*p)];
nc_fun = @(o, r) A0'*C(k)*A0/2./o + r^3*beta1*A0(1)^2*imag(C2(1)*((1-exp(2i*k))*(1-exp(-1i*k)) ...
    - (1-exp(-2i*k))*(1-exp(1i*k))))/o^2 ...
    + r^3*beta2*A0(2)^2*imag(C2(2))/o^2;
om_fun = @(r) omega0 + r.^2*beta2*A0(2)^2*(C0(2)+real(C2(2)))/omega0 ...
    + 3*r.^2*gamma2*A0(2)^4/8/omega0 ...
    + r.^2*beta1*A0(1)^2*real(C2(1)*((1-exp(2i*k))*(1-exp(-1i*k)) ...
    - (1-exp(-2i*k))*(1-exp(1i*k))))/omega0 ...
    + 3*r.^2*gamma1*A0(1)^4*(1-cos(k))^2/2/omega0;
Function = struct('Psi',Psi,'Upsilon',Upsilon,'nc',nc_fun,'om',om_fun);
end