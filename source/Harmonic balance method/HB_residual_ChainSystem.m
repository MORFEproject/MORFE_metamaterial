function [R,dR,Q] = HB_residual_ChainSystem(X,parameter,H,N,varargin)
% Regarding the construction of harmonic residual equation, refer to the work 
% "Wang, T., Touze, C., Li, H.Q., Ding, Q., Nonlinear dispersion relationships 
% and dissipative properties of damped metamaterials embedding bistable attachments, 
% Nonlinear Dynamics, 2024"

%% Set up the parameters
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

%% generating the harmonic coefficients
I0 = 1:2; ID = 2+(1:H*2);
IC = 2+repmat(1:2,1,H) + 2*kron(0:2:2*(H-1),ones(1,2));
IS = IC + 2;
dX = eye(length(X));
Q = zeros(2*(H+1),1); dQ = zeros(size(Q,1),size(dX,2));
Q(I0) = X(I0); dQ(I0,:) = dX(I0,:);
Q(ID) = X(IC)-1i*X(IS); dQ(ID,:) = dX(IC,:)-1i*dX(IS,:);

%% provide some setting
inorm = varargin{1};
a = exp(log(10)*X(end));
da = log(10)*exp(log(10)*X(end))*dX(end,:);

Psi = Q; dPsi = dQ;
Q = Psi*a;
dQ = dPsi*a + Psi*da;

Om = X(end-2);
dOm = dX(end-2,:);

del = X(end-1);
ddel = dX(end-1,:);

Alpha = 2*Om*del;
dAlpha = 2*dOm*del + 2*Om*ddel;

if length(varargin)<2 || isempty(varargin{2})
    fscl = 1;
else
    fscl = varargin{2};
end

Indu = 1:2:2*(H+1);
Indv = Indu+1;

tau = (0:2*pi/N:2*pi-2*pi/N)';
H_iDFT = exp(1i*tau*(0:H));
H_iDFT_1 = exp(1i*(tau+k)*(0:H));
H_iDFT1 = exp(1i*(tau-k)*(0:H));
vn = real(H_iDFT*Q(Indv));
dvn = real(H_iDFT*dQ(Indv,:));
vndot = real(H_iDFT*(1i*Om*(0:H)'.*Q(Indv)));
dvndot = real(H_iDFT*(1i*Om*repmat((0:H)',1,size(dQ,2)).*...
    dQ(Indv,:))) + real(H_iDFT*(1i*(0:H)'.*Q(Indv))*dOm);
vnddot = real(H_iDFT*(-Om^2*((0:H).^2)'.*Q(Indv)));
dvnddot = real(H_iDFT*(-Om^2*repmat(((0:H).^2)',1,size(dQ,2)).*...
    dQ(Indv,:)) + H_iDFT*(-2*Om*((0:H).^2)'.*Q(Indv))*dOm);

un = real(H_iDFT*Q(Indu));
dun = real(H_iDFT*dQ(Indu,:));
undot = real(H_iDFT*(1i*Om*(0:H)'.*Q(Indu)));
dundot = real(H_iDFT*(1i*Om*repmat((0:H)',1,size(dQ,2)).*...
    dQ(Indu,:)) + H_iDFT*(1i*(0:H)'.*Q(Indu))*dOm);
unddot = real(H_iDFT*(-Om^2*((0:H).^2)'.*Q(Indu)));
dunddot = real(H_iDFT*(-Om^2*repmat(((0:H).^2)',1,size(dQ,2)).*...
    dQ(Indu,:)) + H_iDFT*(-2*Om*((0:H).^2)'.*Q(Indu))*dOm);

un_1 = real(H_iDFT_1*Q(Indu));
dun_1 = real(H_iDFT_1*dQ(Indu,:));
undot_1 = real(H_iDFT_1*(1i*Om*(0:H)'.*Q(Indu)));
dundot_1 = real(H_iDFT_1*(1i*Om*repmat((0:H)',1,size(dQ,2)).*...
    dQ(Indu,:)) + H_iDFT_1*(1i*(0:H)'.*Q(Indu))*dOm);

un1 = real(H_iDFT1*Q(Indu));
dun1 = real(H_iDFT1*dQ(Indu,:));
undot1 = real(H_iDFT1*(1i*Om*(0:H)'.*Q(Indu)));
dundot1 = real(H_iDFT1*(1i*Om*repmat((0:H)',1,size(dQ,2)).*...
    dQ(Indu,:)) + H_iDFT1*(1i*(0:H)'.*Q(Indu))*dOm);

Eq1n = (mh+ma)*unddot + ma*vnddot + (2*c1*undot-c1*undot_1 ...
    -c1*undot1) - Alpha*((mh+ma)*undot + ma*vndot) + (2*kappa1*un ...
    -kappa1*un_1-kappa1*un1) + beta1*(un-un_1).^2-beta1*(un-un1).^2 ...
    + gamma1*(un-un_1).^3 + gamma1*(un-un1).^3;
Eq2n = ma*(unddot+vnddot) +c2*vndot-Alpha*(ma*undot+ma*vndot)...
    + kappa2*vn + beta2*vn.^2 + gamma2*vn.^3;

dEq1n = (mh+ma)*dunddot + ma*dvnddot + (2*c1*dundot-c1*dundot_1 ...
    -c1*dundot1) - Alpha*((mh+ma)*dundot + ma*dvndot) + (2*kappa1*dun ...
    -kappa1*dun_1-kappa1*dun1) - ((mh+ma)*undot + ma*vndot)*dAlpha ...
    + 2*beta1*repmat((un - un_1),1,size(dun,2)).*(dun - dun_1) ...
    - 2*beta1*repmat((un - un1),1,size(dun,2)).*(dun - dun1) ...
    + 3*gamma1*repmat((un - un_1).^2,1,size(dun,2)).*(dun - dun_1) ...
    + 3*gamma1*repmat((un-un1).^2,1,size(dun,2)).*(dun-dun1);
dEq2n = ma*(dunddot + dvnddot) + c2*dvndot - Alpha*(ma*dundot + ma*dvndot)...
    +kappa2*dvn + 2*beta2*repmat(vn,1,size(dvn,2)).*dvn...
    +3*gamma2*repmat(vn.^2,1,size(dvn,2)).*dvn - (ma*undot + ma*vndot)*dAlpha;

%% Using the fft to calculate the residual equation
fnl = [Eq1n Eq2n];
Fnlc = fft(fnl(end-N+1:end,:))/N;
Fnl = [real(Fnlc(1,:));2*Fnlc(2:H+1,:)];
Rc = reshape(Fnl.',[],1);
dFnlc1 = fft(dEq1n(end-N+1:end,:))/N;
dFnlc2 = fft(dEq2n(end-N+1:end,:))/N;
dFnl1 = [real(dFnlc1(1,:));2*dFnlc1(2:H+1,:)];
dFnl2 = [real(dFnlc2(1,:));2*dFnlc2(2:H+1,:)];
dRc = zeros(size(Rc,1),size(dFnl1,2));
dRc(1:2:end,:) = dFnl1;
dRc(2:2:end,:) = dFnl2;

%% Scale dynamic force equilibrium (useful for numerical reasons)
Rc = 1/fscl*(Rc);
dRc = 1/fscl*(dRc);

%% Conversion from complex-exponential to sine-cosine representation
R = zeros(size(X,1)-1,1); dR = zeros(size(X,1)-1,size(X,1));
R(I0) = real(Rc(I0)); dR(I0,:) = real(dRc(I0,:));
R(IC) = real(Rc(ID)); dR(IC,:) = real(dRc(ID,:));
R(IS) = -imag(Rc(ID)); dR(IS,:) = -imag(dRc(ID,:));

%% Scale dynamic force equilibrium by modal amplitude
dR(1:end-2,:) = dR(1:end-2,:)/a-R(1:end-2)/a^2*da;
R(1:end-2) = R(1:end-2)/a;

%% Amplitude normalization
Mm = cell(H+1,1);
for ni = 1:H+1
    Mm{ni} = [mh+ma, ma;ma, ma];
end

R(end-1) = real(Psi'*blkdiag(Mm{:})*Psi-1);
dR(end-1,:) = real(2*(Psi'*blkdiag(Mm{:}))*dPsi);

%% Phase normalization
R(end) = (1:H)*imag(Psi(inorm+(2:2:H*2)));
dR(end,:) = (1:H)*imag(dPsi(inorm+(2:2:H*2),:));

%% Added a constraint for the rigid body mode
R(1) = X(1);
dR(1,:) = dX(1,:);
end