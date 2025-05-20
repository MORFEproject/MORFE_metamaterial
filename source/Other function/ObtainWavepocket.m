function [X0,V0] = ObtainWavepocket(X,pha,n,ndof,n_v)
% the X_HB is the wave mode, pha is the wave number, and n denote the
% number of elements with initial condition. 

%% Paraments set
H = ((size(X,1) - 3)/ndof - 1)/2;
I0 = 1:ndof; ID = ndof + (1:H*ndof);
IC = ndof+repmat(1:ndof,1,H)+ndof*kron(0:2:2*(H-1),ones(1,ndof));
IS = IC + ndof;
Q = zeros(ndof*(H+1),1);
Q(I0) = X(I0);
Q(ID) = X(IC)-1i*X(IS);

a = exp(log(10)*X(end));
Q = a*Q;
Om = X(end-2);

X0 = zeros(ndof,n);
V0 = X0;

zeta = -log(1-0.999)/n_v;
Ha_1 = @(n) (1-exp(-zeta*n)) .*(n <= n_v) + (n > n_v);
Ha_2 = @(n) (1-exp(-zeta*(length(n)+1-n))) .*(n >= length(n) - n_v + 1) + (n < length(n) - n_v + 1);

for ni = (1:n)
    H_iDFT = exp(-1i*pha*ni*(0:H)');
    X0(:,ni) = real(reshape(Q,2,[])*H_iDFT);
    V0(:,ni) = real(reshape(Q,2,[])*(H_iDFT.*(1i*Om*(0:H)')));
end
Ha = repmat(Ha_1(1:n),ndof,1).*repmat(Ha_2(1:n),ndof,1);

X0 = reshape(Ha.*X0,[],1);
V0 = reshape(Ha.*V0,[],1);
end
