function Model = build_NonlinearAnnularChain(n, parameter)
% Build dynamic model of bistable nonlinear chain
% The dynamic model of elements is 
% M*\ddotXi + C*\dotXi + K*Xi + Fbl(xmi-1,xmi,xmi+1) + Fbn(xmi-1,xmi,xmi+1)
% + Fn(Xi) = km1*A*cos(omega*t)-cm1*A*omega*sin(omega*t)
% where
%               Xi = [xmi; vni]; dm = [1;0]; dn = [0;1]; vni = xni - xmi
%               M = [mm+mn mn;     D = [0 0;      K = [0 0;
%                        mn mn];                     0 cn];         0 kn1];
%               Fn = dn*kn3*(dn'*Xi)^3
%               Fbl = dm*cm1*dm'*(2\dotXi - \dotXi-1 - \dotXi+1)
%                        + dm*km1*dm'*(2Xi - Xi-1 - Xi+1);
%               Fbn = dm*km3*(dm'*(Xi - Xi-1))^3       
%                         + dm*km3*(dm'*(Xi - Xi+1))^3;
%               A = [A0; 0], if i = 1;
%
% Output: 
%               Model of the chain system
%               Model of the element in bistable nonlinear chain
% Input:    
%               n, denotes the number of element in nonlinear chain
%
%               parameter, The parameter of element of nonlinear chain, it is a
%               structure, and has fields, mm, km1, cm1, km3, mn, kn1, kn3,
%               cn1

%% Set the parameter
n_element = n;

mh = parameter.mh; mn = parameter.ma;
km1 = parameter.kappa1; km2 = parameter.beta1;
km3 = parameter.gamma1; cm1 = parameter.c1;
kn1 = parameter.kappa2; kn2 = parameter.beta2;
kn3 = parameter.gamma2; cn1 = parameter.c2;

dm = [1;0]; dn = [0;1];
%% Element's linear matrix
Mi0 = zeros(2); Mi1 = [mh+mn mn; mn mn]; Mi2 = Mi0;
Mi = [Mi0 Mi1 Mi2];

Di0 = - cm1*(dm*dm'); Di1 = 2*cm1*(dm*dm') + [0 0;0 cn1]; Di2 = - cm1*(dm*dm');
Di = [Di0 Di1 Di2];

Ki0 = - km1*(dm*dm'); Ki1 = 2*km1*(dm*dm') + [0 0;0 kn1]; Ki2 = - km1*(dm*dm');
Ki = [Ki0 Ki1 Ki2];

%% Assemble the linear matrix
M = zeros(2*n_element,2*n_element);
D = M;
K = M;
M(1:2,[end-1:end,1:4]) = M(1:2,[end-1:end,1:4]) + Mi;
D(1:2,[end-1:end,1:4]) = D(1:2,[end-1:end,1:4]) + Di;
K(1:2,[end-1:end,1:4]) = K(1:2,[end-1:end,1:4]) + Ki;
for ni = 2:n_element-1
    M(2*(ni - 1)+1:1:2*ni, 2*(ni - 1)-1:1:2*ni+2) = M(2*(ni - 1)+1:1:2*ni, 2*(ni - 1)-1:1:2*ni+2) + Mi;
    D(2*(ni - 1)+1:1:2*ni, 2*(ni - 1)-1:1:2*ni+2) = D(2*(ni - 1)+1:1:2*ni, 2*(ni - 1)-1:1:2*ni+2) + Di;
    K(2*(ni - 1)+1:1:2*ni, 2*(ni - 1)-1:1:2*ni+2) = K(2*(ni - 1)+1:1:2*ni, 2*(ni - 1)-1:1:2*ni+2) + Ki;
end
M(end-1:end,[end-3:end,1:2]) = M(end-1:end,[end-3:end,1:2]) + Mi;
D(end-1:end,[end-3:end,1:2]) = D(end-1:end,[end-3:end,1:2]) + Di;
K(end-1:end,[end-3:end,1:2]) = K(end-1:end,[end-3:end,1:2]) + Ki;

ndof = size(M,1);

%% Assemble the nonlinear force
if km2 ~= 0
    n_nonlinear = n_element;
    nonlinear_elements_km2 = cell(n_nonlinear,1);
    for ni = 1:n_nonlinear-1
        nonlinear_elements_km2{ni}.islocal = 1;
        nonlinear_elements_km2{ni}.ishysteretic = 0;
        nonlinear_elements_km2{ni}.type = 'quadraticSpring';
        force_direction = zeros(ndof,1); force_direction(2*(ni-1)+1:2*(ni+1),1) = [-dm;dm];
        nonlinear_elements_km2{ni}.force_direction = force_direction;
        nonlinear_elements_km2{ni}.stiffness = km2;
    end
    nonlinear_elements_km2{end}.islocal = 1;
    nonlinear_elements_km2{end}.ishysteretic = 0;
    nonlinear_elements_km2{end}.type = 'quadraticSpring';
    force_direction = zeros(ndof,1); force_direction([end-1:end,1:2],1) = [-dm;dm];
    nonlinear_elements_km2{end}.force_direction = force_direction;
    nonlinear_elements_km2{end}.stiffness = km2;
else
    nonlinear_elements_km2 = {};
end

if km3 ~= 0
    n_nonlinear = n_element;
    nonlinear_elements_km3 = cell(n_nonlinear,1);
    for ni = 1:n_nonlinear-1
        nonlinear_elements_km3{ni}.islocal = 1;
        nonlinear_elements_km3{ni}.ishysteretic = 0;
        nonlinear_elements_km3{ni}.type = 'cubicSpring';
        force_direction = zeros(ndof,1); force_direction(2*(ni-1)+1:2*(ni+1),1) = [-dm;dm];
        nonlinear_elements_km3{ni}.force_direction = force_direction;
        nonlinear_elements_km3{ni}.stiffness = km3;
    end
    nonlinear_elements_km3{end}.islocal = 1;
    nonlinear_elements_km3{end}.ishysteretic = 0;
    nonlinear_elements_km3{end}.type = 'cubicSpring';
    force_direction = zeros(ndof,1); force_direction([end-1:end,1:2],1) = [-dm;dm];
    nonlinear_elements_km3{end}.force_direction = force_direction;
    nonlinear_elements_km3{end}.stiffness = km3;
else
    nonlinear_elements_km3 = {};
end

if kn2 ~= 0
    n_nonlinear = n_element;
    nonlinear_elements_kn2 = cell(n_nonlinear,1);
    for ni = 1:n_nonlinear
        nonlinear_elements_kn2{ni}.islocal = 1;
        nonlinear_elements_kn2{ni}.ishysteretic = 0;
        nonlinear_elements_kn2{ni}.type = 'quadraticSpring';
        force_direction = zeros(ndof,1); force_direction(2*(ni-1)+1:2*ni,1) = dn;
        nonlinear_elements_kn2{ni}.force_direction = force_direction;
        nonlinear_elements_kn2{ni}.stiffness = kn2;
    end
else
    nonlinear_elements_kn2 = {};
end

if kn3 ~= 0
    n_nonlinear = n_element;
    nonlinear_elements_kn3 = cell(n_nonlinear,1);
    for ni = 1:n_nonlinear
        nonlinear_elements_kn3{ni}.islocal = 1;
        nonlinear_elements_kn3{ni}.ishysteretic = 0;
        nonlinear_elements_kn3{ni}.type = 'cubicSpring';
        force_direction = zeros(ndof,1); force_direction(2*(ni-1)+1:2*ni,1) = dn;
        nonlinear_elements_kn3{ni}.force_direction = force_direction;
        nonlinear_elements_kn3{ni}.stiffness = kn3;
    end
else
    nonlinear_elements_kn3 = {};
end

nonlinear_elements = [nonlinear_elements_km2; nonlinear_elements_km3; ...
    nonlinear_elements_kn2; nonlinear_elements_kn3];
%% Get the force vertor
Fext = zeros(ndof,1); Fext(1) = 1;

%% Construct the dynamic model of finite nonlinear chain

Model = struct('ndof',ndof,'M',M,'D',D,'K',K,'Fext',Fext);
Model.nonlinear_elements = nonlinear_elements;
OdeFun = @(t,X) ConstructFun(Model,t,X);
Model.Fun = OdeFun;
end

function dX = ConstructFun(Model,t,X)
%% Set the parameter
ndof = Model.ndof;
M = Model.M;
D = Model.D;
K = Model.K;

nonlinear_elements = Model.nonlinear_elements;

%% Get the matrix of the state-space model
A = [eye(ndof) zeros(ndof);
    zeros(ndof) M];
B = [zeros(ndof) eye(ndof);
    -K -D];

Fn = zeros(2*ndof,1);
for i = 1:length(nonlinear_elements)
    dn = nonlinear_elements{i}.force_direction;
    kn = nonlinear_elements{i}.stiffness;
    switch nonlinear_elements{i}.type
        case 'quadraticSpring'
            Fn = Fn + [zeros(ndof,1);-kn*dn*(dn'*X(1:ndof))^2];
        case 'cubicSpring'
            Fn = Fn + [zeros(ndof,1);-kn*dn*(dn'*X(1:ndof))^3];
    end
end

dX = A\(B*X + Fn);
end
