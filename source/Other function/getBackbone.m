function BB = getBackbone(Model, outdof, Opts)
% Obtain the backbone
np = length(Model.f);
Rgamma = zeros(np,1);
Igamma = zeros(np,1);
for ni = 1:2:np
    Rgamma(ni) = real(Model.f{ni}(1,(ni-1)/2 + 1));
    Igamma(ni) = imag(Model.f{ni}(1,(ni-1)/2 + 1));
end
Omega_fun = @(rho) BBfrequency(rho, Igamma);
Damping_fun = @(rho) BBDamping(rho, Rgamma);

Psi = Model.Psi;
Upsilon = Model.Upsilon;
I = Model.I;

rho = Opts.rho;
theta = Opts.theta;

omega = Omega_fun(rho);
nc = Damping_fun(rho)./omega;
z = rho'*exp(theta*1i)/2;
zc = rho'*exp(-theta*1i)/2;

N = size(Psi{1},1);
xoutdof = cell(length(outdof),1);
xmax = zeros(length(rho),length(outdof));
for nd = 1:length(outdof)
xoutdof{nd} = zeros(size(zc));
for ni = 1:length(rho)
    for p_i = 1:length(Psi)
        for l_pi = 1:size(Psi{p_i},2)
            if outdof(nd) < N+1
                xoutdof{nd}(ni, :) = xoutdof{nd}(ni, :) + Psi{p_i}(outdof(nd),l_pi)*...
                    prod([z(ni,:);zc(ni,:)].^I{p_i}(:,l_pi),1);
            else
                xoutdof{nd}(ni, :) = xoutdof{nd}(ni, :) + Upsilon{p_i}(outdof(nd)-N,l_pi)*...
                    prod([z(ni,:);zc(ni,:)].^I{p_i}(:,l_pi),1);
            end
        end
    end
end
xmax(:,nd) = max(abs(xoutdof{nd}),[],2);
end
BB = struct('Rgamma',Rgamma,'Igamma',Igamma,'Omegafun',Omega_fun,...
    'Dampingfun',Damping_fun,'omega',omega,'nc',nc,'xmax',xmax);
BB.xoudof = xoutdof;
end

function omega = BBfrequency(rho, Igamma)
% Nonlinear frequency
np = length(Igamma);
omega = zeros(size(rho));
for ni = 1:2:np
    omega = omega + Igamma(ni)*rho.^(ni-1)/(2.^(ni-1));
end
end

function damping = BBDamping(rho, Rgamma)
% Nonlinear damping
np = length(Rgamma);
damping = zeros(size(rho));
for ni = 1:2:np
    damping = damping - Rgamma(ni)*rho.^(ni-1)/(2.^(ni-1));
end
end