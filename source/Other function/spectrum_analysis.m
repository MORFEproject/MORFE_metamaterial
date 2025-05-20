function Model = spectrum_analysis(Model,n,Flag)
% Providing the eigenvalue and eigenvectors
switch Flag
    case 'numerical'
        % Solving the numerical eigenproblem
        M = Model.M;
        C = Model.C;
        K = Model.K;
        
        A = [C M;M zeros(size(M))];
        B = blkdiag(-K,M);
        
        [PhTol,Lamb] = eigs(B,A);
        
        omega = diag(Lamb);
        [~,so] = sort(imag(omega),'ascend');
        Omega = omega([so(n+1:end);flip(so(1:n))]);
        Ph = PhTol(:,[so(n+1:end);flip(so(1:n))]);
        if rank(Ph) >= 4
            Y = Ph/diag(sqrt(diag(Ph.'*A*Ph)));
            L = diag(diag(Y.'*B*Y));
        else
            Y = Ph;
            L = diag(Omega);
        end
    case 'symbolic'
        % Solving the symbolic eigenproblem
        zero_matrix = zeros([2,2]);
        is_all_zero = all(all(Model.Cn{1} == zero_matrix));
        if is_all_zero
            % For undamped case
            syms v1 v2 om
            V = [v1 v1*1i;
                v2 v2*1i];                                          % Describing the wave form
            Y = [V(:,1) V(:,2);
                om*1i*V(:,1) -om*1i*V(:,2)];          % Eigenvectors of concerned dispersion branch
            L = [om*1i,0;
                0,-om*1i];                                          % Eigenvalues of concerned dispersion branch
        else
            % For underdamped case
            syms v1 v2 v1c v2c lamb1 lamb1c
            V = [v1 v1c;
                v2 v2c];                                              % Describing the wave form
            Y = [V(:,1) V(:,2);
                lamb1*V(:,1) lamb1c*V(:,2)];           % Eigenvectors of concerned dispersion branch
            L = [lamb1,0;
                0,lamb1c];                                          % Eigenvalues of concerned dispersion branch
        end
end

spectrum = struct('L',L,'Y',Y,'Lambda',diag(L));
Model.spectrum = spectrum;
end