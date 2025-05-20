function R = getRInd(Psi, f, Indp, I, Flag, k)
% Calculating the mu^(p,q) and nu^(p,q)
N = size(Psi{1},1);
n = size(f{1},1);

R = zeros(N, 1);
Ind = Indp;
p = sum(I);

for s = 1:n
    for p_1 = 2:p-1
        p_2 = p - p_1 + 1;
        mp1 = size(Ind{p_1},2);
        for k_1 = 1:mp1
            I1 = Ind{p_1}(:,k_1);
            I2 = I-I1; I2(s) = I2(s)+1;
            if min(I2) >= 0
                k_2 = find_index_VinA(I2, Ind{p_2});
                if isempty(k_2)
                    continue
                end
                Rskn = I2(s)*Psi{p_2}(:,k_2)*f{p_1}(s,k_1);
                if Flag == 0
                    R = R + Rskn;
                else
                    R = R + Rskn*exp(Flag*1i*k*(I2(2)-I2(1)));
                end
            else
                continue
            end
        end
    end
end

end
