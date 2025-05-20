function FnI = getFnI(Fn, Psi, Indp, I, order, Flag, k)
% Calculating the G^(p,q) and H^(p,q)
N = size(Psi{1}, 1);

FnI = zeros(N, 1);
Ind = Indp;

p = sum(I);
if (order == 2)&&(~isempty(Fn{order}))
    for p_1 = 1:p-1
        p_2 = p-p_1;
        mp1 = size(Ind{p_1},2);
        for k_1 = 1:mp1
            I1 = Ind{p_1}(:,k_1);
            I2 = I - I1;
            if min(I2) >= 0
                k_2 = find_index_VinA(I2,Ind{p_2});
                if isempty(k_2)
                    continue
                end
                Index_Fn = Fn{order}.I;
                Fn_vector = Fn{order}.vector;
                FnIkln = zeros(N, 1);
                for nm = 1:size(Index_Fn,2)
                    FnIkln = FnIkln + Fn_vector(:,nm)*Psi{p_1}(Index_Fn(1,nm),k_1)...
                        *Psi{p_2}(Index_Fn(2,nm),k_2);
                end
                if Flag == 0
                    FnI = FnI + FnIkln;
                else
                    FnI = FnI - Flag*FnIkln*(1-exp(Flag*1i*k*(I1(2)-I1(1))))*(1-exp(Flag*1i*k*(I2(2)-I2(1))));
                end
            else
                continue
            end
        end
    end

elseif (order == 3)&&(~isempty(Fn{order}))
 for p_1 = 1:p-2
     for p_2 = 1:p-p_1-1
        p_3 = p-p_1-p_2;
        mp1 = size(Ind{p_1},2);
        mp2 = size(Ind{p_2},2);
        for k_1 = 1:mp1
            for k_2 = 1:mp2
                I1 = Ind{p_1}(:,k_1);
                I2 = Ind{p_2}(:,k_2);
                I3 = I - I1 - I2;
                if min(I3) >= 0
                    k_3 = find_index_VinA(I3,Ind{p_3});
                    if isempty(k_3)
                        continue
                    end
                    Index_Fn = Fn{order}.I;
                    Fn_vector = Fn{order}.vector;
                    FnIkln = zeros(N, 1);
                    for nm = 1:size(Index_Fn,2)
                        FnIkln = FnIkln + Fn_vector(:,nm)*Psi{p_1}(Index_Fn(1,nm),k_1)...
                            *Psi{p_2}(Index_Fn(2,nm),k_2)*Psi{p_3}(Index_Fn(3,nm),k_3);
                    end
                    if Flag == 0
                        FnI = FnI + FnIkln;
                    else
                        FnI = FnI + FnIkln*(1-exp(Flag*1i*k*(I1(2)-I1(1))))...
                            *(1-exp(Flag*1i*k*(I2(2)-I2(1))))*(1-exp(Flag*1i*k*(I3(2)-I3(1))));
                    end
                else
                    continue
                end
            end
        end
     end
 end
end
end
