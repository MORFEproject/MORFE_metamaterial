function Dispersion = getDispersion(Results, pha, Amax, Flag)
% This is a founction to get the dispertion curve through calculate the
% cross-section of the 2D manifold constructed by dNNMs under different
% wave number k
if nargin < 4
    Flag = 0;
end
om = Results.om;
um = Results.um;
Dispersion = cell(length(om),1);
Ninter = 2^18;
for ni = 1:length(om)
    Dispersion{ni} =  cell(length(om{ni}),1);
    for nj = 1:length(om{ni})
        A = [om{ni}{nj};um{ni}{nj}.'];
        Ind = 1:length(om{ni}{nj});
        Ind_int = linspace(1,length(om{ni}{nj}), Ninter);
        Aint = interp1(Ind,A',Ind_int)';
        Aselected = findPoints(Aint,Amax,1e-3);
        if (Flag == 1) & (size(Aselected,2) > 1)
            Aselected(:,2:end) = [];
        end
        Dispersion{ni}{nj}.om = Aselected(1,:);
        Dispersion{ni}{nj}.k = pha(nj)*ones(size(Aselected(1,:)));
    end
end
end