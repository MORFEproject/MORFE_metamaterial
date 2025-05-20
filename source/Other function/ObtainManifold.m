function manifold = ObtainManifold(manifold, outdof)
Out = cell(length(outdof),1);
for ni = 1:length(manifold)
    for nj = 1:length(outdof)
        Out{nj} = zeros(size(manifold{ni}.Orbits{1},1),length(manifold{ni}.Orbits));
        for nk = 1:length(manifold{ni}.Orbits)
            Out{nj}(:,nk) = manifold{1}.Orbits{nk}(:,outdof(nj));
        end
    end
    manifold{ni}.out = Out;
end
end