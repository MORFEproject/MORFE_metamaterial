function ReducedModelFun = getReducedModelFun(ReducedModel,varargin)
if (nargin > 1)&&strcmp(varargin(1),'real')
    f = ReducedModel.f;
    ReducedModelFun = @(t,z) RealReducedModelFunction(t,z,f);
else
    I = ReducedModel.I;
    f = ReducedModel.f;
    ReducedModelFun = @(t, z) ReducedModelFunction(t, z, f, I);
end
end

function dz = ReducedModelFunction(t, z, f, I)
fTotal = [];
ITotal = [];
for ni = 1:length(f)
    fTotal = [fTotal f{ni}];
    ITotal = [ITotal I{ni}];
end
dz = fTotal*prod(z.^ITotal,1).';
end

function dz = RealReducedModelFunction(t, z, f)
fTotal = [];
ITotal = [];
for ni = 1:length(f)
    fTotal = [fTotal f{ni}.vector];
    ITotal = [ITotal f{ni}.I];
end
dz = fTotal*prod(z.^ITotal,1).';
end