function ReducedModel = ObtainDispersion_Analytic_CNF_BP(parameters)
% Truncation order limited to 3
p = 3;

% Build dynamic model for symbolic DPIM calculation
model_DPIM = Build_cell_dpim(parameters,parameters.k);

% Providing the linear dispersion spectrum
model_DPIM = spectrum_analysis(model_DPIM,[],'symbolic');

% Providing the multi index for n DOF variables up to order p
Indp = getIndp(2, p);
model_DPIM.Indp = Indp;

% Performing the DPIM to calculate the CNF-BP
tic
ReducedModel = InvariantManifold_Analytic_CNF_BP(model_DPIM);
toc
end
