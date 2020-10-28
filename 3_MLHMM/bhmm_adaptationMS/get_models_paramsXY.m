function [Tij, muXY, sigmaXY] = get_models_paramsXY(models)
% Extract model parameters, returning them as a list of matrices.
%
% Adapted from: [Tij, mu, sigma] = extract_model_parameters(model)

Tij = models{1}.Tij;
nstates = models{1}.nstates;
nmodels = size(models,1);
muXY = cell(nmodels,1);
sigmaXY = cell(nmodels,1);
for n = 1:nmodels
    muXY{n} = zeros(nstates,2); % 2 because of XY
    sigmaXY{n} = zeros(nstates,2);
    for i = 1:nstates
        muXY{n}(i,:) = models{n}.states{i}.mu;
        sigmaXY{n}(i,:) = models{n}.states{i}.sigma;
    end
end
muXY = cat(1,muXY{:});
sigmaXY = cat(1,sigmaXY{:});

return

  