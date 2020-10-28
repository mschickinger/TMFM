function [Tij, mu, sigma] = get_models_params(models)
% Extract model parameters, returning them as a list of matrices and vectors.
%
% Adapted from: [Tij, mu, sigma] = extract_model_parameters(model)

Tij = models{1}.Tij;
nstates = models{1}.nstates;
nmodels = size(models,1);
mu = zeros(nstates,nmodels);
sigma = zeros(nstates,nmodels);
for i = 1:nstates
    for j = 1:nmodels
        mu(i,j) = models{j}.states{i}.mu;
        sigma(i,j) = models{j}.states{i}.sigma;
    end
end
mu = reshape(mu,[],1);
sigma = reshape(sigma,[],1);

return

  