function model = generate_initial_modelXY(data, nstates, options)
% Generate a (poor) initial model for hidden Markov model (HMM).
%
% model = generate_initial_model(data, nstates, options)
%
% ARGUMENTS
%   data (cell array of 1D arrays) - observed trajectories of some real-valued signal
%   nstates - number of states
%   options - options data structure
%
% RETURNS
%   model (structure) - initial model parameters
%
% NOTES
%   A (poor) guess at the initial model parameters is generated.
%   This guess must be subsequently refined by iterations of either the maximum-likeihood or Bayesian HMM algorithms.
%
% TODO
%   Implement initial guess via EM algorithm
%   Guess number of states by multiple-normal density fitting?

if (options.verbosity >= 1)
    disp('******************************************************************************');
    disp('Bayesian HMM initial model generation');
end

% Make sure that data is a cell array (even if it's just one trajectory)
if ~iscell(data)
    data = {data};
end

% Create initial model from Gaussian mixture model.
%model = guess_paramsXY(data, nstates, options);
%model = em_gaussian_mixtureXY(data, nstates, options);
model.nstates = nstates;

% Generate initial transition matrix.
% Compute expected transition counts using only fractional state assignemnts.
Nij = zeros(nstates, nstates); % Nij(i,j) is fractional expected transition counts from i to j
% Construct transition count matrix from state trajectories.
mu = zeros(nstates, 1);
sigma = zeros(nstates, 1);
Pi = zeros(nstates, 1);
% for i = 1:nstates
%     mu(i) = model.states{i}.mu;
%     sigma(i) = model.states{i}.sigma;
%     Pi(i) = model.Pi(i);
% end
tmp_options = options;
tmp_options.verbosity = 0;
counterT = 0;
observations = zeros(0,1);
all_p_ti = zeros(0,2);

for trajectory_index = 1:length(data)
    if counterT >= 2e5
        break
    end
    % Extract trajectory of observables.
    o_t = RMSfilt2d(data{trajectory_index}',11)';
    observations = [observations reshape(data{trajectory_index}',[],1)];
    T = length(o_t);
    counterT = counterT + T;
    disp('...')
    tmp_model = em_gaussian_mixture({o_t}, nstates, tmp_options);
    for i = 1:nstates
        mu(i) = tmp_model.states{i}.mu;
        sigma(i) = tmp_model.states{i}.sigma;
        Pi(i) = tmp_model.Pi(i);
    end
    % Compute fractional assignment of samples to states
    p_ti = computeStateProbabilities(o_t, mu, sigma, Pi);
    all_p_ti = [all_p_ti;p_ti;p_ti];
    % Accumulate transition counts from this trajectory.
    for t = 1:(T-1)
        Nij = Nij + p_ti(t,:)' * p_ti(t+1,:);
    end
end
% Update transition matrix estimate.
model.Tij = transition_matrix_mle(Nij, options);
clear mu sigma Pi Nij;

% Compute stationary probability.
model.Pi = mean(p_ti);
%model.Pi = stationary_probability(model.Tij);
model.logPi = log(model.Pi);

% Update Gaussian parameters.

for i = 1:nstates
    model.states{i}.mu = sum(all_p_ti(:,i) .* observations) / sum(all_p_ti(:,i));
    model.states{i}.sigma = sqrt(sum(all_p_ti(:,i) .* (observations - model.states{i}.mu).^2) / sum(all_p_ti(:,i)));
end

% Assign samples based on probability alone.
%disp('Generating initial guess of state trajectory...');
%mloptions = options;
%%mloptions.updateMethod = 'maximum-likelihood';
%model = update_state_trajectories(data, model, mloptions);

if (options.verbosity >= 1)
    disp('******************************************************************************');
end

return
