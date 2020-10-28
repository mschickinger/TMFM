function model = em_gaussian_mixtureXY(data, nstates, options)
% Fit a Gaussian mixture model to observed data using the EM algorithm.
%
% model = em_gaussian_mixture(data, nstates, options)
%
% ARGUMENTS
%   data (cell array of 2D arrays) - observed trajectories of some real-valued signal
%   nstates - number of Gaussian components to fit
%   options - options data structure
%
% RETURNS
%   model (structure) - initial model parameters
%
% NOTES
%
%
% REFERENCES
%   

if (options.verbosity >= 1)
  disp('******************************************************************************');
  disp('Gaussian mixture model fitting by expectation maximization (EM) algorithm');
end

% DEBUG
options.enforce_state_ordering = true; % enforce ordering of states at this stage
RELATIVE_TOLERANCE = 1e-3;%options.convergenceTolerance;

% Create model structure.
model = struct();

% Store number of states.
model.nstates = nstates;

% Determine number of trajectories.
nTrajectories = length(data);

% Concatenate all observations.
observations = [];
for trajectoryIndex = 1:nTrajectories
    % Extract trajectory.
    o_t = data{trajectoryIndex};
    % Make into a row vector / horizontal array.
    o_t = reshape(o_t, [], length(o_t));
    % Concatenate observations.
    observations = [observations o_t];  
end
nobservations = length(observations);


% Generate initial guess by dividing data into two normal distributions
% centered around zero but with different widths.

w = 100;
W = min(nobservations,10000*w);
S = zeros(W-w,size(observations,1)); % MS 170107: for bivariate observable. (e.g. x/y)
for i = 1:length(S)
    for coord = 1:size(observations,1)
        S(i,coord) = std(observations(coord,i:i+w));
    end
end

% Initial model parameters.
mu = zeros(nstates,1); % mu(i) is the Gaussian mean for state i
sigma = [sqrt(mean(min(S).^2));sqrt(mean(max(S).^2))]; % sigma(i) is the Guassian standard deviation for state i
%Pi = zeros(nstates,1); % Pi(i) is the normalized weight for Gaussian state i
obsW = observations(:,1:W);
Pi(2) = (var(obsW(:))-sigma(1)^2)/(sigma(2)^2-sigma(1)^2);
Pi(1) = 1-Pi(2);

if (options.verbosity >= 1)
  disp('initial guess:');
  for i = 1:nstates
    disp(sprintf('state %5d : Pi = %8.6f mu = %8.3f  sigma = %8.3f', i, Pi(i), mu(i), sigma(i)));
  end
  disp('******************************************************************************');
end

% DEBUG: Initialize for example to give segmentation model best chance of working.
%mu(1) = 3;
%mu(2) = 4.7;
%mu(3) = 5.6;
%sigma(1) = 1.0;
%sigma(2) = 0.3;
%sigma(3) = 0.2;

%
% Carry out EM iterations.
%

niterations = 1000; % TODO: Set this to some large maximum once convergence assessment is added.

for iteration = 1:niterations
  if (options.verbosity >= 2)
    disp(sprintf('EM iteration %d', iteration));
  end

  % Store old model.
  old_mu = mu;
  old_sigma = sigma;

  % Compute p(l | x_i, mu, sigma) for each sample.
  p_ti = computeStateProbabilitiesXY(observations, mu, sigma, Pi); %MS 170107: wrote a Matlab function for this computation.

  % Update Gaussian parameters.
    for i = 1:model.nstates
        Pi(i) = mean(p_ti(:,i));
        mu(i) = 0;
        sigma(i) = 0;
        for coord = 1:2
            mu(i) = mu(i) + sum(p_ti(:,i)' .* observations(coord,:));
            sigma(i) = sigma(i) + sum(p_ti(:,i)' .* (observations(coord,:) - mu(i)).^2);
        end
        mu(i) = mu(i) / sum(p_ti(:,i)) / 2;
        sigma(i) = sqrt(sigma(i) / sum(p_ti(:,i)) / sum(p_ti(:,i)) / 2);
    end

  % Reorder states if necessary.
  if options.enforce_state_ordering
    % Sort by state widths.
    [~ , indices] = sort(sigma);
    % Permute states.
    Pi = Pi(indices);
    mu = mu(indices);
    sigma = sigma(indices);
  end
    
  if (options.verbosity >= 2)
    for i = 1:nstates
      disp(sprintf('state %5d : Pi = %8.6f mu = %8.3f  sigma = %8.3f', i, Pi(i), mu(i), sigma(i)));
    end
    disp('******************************************************************************');
  end

  % Assess convergence.
    if (norm((mu - old_mu)./old_mu) < 1e-2) && (norm((sigma - old_sigma)./old_sigma) < 1e-3)
        if (options.verbosity >= 2)
            fprintf('state means and variances converged to specified relative tolerances %e  and  %e\n', 1e-2, 1e-3);
        end
        break;
    else
        if (options.verbosity >= 2)
            fprintf('relative change in state means: %8.6f  , in state variances %8.6f\n', norm((mu - old_mu)./old_mu), norm((sigma - old_sigma)./old_sigma));
        end
    end
end

if (options.verbosity >= 1)
  disp('Converged results:');
  for i = 1:nstates
    disp(sprintf('state %5d : Pi = %8.6f mu = %8.3f  sigma = %8.3f', i, Pi(i), mu(i), sigma(i)));
  end
  disp('******************************************************************************');
end


% Store in model.
model.states = cell(nstates,1);
model.Pi = Pi;
for i = 1:nstates
  % Compute parameters that define emission probabilities of this state.
  state = struct();
  state.mu = mu(i);
  state.sigma = sigma(i);
  % Store this state.
  model.states{i} = state;  
end

return
