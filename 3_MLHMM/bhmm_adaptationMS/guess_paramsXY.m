function model = guess_paramsXY(data, nstates, options)
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
    disp('Guessing gaussian mixture model parameters');
end

% DEBUG
options.enforce_state_ordering = true; % enforce ordering of states at this stage

% Create model structure.
model = struct();

% Store number of states.
model.nstates = nstates;

% Determine number of trajectories.
if ~iscell(data)
    data = {data};
end
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

w = 10;
W = min(nobservations,10000*w);
S = zeros(size(observations,1),W-w); % MS 170107: for bivariate observable. (e.g. x/y)
S(:,1) = sum(observations(:,1:1+w).^2,2);
for i = 2:size(S,2)
    S(:,i) = S(:,i-1) - observations(:,i-1).^2 + observations(:,i+w).^2;
end
S = sqrt(S./(w-1));

% Initial model parameters.
mu = zeros(nstates,1); % mu(i) is the Gaussian mean for state i
sigma = [sqrt(mean(min(S').^2));sqrt(mean(max(S').^2))]; % sigma(i) is the Guassian standard deviation for state i
%Pi = zeros(nstates,1); % Pi(i) is the normalized weight for Gaussian state i
obsW = observations(:,1:W);
Pi(2) = (var(obsW(:))-sigma(1)^2)/(sigma(2)^2-sigma(1)^2);
Pi(1) = 1-Pi(2);

if (options.verbosity >= 1)
    disp('initial guess:');
    for i = 1:nstates
    fprintf('state %5d : Pi = %8.6f mu = %8.3f  sigma = %8.3f\n', i, Pi(i), mu(i), sigma(i));
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
