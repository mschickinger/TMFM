function [mlmodels, state_trajectory] = mlhmmXYsegments(data, segments, nstates, options, varargin)
% Construct maximum-likelihod hidden Markov model (HMM), useful for initializing a Bayesian HMM.
%
% mlmodel = mlhmmXY(data, nstates, options, [initial_model])
%
% ARGUMENTS
%   data (2D array) - observed XY-trajectory from TPM data
%   segments (Kx2 array) - start and end frames of K segments
%   [segments(k+1,1) must equal segments(k,2)+1]
%   nstates (int) - number of states to fit
%   options - options structure to control model generation
%
% OPTIONAL ARGUMENTS
%   initial_model - model to initialized MLHMM procedure
%
% RETURNS
%   mlmodel (structure) - maximum likelihood model after fit to data
%
% NOTES
%   Trajectories provided in 'data' cell array may be of different lengths. 
%   The expectation-maximization (EM) algorithm is used for updates; this is a variation on Baum-Welch.
%   This may not produce the global maximum-likelihood model.
%
% TODO
%   Add ability for user to pass in optional 'options' structure.

% Parse input
p = inputParser;
addRequired(p,'data')
addRequired(p,'segments')
addRequired(p,'nstates')
addRequired(p,'options')
addOptional(p,'sigmaIn',[])
addOptional(p,'initial_model',[])

parse(p, data, segments, nstates, options, varargin{:})

data = p.Results.data;
segments = p.Results.segments;
nstates = p.Results.nstates;
options = p.Results.options;

% Import fast Java helper code.
import bhmm_helper

% DEBUG
RELATIVE_TOLERANCE = options.convergenceTolerance;

if isempty(p.Results.initial_model)
    % Generate a valid initial model (which may be a poor guess).
    model = generate_initial_modelXY(data, nstates, options);
elseif ~iscell(p.Results.initial_model)
    % Use specified initial model.
    if (options.verbosity >= 1)
        disp('Using supplied initial model');
    end
    model = p.Results.initial_model;
else
    % Use specified initial models cell array.
    if (options.verbosity >= 1)
        disp('Using supplied initial models cell array');
    end
    models = p.Results.initial_model;
end
if ~exist('models','var')
    if ~isempty(p.Results.sigmaIn)
        % assign external starting values for X/Y widths
        for i = 1:nstates
            model.states{i}.mu = 0;
            model.states{i}.sigma = p.Results.sigmaIn(i);
        end
    end
    
    % Make model emission parameters 2-dimensional
    for i = 1:nstates
        model.states{i}.mu = model.states{i}.mu*ones(1,2);
        model.states{i}.sigma = model.states{i}.sigma*ones(1,2);
    end

    % Create a models cell array - one field for each segment
    
    models = cell(size(segments,1),1);
    for n = 1:length(models)
        models{n} =  model;
    end
end

if (options.verbosity >= 1)
    disp('******************************************************************************');
    disp('Initial model for segmented HMM:');
    %log_likelihood
    show_models(models);
    disp('******************************************************************************');
end

if (options.verbosity >= 1)
    disp('Fitting maximum-likelihood HMM with EM algorithm');
end

% Run cycles of maximum-likelihood updating.
for iteration = 1:options.maximumIterations

    if (options.verbosity >= 2)
        disp('******************************************************************************');
        fprintf('EM iteration %d\n', iteration);
    end

    % Store old model.
    [old_Tij, old_muXY, old_sigmaXY] = get_models_paramsXY(models);

    % Update emission probability functions associated with each state.
    for n = 1:numel(models)
        for i = 1:nstates
            state = models{n}.states{i};
            % Define normal density for emission probability.
            % models{n}.states{i}.emission_probability = @(o) normpdf(o(:,2),state.mu(1),state.sigma(1)).*normpdf(o(:,2),state.mu(2),state.sigma(2));
            if any(state.sigma==0 | isnan(state.sigma))
                models{n}.states{i}.log_emission_probability = @(o) -Inf.*ones(length(o),1); % Prevents problems in intervals without any transitions
            else
                models{n}.states{i}.log_emission_probability = @(o) -log(2*pi) - log(prod(state.sigma)) - ...
                                                                (1/2)*(((o(:,1)-state.mu(1))./state.sigma(1)).^2+((o(:,2)-state.mu(2))./state.sigma(2)).^2);
            end
        end
    end

    % Compute expected transition and state-occupation quantities.  
    Nij = zeros(nstates, nstates);
    o_n = data';
    
    % Compute log probabilities of transitions (log_xi_tij) and observations (log_gamma_ti) with Baum-Welch.
    % perform baum-welch-algorithm for bivariate probability distribution
    % (xy) and in segments
    [log_xi_tij, log_gamma_ti] = baum_welchXYsegments(o_n, segments, models, options);    

    % Compute expected transition counts.
    Nij = Nij + squeeze(sum(exp(log_xi_tij),1));

    % Compute weights for observations for each state.
    w_ni = exp(log_gamma_ti); % w_ni(n,i) is probability observation n came from state i

    %o_n = reshape(o_n,[],1);

    % Determine maximum likelihood transition matrix using fractional transition counts.
    models{1}.Tij = transition_matrix_mle(Nij, options);
    models{1}.logTij = log(models{1}.Tij);
    models{1}.Pi = stationary_probability(models{1}.Tij);
    models{1}.logPi = log(models{1}.Pi);

    % Update state observation probabilities.
    for n = 1:numel(models)
        tmpS = segments(n,:);
        for i = 1:nstates
            state = models{n}.states{i};

            % Extract weights for this state.
            w_n = w_ni(tmpS(1):tmpS(2),i);

            % Determine weighted sample statistics.
            for coord = 1:2
                state.mu(coord) = sum(w_n .* o_n(tmpS(1):tmpS(2),coord)) / sum(w_n);
                state.sigma(coord) = sqrt(sum(w_n .* (o_n(tmpS(1):tmpS(2),coord) - state.mu(coord)).^2) / sum(w_n));
                if isnan(state.sigma(coord))
                    state.mu(coord) = 0;
                    state.sigma(coord) = 0;
                end
            end
            
            models{n}.states{i} = state;
        end
    end

    if (options.verbosity >= 2)
        show_models(models);
    end

    % Check convergence criteria by computing relative change.
    [Tij, muXY, sigmaXY] = get_models_paramsXY(models);  
    if (norm(Tij - old_Tij) / norm(Tij) < RELATIVE_TOLERANCE) && (norm(muXY - old_muXY) / norm(muXY) < RELATIVE_TOLERANCE) && (norm(sigmaXY - old_sigmaXY) / norm(sigmaXY) < RELATIVE_TOLERANCE)
        if (options.verbosity >= 1)
            fprintf('Relative convergence tolerance of %e achieved.\n', RELATIVE_TOLERANCE);
        end
        break;
    end
end
clear old_model;

% Update emission probability functions associated with each state.
for n = 1:numel(models)
    for i = 1:nstates
        state = models{n}.states{i};
        % Define normal density for emission probability.
        % models{n}.states{i}.emission_probability = @(o) normpdf(o(:,2),state.mu(1),state.sigma(1)).*normpdf(o(:,2),state.mu(2),state.sigma(2));
        if any(state.sigma==0 | isnan(state.sigma))
            models{n}.states{i}.log_emission_probability = @(o) -Inf.*ones(length(o),1); % Prevents problems in intervals without any transitions
        else
            models{n}.states{i}.log_emission_probability = @(o) -log(2*pi) - log(prod(state.sigma)) - ...
                                                            (1/2)*(((o(:,1)-state.mu(1))./state.sigma(1)).^2+((o(:,2)-state.mu(2))./state.sigma(2)).^2);
        end
    end
end

if (options.verbosity >= 1)
    disp('******************************************************************************');
    disp('Final converged results:');
    %log_likelihood
    show_models(models);
    disp('******************************************************************************');
end

% Determine maximum-likelihood state trajectories.
if isfield(options, 'assign_states')
    assign_states = options.assign_states;
else
    assign_states = 1;
end
if assign_states == 1
    if (options.verbosity >= 1)
        disp('******************************************************************************');
        disp('State assignment by Viterbi Algorithm...');
    end

    o_t = data'; % o_t(t) is the observation at time t \in 1...T.
    [state_trajectory, ~ ] = viterbiXYsegments(o_t, segments, models);

    if (options.verbosity >= 1)
        disp('Done.');
        disp('******************************************************************************');
    end
end
  
% Return maximum-likelihood model.
mlmodels = models;

return


