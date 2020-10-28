function [log_xi_tij, log_gamma_ti] = baum_welchXYsegments(o_t, segments, models, options)
% Baum-Welch algorithm for computing expected transition and state observation counts.
% 
% ARGUMENTS
%   o_t (T) [real array] - 2D-trajectory of real-valued observations
%   model [structure] - model structure
%   options [structure] - options structure
%
% RETURNS
%   log_xi_tij (T,N,N) [real array] - log-probability of observing transition  - xi_tij(t,i,j) = P(X_t = i, X_{t+1} = j | O, model)
%   log_gamma_ti (T,N) [real array] - log-probability observing state i at time t - gamma_ti(t,i) = P(X_t = i | O, model)

% Constants.
T = length(o_t); % number of observations in this trajectory
N = models{1}.nstates; % number of states

if options.equilibrium
    % Experiment starts in equilibrium.
    log_rho_i = models{1}.logPi;
elseif isfield(options, 'log_rho_i')
    % Use specified nonequilibrium starting conditions.
    log_rho_i = options.log_rho_i;
else
    % Use uniform initial starting conditions, which should be equivalent to complete ignorance.
    log_rho_i = models{1}.logPi * 0.0;
    log_rho_i = log_rho_i - log(sum(exp(log_rho_i)));
end

% Don't use Java code.
% if options.use_java == 42
%   % Prepare data for fast Java helper code.
%   mu = zeros(model.nstates,1);
%   sigma = zeros(model.nstates,1);
%   for i = 1:model.nstates
%     mu(i) = model.states{i}.mu;
%     sigma(i) = model.states{i}.sigma;
%   end  
%   % Update model parameters.
%   model.logTij = log(model.Tij); % elementwise logarithm - can we use a more numerically stable approach later?
%   model.logPi = log(model.Pi);
% 
%   % Perform Baum-Welch.
%   % log-probability of observing transition i to j at time t : xi_tij(t,i,j) = P(X_t = i, X_{t+1} = j | O, model)
%   log_alpha_ti = bhmm_helper.forwardAlgorithm(o_t, mu, sigma, model.logPi, model.logTij, log_rho_i);  
%   log_beta_ti = bhmm_helper.backwardAlgorithm(o_t, mu, sigma, model.logPi, model.logTij, log_rho_i);  
%   log_xi_tij = bhmm_helper.baumWelch_xi(o_t, mu, sigma, model.logPi, model.logTij, log_alpha_ti, log_beta_ti);
%   log_gamma_ti = bhmm_helper.baumWelch_gamma(o_t, mu, sigma, model.logPi, model.logTij, log_alpha_ti, log_beta_ti);
%   return
% end

% Compute forward-backward parameters.
[log_alpha_ti, log_beta_ti] = fwd_bwd_XYsegments(o_t, segments, models, log_rho_i);

% Make sure model log transition matrix is current.
models{1}.logTij = log(models{1}.Tij);

% First, compute log-probability of observation sequence o_t given model (used subsequently as a normalizing constant)
log_O = logsum(log_alpha_ti(T,:));

% Compute desired log-probabilities by Baum-Welch.
% TODO This could potentially be speeded up by using vector notation
log_xi_tij = zeros(T-1,N,N); % log-probability of observing transition i to j at time t : xi_tij(t,i,j) = P(X_t = i, X_{t+1} = j | O, model)
log_gamma_ti = zeros(T,N); % log-probability observing state i at time t : gamma_ti(t,i) = P(X_t = i | O, model)

for i = 1:N
    for j = 1:N
        % calculate segmented log-emission-probability vector
        tmp_log_em_prob = zeros(T,1);
        for k = 1:size(segments,1)
            tmpS = segments(k,:);
            tmp_log_em_prob(tmpS(1):tmpS(2)) = models{k}.states{j}.log_emission_probability(o_t(tmpS(1):tmpS(2),:));
        end
        log_xi_tij(:,i,j) = log_alpha_ti(1:end-1,i) ...
                            + tmp_log_em_prob(2:end) ...
                            + models{1}.logTij(i,j) ...
                            + log_beta_ti(2:end,j) ...
                            - log_O;
    end
%     for t = 1:(T-1)
%         log_gamma_ti(t,i) = logsum(log_xi_tij(t,i,:));
%     end
%     log_gamma_ti(T,i) = log_alpha_ti(T,i) + log_beta_ti(T,i) - log_O;
    log_gamma_ti(:,i) = log_alpha_ti(:,i) + log_beta_ti(:,i) - log_O;
end


return
