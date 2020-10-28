function [log_alpha_ti, log_beta_ti] = fwd_bwd_XYsegments(o_t, segments, models, log_rho_i)
% Forward-backward algorithm to compute partial trajectory probabilities.
% 
% ARGUMENTS
%   o_t (T) [real array] - trajectory of real-valued observations
%   model [structure] - model structure
%   log_rho_i (N) [real array] - log probability of initial state
%
% RETURNS
%   log_alpha_ti (T,N) [real array] - logarithms of forward probabilities  - alpha_ti(t,i) = P(O_1 = o_1, ..., O_t = o_t, X_t = i | model)
%   log_beta_ti (T,N) [real array] - logarithms of backward probabilities - beta_ti(t,i)  = P(O_{t+1} = o_{t+1}, ..., O_T = o_t | X_t = i, model)
%
% TODO
%   * This could be speeded up somewhat if we could compute a vector of log-emission probabilities at once
%   * Can these be written in vectorized form?

% Constants.
T = length(o_t); % number of observations in this trajectory
N = models{1}.nstates; % number of states
K = size(segments,1);

% Make sure model log probability and transition matrix is current.
models{1}.logTij = log(models{1}.Tij);
models{1}.logPi = log(models{1}.Pi);

% Compute all log_emission_probabilities up front:
log_emission_probabilities = zeros(T,N);
for k = 1:K
    tmpS = segments(k,:);
    for j = 1:N
        log_emission_probabilities(tmpS(1):tmpS(2),j) = models{k}.states{j}.log_emission_probability(o_t((tmpS(1):tmpS(2)),:));
    end
end

%
% Forward part:
% Compute alpha_ti = P(O_1 = o_1, ..., O_t = o_t, X_t = i | model)
%

log_alpha_ti = zeros(T,N);

% Compute probabilities for t = 1.
for i = 1:N
    log_alpha_ti(1,i) = log_rho_i(i) + log_emission_probabilities(1,i);  
end

% Use recurrence relation to compute log probabilities for times 2..T.
for t = 2:T
    for j = 1:N
        log_alpha_ti(t,j) = logsum(log_alpha_ti(t-1,:) + models{1}.logTij(:,j)') + log_emission_probabilities(t,j);    
    end
end

%
% Backward part:
% Compute beta_it = P(O_{t+1} = o_{t+1}, ..., O_T = o_t | X_t = i, model)
%

log_beta_ti = zeros(T,N);

% Compute probabilities for t = T.
for i = 1:N
    log_beta_ti(T,i) = 0;
end

% Use recurrence relation to compute log probabilities for times (T-1)...1.
for t = (T-1):-1:1
    for i = 1:N
        log_beta_ti(t,i) = logsum(models{1}.logTij(i,:) ...
                                    + log_emission_probabilities(t+1,:) ...
                                    + log_beta_ti(t+1,:));
    end
end

return
