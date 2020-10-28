function p_ti = computeStateProbabilities( o_t, mu, sigma, Pi )
%    Calculation of emission probabilities for each state.
%
%    o_t    trajectory of observations
%    mu     mu[i] is mean of emission model of state i
%    sigma  sigma[i] is standard deviation of emission model of state i
%    Pi     Pi[i] is the normalized weight or equilibrium probability of state i
%    p_ti   p_it[t][i] is the probability that o_t[t] came from state i

T = length(o_t);
nstates = length(mu);
log_Pi = log(Pi);
log_p_ti(T,nstates) = 0;

for t = 1:T
    for i = 1:nstates
        log_p_ti(t,i) = log_Pi(i) + logEmissionProbability(o_t(t), mu(i), sigma(i));
    end
    log_denom = logsum(log_p_ti(t,:));
    log_p_ti(t,:) = log_p_ti(t,:) - log_denom;
end

p_ti(T,nstates) = 0;

for t = 1:T
    for i = 1:nstates
        p_ti(t,i) = exp(log_p_ti(t,i));
    end
    denom = sum(p_ti(t,:));
    p_ti(t,:) = p_ti(t,:)./denom;
end

    function log_P = logEmissionProbability(o, mu, sigma)
        log_P = - 1/2*log(2*pi) - log(sigma) - 1/2*((o-mu)/sigma)^2;
    end


end

