function p_ti = computeStateProbabilitiesXY( o_t, mu, sigma, Pi )
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

for i = 1:nstates
    log_p_ti(:,i) = log_Pi(i) + logEmissionProbability(o_t, mu(i), sigma(i));
end
for t = 1:T
    log_denom = logsum(log_p_ti(t,:));
    log_p_ti(t,:) = log_p_ti(t,:) - log_denom;
end


p_ti = exp(log_p_ti);
denom = sum(p_ti,2);
for i = 1:nstates
    p_ti(:,i) = p_ti(:,i)./denom;
end

    function log_P = logEmissionProbability(o, mu, sigma)
        log_P = - log(2*pi) - 2*log(sigma) - ...
                    (1/2)*(((o(1,:)-mu)./sigma).^2 + ((o(2,:)-mu)./sigma).^2);
    end


end

