function [ khat, tauhat ] = get_corrected_rates( lts, Tmin, Tmax )
%
expcdf_mod = @(t,tau,Tmin,Tmax)(exp(-t./tau)-exp(-Tmin./tau))/(exp(-Tmax./tau)-exp(-Tmin./tau));
exppdf_mod = @(t,tau,Tmin,Tmax)1./tau.*exp(-t./tau)./(exp(-Tmin./tau)-exp(-Tmax./tau));

tauhat = [0 0];
pci = zeros(2);
for i = 1:2
    testpdf = @(t,tau)exppdf_mod(t,tau,Tmin(i),Tmax);
    testcdf = @(t,tau)expcdf_mod(t,tau,Tmin(i),Tmax);
    [tauhat(i), pci(i,:)] = mle(lts{i},'pdf',testpdf,'start',mean(lts{i}),'cdf',testcdf);
end

%
optimfun = @(k) ...
    (exp(-Tmin(2)*k(2))-1/(k(1)*tauhat(1))).^2 + ...
    (exp(-Tmin(1)*k(1))-1/(k(2)*tauhat(2))).^2;
              
kstart = 1./tauhat;
khat = fminsearch(optimfun,kstart);

end

