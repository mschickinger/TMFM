function [simStraj, simXY, simRMS, simT] = sim_twostate_XY_RMS(params)

tau = params.tau;
tpf = params.tpf*2/1000;
L = params.L;
simMu = params.mu;
simSigma = params.sigma;
wSize = params.wSize;

% start state
Pcrit = tau(1)/sum(tau);
tmp = random('Uniform',0,1e6);
if tmp/1e6 <= Pcrit
    S1 = 1;
else
    S1 = 2;
end

% simulate trajectories
simStraj = zeros(1,L);
simXY = zeros(2,L);
allT = zeros(ceil(2*tpf*L/sum(tau)),1);

sumT = 0;
st = S1;
counter = 0;
while sumT < (L-0.5)*tpf
    counter = counter+1;
    allT(counter) = exprnd(tau(st));
    tmpF = ceil(sumT/tpf+0.5);
    sumT = sumT + allT(counter);
    tmpE = min(L,round(sumT/tpf));
    if tmpE>=tmpF
        simStraj(tmpF:tmpE) = st;
        simXY(:,tmpF:tmpE) = random('norm',simMu,simSigma(st),2,tmpE-tmpF+1);
    end
    st = mod(st,2)+1;
end
allT = allT(1:counter-1);
simRMS = RMSfilt2d(simXY,wSize);

% 'real' lifetime distributions
for i = 2:-1:1
    simT{1,i} = allT((1+(S1~=i)):2:end);
end
simT{S1}(1) = [];
    
return