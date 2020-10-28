expcdf = @(t,tau)1-exp(-t./tau);
expcdf_mod = @(t,tau,Tmin,Tmax)(exp(-t./tau)-exp(-Tmin./tau))/(exp(-Tmax./tau)-exp(-Tmin./tau));
exppdf = @(t,tau)1./tau.*exp(-t./tau);
exppdf_mod = @(t,tau,Tmin,Tmax)1./tau.*exp(-t./tau)./(exp(-Tmin./tau)-exp(-Tmax./tau));

%%
%tau0 = [200 300];
Tmean = 10*[mean(vertcat(simT{:,1})) mean(vertcat(simT{:,2}))];
tau0 = simParams.tau/(2*simParams.tpf/1000);
tauSYN = [0 0];
for s = 1:2
    tmp = vertcat(simT{:,s});
    tmp(tmp<0.1)=[];
    synpdf = @(t,tau)exppdf_mod(t,tau,0.1,max(tmp)+0.1);
    syncdf = @(t,tau)expcdf_mod(t,tau,0.1,max(tmp)+0.1);
    tauSYN(s) = mle(tmp,'pdf',synpdf,'start',mean(tmp),'cdf',syncdf);
end
tauSYN = 10*tauSYN;

if ~isempty(T_missed{1})
    t1 = 1:ceil(max(T_missed{1}));
else
    t1 = 1:20;
end
if ~isempty(T_missed{2})
    t2 = 1:ceil(max(T_missed{2}));
else
    t2 = 1:20;
end
[t1,t2] = meshgrid(t1,t2);
RES = zeros(size(t1));
FracPost = zeros(size(RES));
FracMissed = zeros(size(RES));
Nmissed = numel(T_missed{1}) + numel(T_missed{2});
TAU = cell(size(t1));
TAUHAT = cell(size(t1));
h = waitbar(0,'');
for i = 1:numel(t1)
    waitbar(i/numel(t1),h,['Iteration ' num2str(i) ' of ' num2str(numel(t1))]);
    theta = [t1(i),t2(i)];
    [tmpSt, Npost] = lts_strict_cutoff(state_trajectories,theta);
    %disp(Npost)
    taus = [0 0];
    %DenomPost = 0;
    for s = 1:2
        %taus(s) = mean(tmpSt{s});
        testpdf = @(t,tau)exppdf_mod(t,tau,theta(s),max(tmpSt{s})+1);
        testcdf = @(t,tau)expcdf_mod(t,tau,theta(s),max(tmpSt{s})+1);
        taus(s) = mle(tmpSt{s},'pdf',testpdf,'start',mean(tmpSt{s}),'cdf',testcdf);
        %FracPost(i) = FracPost(i) + sum(T_post{s}<theta(s));
        %DenomPost = DenomPost + sum(T_post{s}<theta(s)) + sum(T_missed{s}<theta(s));
        FracMissed(i) = FracMissed(i) + sum(T_missed{s}<theta(s));       
    end
    %FracPost(i) = FracPost(i)/DenomPost;
    FracPost(i) = sum(Npost)/(sum(Npost)+FracMissed(i));
    FracMissed(i) = FracMissed(i)/Nmissed;
    TAU{i} = taus;
    kstart = 1./taus;
    optimfun = @(k) ...
                (exp(-theta(2)*k(2))-1/(k(1)*taus(1))).^2 + ...
                (exp(-theta(1)*k(1))-1/(k(2)*taus(2))).^2;
    kHat = fminsearch(optimfun,kstart);
    TAUHAT{i} = 1./kHat;
    RES(i) = sqrt(sum((kHat.*tauSYN - [1 1]).^2)); %relative deviations in k
    %RES(i) = sqrt(sum((TAUHAT{i}-tau0).^2)); %absolute deviations in tau
end
close(h)

%%
FracMissedMin = 0.8;
FracPostMax = 0.5;
RESmax = 0.01;
INDICES = find(FracPost(:)<=FracPostMax & FracMissed(:)>=FracMissedMin & RES(:)<=RESmax);
[~,findex] = min(FracPost(INDICES));
findex = INDICES(findex);
display('Best result for auto-cutoff:')
fprintf('Times:\t%i\t%i\n',t1(findex), t2(findex))
fprintf('tauIN\t%.2f\t%.2f\n',simParams.tau./(2*simParams.tpf/1000))
fprintf('tauSYN\t%.2f\t%.2f\n',tauSYN)
fprintf('TAUmle\t%.2f\t%.2f\n',TAU{findex})
fprintf('TAUcor\t%.2f\t%.2f\n',TAUHAT{findex})
fprintf('RESidx\t%.4f\n',RES(findex))
fprintf('Missed\t%.4f\n',FracMissed(findex))
fprintf('ProPost\t%.4f\n',FracPost(findex))
fprintf('* * * * *\n\n')

tmp = 'No';
%tmp = questdlg('Plot?','Plot','No');
if strcmp(tmp,'Yes')
    figure('Units', 'normalized', 'Position', [0 0 1 1]) 
    scatter3(FracPost(:),FracMissed(:),RES(:),'.')
    hold on
    scatter3(FracPost(findex),FracMissed(findex),RES(findex),'ro')
    xlabel('FracPost','FontSize',12)
    ylabel('FracMissed','FontSize',12)
    zlabel('RES','FontSize',12)
    set(gca,'Zscale','log')
    %zlim([1e-3 sqrt(2)*1e-2])

    figure('Units', 'normalized', 'Position', [0 0 1 1]) 
    surf(t1,t2,RES)
    set(gca,'ZScale','log')
end

save opTmin.mat t1 t2 RES FracPost FracMissed tau0 TAU TAUHAT INDICES findex

%% values for manual cutoffs, e.g. T50
display('Values for manual cutoff:')
times_lookup = [13 7];%T50?;
lindex = find(t1==times_lookup(1) & t2==times_lookup(2));
fprintf('Times:\t%i\t%i\n',t1(lindex), t2(lindex))
fprintf('tauIN\t%.2f\t%.2f\n',simParams.tau./(2*simParams.tpf/1000))
fprintf('tauSYN\t%.2f\t%.2f\n',tauSYN)
fprintf('TAUmle\t%.2f\t%.2f\n',TAU{lindex})
fprintf('TAUcor\t%.2f\t%.2f\n',TAUHAT{lindex})
fprintf('RESidx\t%.4f\n',RES(lindex))
fprintf('Missed\t%.4f\n',FracMissed(lindex))
fprintf('ProPost\t%.4f\n',FracPost(lindex))
fprintf('* * * * *\n\n')