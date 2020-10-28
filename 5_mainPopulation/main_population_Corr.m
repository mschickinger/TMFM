%% Startup
clear variables
close all
run('my_prefs.m')

%% Extract relevant data from hop struct
Nsamples = str2double(inputdlg('How many samples?', 'Nsamples', 1, {'1'}));
Nspecies = str2double(inputdlg('How many species?', 'Nspecies', 1, {'1'}));
isnost = strcmp(questdlg('Before or after non-stick correction?','Nost?','Before','After','Before'),'After');
if isnost
    filedef = '*dataPostHMM_nost.mat';
else
    filedef = '*dataPostHMM.mat';
end
filepaths = cell(Nsamples,1);
hops = cell(Nsamples,1);
if Nspecies > 1
    % Names of species
    nameprompt = cell(Nspecies,1);
    defnames = cell(size(nameprompt));
    indexprompt = cell(Nspecies,1);
    for i = 1:Nspecies
        nameprompt{i} = ['Name of species ' num2str(i) ':'];
        defnames{i} = ['species ' num2str(i)];
        indexprompt{i} = ['Movies for species #' num2str(i) ':'];
    end
    Names = inputdlg(nameprompt, 'Names of species', 1, defnames);
else
    Names = {''};
end
% Movie indices in species
movInds = cell(Nsamples,Nspecies);
for i = 1:Nsamples
    [filename, pathname] = uigetfile([filesep filedef],['Pick sample file number ' num2str(i)]);
    filepaths{i} = [pathname filename];
    if isnost
        load(filepaths{i},'outputPostHMM_nost');
        hops{i} = outputPostHMM_nost.hop;
    else
        load(filepaths{i},'outputPostHMM');
        hops{i} = outputPostHMM.hop;
    end
    if Nspecies>1
        tmp = inputdlg(indexprompt, 'Movie indices',1);
        for j = 1:Nspecies
            movInds{i,j} = str2num(tmp{j});
        end
    else
        movInds{i} = unique(hops{i}.indices(:,1));
    end
end
clear tmp

%%
Ntraj = zeros(Nsamples,Nspecies);
for i = 1:Nsamples
    for j = 1:Nspecies
        Ntraj(i,j) = sum(ismember(hops{i}.indices(:,1),movInds{i,j}));
    end
end

%% set parameters for missed-event correction
Fmin = [13 7];
Tmin = round(0.1*Fmin,1);
Tmax = 0;

%%
resultsALL = cell(Nspecies,1);
tpfALL = cell(Nspecies,1);
lifetimes = cell(Nspecies,1);
for j = 1:Nspecies
    resultsALL{j} = cell(sum(Ntraj(:,j)),1);
    tpfALL{j} = zeros(sum(Ntraj(:,j)),1);
    lifetimes{j}.MEAN = zeros(sum(Ntraj(:,j)),2);
    lifetimes{j}.SUM = zeros(size(lifetimes{j}.MEAN));
    lifetimes{j}.N = lifetimes{j}.SUM;
    lifetimes{j}.ALL = cell(size(lifetimes{j}.MEAN));
    lifetimes{j}.ScatDat = zeros(sum(Ntraj(:,j)),6);
    counter = 0;
    for i = 1:Nsamples
        tmp = cell(size(movInds{i,j}));
        for m = 1:length(tmp)
            tmp{m} = hops{i}.tpf(movInds{i,j}(m))*ones(length(hops{i}.results{movInds{i,j}(m)}),1);
        end
        tpfALL{j}(counter + (1:Ntraj(i,j))) = vertcat(tmp{:});
        lifetimes{j}.ScatDat(counter + (1:Ntraj(i,j)),5) = ones(Ntraj(i,j),1)*i;
        lifetimes{j}.ScatDat(counter + (1:Ntraj(i,j)),6) = hops{i}.indices(ismember(hops{i}.indices(:,1),movInds{i,j}),1);
        resultsALL{j}(counter + (1:Ntraj(i,j))) = vertcat(hops{i}.results{movInds{i,j}});        
        counter = sum(Ntraj(1:i,j));
    end
    % fill in lifetimes
    for i = 1:size(lifetimes{j}.ALL,1)
        if ~isempty(resultsALL{j}{i}.state_trajectory)
            lifetimes{j}.ALL(i,:) = lts_strict_cutoff({resultsALL{j}{i}.state_trajectory},Fmin,{resultsALL{j}{i}.ex_int});
        end
        for k = 1:2
            lifetimes{j}.ALL{i,k} = lifetimes{j}.ALL{i,k}.*2*tpfALL{j}(i)/1000;
            lifetimes{j}.SUM(i,k) = sum(lifetimes{j}.ALL{i,k});
            lifetimes{j}.N(i,k) = length(lifetimes{j}.ALL{i,k});
            if ~isempty(lifetimes{j}.ALL{i,k})
                Tmax = max(Tmax,max(lifetimes{j}.ALL{i,k}));
            end
        end
        if ~isempty(resultsALL{j}{i}.hi) && ~isempty(resultsALL{j}{i}.lo)
            lifetimes{j}.ScatDat(i,1:2) = [mean(resultsALL{j}{i}.lo(:,2)) mean(resultsALL{j}{i}.hi(:,2))]*2*tpfALL{j}(i)/1000;
            lifetimes{j}.ScatDat(i,3:4) = [std(resultsALL{j}{i}.lo(:,2))/sqrt(size(resultsALL{j}{i}.lo,1)) ...
                                    std(resultsALL{j}{i}.hi(:,2))/sqrt(size(resultsALL{j}{i}.hi,1))]*2*tpfALL{j}(i)/1000;
        end
    end
    lifetimes{j}.MEAN = lifetimes{j}.SUM./lifetimes{j}.N;
end
%%
clear resultsALL tpfALL
%close all

mode = 3;
naive = 0;

lim = [1e-6 1; 1e-5 1; 1e-4 1];
lim(:,2) = 1-lim(:,1);
Nmin = 1;

tauRange = cell(Nspecies,2);
nRange = cell(Nspecies,size(lim,1));
INitial = cell(Nspecies,size(lim,1));
fINal = INitial;
taus = zeros(Nspecies,2);
tauhats = zeros(Nspecies,2);
Ns = cell(size(lim,1),1);
stDevs = cell(size(lim,1),1);
SEMs = cell(size(lim,1),1);
for l = 1:size(lim,1)
    Ns{l} = zeros(Nspecies,2);
    stDevs{l} = zeros(Nspecies,2);
    SEMs{l} = zeros(Nspecies,2);
end
for j = 1:Nspecies
    for l = 1:size(lim,1)
        display('Determining initial INdices...')
        for s = 1:2
            tmp = sort(vertcat(lifetimes{j}.ALL{:,s}));
            %tauRange{j,s} = linspace(tmp(ceil(0.01*length(tmp))),tmp(floor(0.99*length(tmp))),100);
            tauRange{j,s} = linspace(max(1.45*Tmin(s),tmp(ceil(0.01*length(tmp)))),tmp(floor(0.99*length(tmp))),100);
        end
        iRange = cell(length(tauRange{j,1}),length(tauRange{j,1}));
        nRange{j,l} = zeros(length(tauRange{j,1}),length(tauRange{j,1}));
        for i1 = 1:length(tauRange{j,1})
            for i2 = 1:length(tauRange{j,2})
                tmpI = zeros(size(lifetimes{j}.SUM,1),1);
                tmpTaus = [tauRange{j,1}(i1) tauRange{j,2}(i2)];
                Nmult = exp(Tmin./tmpTaus);
                alpha = Nmult - [1 1];
                meanTmissed = tmpTaus - exp(-Tmin./tmpTaus).*(Tmin+tmpTaus);
                n = 0;
                for k = 1:size(lifetimes{j}.SUM,1)
                    tmp = [0 0];
                    switch mode
                        case 0 % uncorrected
                            tmpN = lifetimes{j}.N(k,:);
                            tmpS = lifetimes{j}.SUM(k,:);
                        case 1 % 'traditional'
                            tmpN = [Nmult(1)*lifetimes{j}.N(k,1)+(1-1/Nmult(2))*lifetimes{j}.N(k,2) , Nmult(2)*lifetimes{j}.N(k,2)+(1-1/Nmult(1))*lifetimes{j}.N(k,1)];
                            tmpS = lifetimes{j}.SUM(k,:);
                        case 2 % corrected 'other-state-N'
                            tmpN = Nmult.*lifetimes{j}.N(k,:) + fliplr((Nmult-1).*lifetimes{j}.N(k,:));
                            tmpS = lifetimes{j}.SUM(k,:);
                        case 3 % correct both Ns and SUMs
                            if naive
                                tmpN = Nmult.*lifetimes{j}.N(k,:) + fliplr((Nmult-1).*lifetimes{j}.N(k,:));
                            else
                                tmpN = [(lifetimes{j}.N(k,1) + alpha(2)*lifetimes{j}.N(k,2)) , ...
                                        (lifetimes{j}.N(k,2) + alpha(1)*lifetimes{j}.N(k,1)) ].*(1+alpha)/(1-prod(alpha));
                            end
                            tmpS = lifetimes{j}.SUM(k,:) + (tmpN-lifetimes{j}.N(k,:)).*meanTmissed - fliplr((tmpN-lifetimes{j}.N(k,:)).*meanTmissed);
                        case 4 % correct one N and SUM
                            if naive
                                tmpN = lifetimes{j}.N(k,:) + fliplr((Nmult-1).*lifetimes{j}.N(k,:));         
                            else
                                tmpN = [(lifetimes{j}.N(k,1) + alpha(2)*lifetimes{j}.N(k,2)) , ...
                                        (lifetimes{j}.N(k,2) + alpha(1)*lifetimes{j}.N(k,1)) ]./(1-prod(alpha));
                            end
                            tmpS = lifetimes{j}.SUM(k,:) - fliplr((tmpN-lifetimes{j}.N(k,:)).*meanTmissed);
                    end
                    for state = 1:2
                        tmp(state) = erlangcdf(tmpS(state),1/tmpTaus(state),tmpN(state));
                    end
                    if min(tmp)>=lim(l,1) && max(tmp)<=lim(l,2)
                        n = n+1;
                        tmpI(k) = 1;
                    end
                end
                iRange{i1,i2} = find(tmpI==1);
                nRange{j,l}(i1,i2) = n;
            end
        end
        %
        tmpI = find(nRange{j,l}==max(nRange{j,l}(:)));
        INitial{j,l} = [];
        for i = reshape(tmpI,1,[])
            INitial{j,l} = union(INitial{j,l},iRange{i});
        end

        %
        
        display([num2str(length(INitial{j,l})) ' spots, ' num2str(length(vertcat(lifetimes{j}.ALL{INitial{j,l},1}))) ' bound lifetimes, ' ...
                num2str(length(vertcat(lifetimes{j}.ALL{INitial{j,l},2}))) ' unbound lifetimes.'])

        %

        display('Re-evaluation of taus and INdices...')
        fINal{j,l} = INitial{j,l};
        go_on = 1;
        while go_on
            INdicesOld = fINal{j,l};
            Tmax = max(vertcat(lifetimes{j}.ALL{:}));
            [khat, tauhats(j,:)] = get_corrected_rates({vertcat(lifetimes{j}.ALL{fINal{j,l},1}) vertcat(lifetimes{j}.ALL{fINal{j,l},2})},Tmin,Tmax);
            taus(j,:) = 1./khat;
            meanTmissed = taus(j,:) - exp(-Tmin./taus(j,:)).*(Tmin+taus(j,:));
            display(['Bound tau: ' num2str(taus(j,1)) ', unbound tau: ' num2str(taus(j,2))])
            Nmult = exp(Tmin.*khat);
            alpha = Nmult - [1 1];
            tmpI = zeros(size(lifetimes{j}.SUM,1),1);
            for k = 1:size(lifetimes{j}.SUM,1)
                tmp = [0 0];
                switch mode
                    case 0 % uncorrected
                        tmpN = lifetimes{j}.N(k,:);
                        tmpS = lifetimes{j}.SUM(k,:);
                    case 1 % 'traditional'
                        tmpN = [Nmult(1)*lifetimes{j}.N(k,1)+(1-1/Nmult(2))*lifetimes{j}.N(k,2) , Nmult(2)*lifetimes{j}.N(k,2)+(1-1/Nmult(1))*lifetimes{j}.N(k,1)];
                        tmpS = lifetimes{j}.SUM(k,:);
                    case 2 % corrected 'other-state-N'
                        tmpN = Nmult.*lifetimes{j}.N(k,:) + fliplr((Nmult-1).*lifetimes{j}.N(k,:));
                        tmpS = lifetimes{j}.SUM(k,:);
                    case 3 % correct both Ns and SUMs
                        if naive
                            tmpN = Nmult.*lifetimes{j}.N(k,:) + fliplr((Nmult-1).*lifetimes{j}.N(k,:));
                        else
                            tmpN = [(lifetimes{j}.N(k,1) + alpha(2)*lifetimes{j}.N(k,2)) , ...
                                    (lifetimes{j}.N(k,2) + alpha(1)*lifetimes{j}.N(k,1)) ].*(1+alpha)/(1-prod(alpha));
                        end
                        tmpS = lifetimes{j}.SUM(k,:) + (tmpN-lifetimes{j}.N(k,:)).*meanTmissed - fliplr((tmpN-lifetimes{j}.N(k,:)).*meanTmissed);
                    case 4 % correct one N and SUM
                        if naive
                            tmpN = lifetimes{j}.N(k,:) + fliplr((Nmult-1).*lifetimes{j}.N(k,:));         
                        else
                            tmpN = [(lifetimes{j}.N(k,1) + alpha(2)*lifetimes{j}.N(k,2)) , ...
                                    (lifetimes{j}.N(k,2) + alpha(1)*lifetimes{j}.N(k,1)) ]./(1-prod(alpha));
                        end
                        tmpS = lifetimes{j}.SUM(k,:) - fliplr((tmpN-lifetimes{j}.N(k,:)).*meanTmissed);
                end
                if any(lifetimes{j}.N(k,:)<Nmin)
                    tmp = [0 1];
                else
                    for state = 1:2
                        tmp(state) = erlangcdf(tmpS(state),khat(state),tmpN(state));
                    end
                end
                if min(tmp)>=lim(l,1) && max(tmp)<=lim(l,2)
                    tmpI(k) = 1;
                else
                    display(['Removing index ' num2str(k) ' from analysis.'])
                    display(['P_state1 = ' num2str(tmp(1),12) ', P_state2 = ' num2str(tmp(2),12)])
                end
            end
            fINal{j,l} = find(tmpI==1);
            for k = 1:2
                Ns{l}(j,k) = length(vertcat(lifetimes{j}.ALL{fINal{j,l},k}));
                stDevs{l}(j,k) = std(vertcat(lifetimes{j}.ALL{fINal{j,l},k}));
                SEMs{l}(j,k) = stDevs{l}(j,k)/sqrt(Ns{l}(j,k));
            end
            display([num2str(length(fINal{j,l})) ' spots, ' num2str(Ns{l}(j,1)) ' bound lifetimes, ' num2str(Ns{l}(j,2)) ' unbound lifetimes.'])
            go_on = ~isequal(fINal{j,l},INdicesOld);
        end
    end
end

% Correct MEAN lifetime values
display('Correcting mean values in individual particles...')
for j = 1:Nspecies
    Nmult = exp(Tmin./taus(j,:));
    alpha = Nmult - [1 1];
    meanTmissed = taus(j,:) - exp(-Tmin./taus(j,:)).*(Tmin+taus(j,:));
    lifetimes{j}.corrMEAN = zeros(size(lifetimes{j}.MEAN));
    for i = 1:size(lifetimes{j}.MEAN,1)
        %Ndiv = [Nmult(1)*lifetimes{j}.N(i,1)+(1-1/Nmult(2))*lifetimes{j}.N(i,2) ...
         %       Nmult(2)*lifetimes{j}.N(i,2)+(1-1/Nmult(1))*lifetimes{j}.N(i,1)];
        if naive
            tmpN = Nmult.*lifetimes{j}.N(i,:) + fliplr((Nmult-1).*lifetimes{j}.N(i,:));
        else
            tmpN = [(lifetimes{j}.N(i,1) + alpha(2)*lifetimes{j}.N(i,2)) , ...
                        (lifetimes{j}.N(i,2) + alpha(1)*lifetimes{j}.N(i,1)) ].*(1+alpha)./(1-prod(alpha));
        end
        tmpS = lifetimes{j}.SUM(i,:) - fliplr((tmpN-lifetimes{j}.N(i,:)).*meanTmissed);
        %tmpS = lifetimes{j}.SUM(i,:);
        lifetimes{j}.corrMEAN(i,:) = tmpS./tmpN;
    end
end
display('Done.')

% Surface and scatter plots with IN/OUT stats
%close all
surfig = figure('Units','normalized','Position',[0 0 1 1]);
for j = 1:Nspecies
    for l = 1:size(lim,1)
        %figure(surfig)
        subplot(size(lim,1),Nspecies,(l-1)*Nspecies+j)
        surf(tauRange{j,1},tauRange{j,2},nRange{j,l})
        title(['Species ' num2str(j) ' (' Names{j} '), Limit: ' num2str(lim(l,1)) '.'],'FontSize',12)
    end
end
dotfig = figure('Units','normalized','Position',[0 0 1 1]);
oCols = {'k','c','r'};
for j = 1:Nspecies
    %figure(dotfig)
    
    %INitial
    subplot(2,Nspecies,j)
    loglog(lifetimes{j}.corrMEAN(:,2),lifetimes{j}.corrMEAN(:,1),'.')
    hold on
    legendary = cell(size(lim,1),2);
    legendary = [[Names(j) Names(j)] ; legendary];
    for l = 1:size(lim,1)
        loglog(lifetimes{j}.corrMEAN(INitial{j,l},2),lifetimes{j}.corrMEAN(INitial{j,l},1),'o','Color',oCols{l})
        legendary{l+1,1} = [num2str(lim(l,1)) ' , ' num2str(length(INitial{j,l}))];
    end
    xlabel('mean unbound dwell times (s)','FontSize',12)
    ylabel('mean bound dwell times (s)','FontSize',12)
    title(['Species ' num2str(j) ' (' Names{j} '), Limit / Initial nIN: Legend'],'FontSize',14)
    legend(legendary(:,1),'FontSize',12)
    
    %fINal
    subplot(2,Nspecies,Nspecies+j)
    loglog(lifetimes{j}.corrMEAN(:,2),lifetimes{j}.corrMEAN(:,1),'.')
    hold on
    for l = 1:size(lim,1)
        loglog(lifetimes{j}.corrMEAN(fINal{j,l},2),lifetimes{j}.corrMEAN(fINal{j,l},1),'o','Color',oCols{l})
        legendary{l+1,2} = [num2str(lim(l,1)) ' , ' num2str(length(fINal{j,l}))];
    end
    xlabel('mean unbound dwell times (s)','FontSize',12)
    ylabel('mean bound dwell times (s)','FontSize',12)
    title(['Species ' num2str(j) ' (' Names{j} '), Limit / Final nIN: Legend'],'FontSize',14)
    legend(legendary(:,2),'FontSize',12)
end

% Choose INdices from candidates
L = 1; 
INdices = fINal(:,L);

%% Save data
if ~exist('pathname','var')
    pathname = cd;
end
cd(pathname)
display('Choose folder for results in save dialog')
[~,savepath] = uiputfile;
cd(savepath)
save lts_main_pop.mat filepaths hops Names Fmin Tmin Tmax lifetimes INdices taus tauhats Ns stDevs SEMs lim Nmin L Nspecies Nsamples
savefig(surfig,'surfig.fig')
savefig(dotfig,'dotfig.fig')
fileID = fopen('taus_main_pop.txt', 'w');
fprintf(fileID, 'tau_b\ttau_u\tN_b\tN_u\tstDev_b\tstDev_u\n');
for j = 1:Nspecies
    fprintf(fileID, [num2str(taus(j,1),'%.1f') '\t' num2str(taus(j,2),'%.1f') '\t']);
    fprintf(fileID, [num2str(Ns{L}(j,1),'%d') '\t' num2str(Ns{L}(j,2),'%d') '\t']);
    fprintf(fileID, [num2str(stDevs{L}(j,1),'%.1f') '\t' num2str(stDevs{L}(j,2),'%.1f') '\n']);
end
fclose(fileID);

%% Plots
close all
exppdf_mod = @(t,tau,Tmin,Tmax)1./tau.*exp(-t./tau)./(exp(-Tmin./tau)-exp(-Tmax./tau));
expcdf_mod = @(t,tau,Tmin,Tmax)(exp(-t./tau)-exp(-Tmin./tau))/(exp(-Tmax./tau)-exp(-Tmin./tau));
colors = {[204 0 0]/255,[0 102 153]/255};
for j = 1:Nspecies
    figure('Units', 'normalized', 'Position', [0 0 1 1], 'PaperPositionMode', 'auto');
    for state = 1:2
        lt_state = sort(vertcat(lifetimes{j}.ALL{INdices{j},state}));
        Tmax = max(lt_state)+0.1;
        %
        centers = logspace(-2,log10(Tmax),2e3);
        cumcts = zeros(size(centers));
        for i = 1:length(centers)
            cumcts(i) = sum(lt_state<=centers(i));
        end
        cumcts = cumcts/cumcts(end);
        tmax = centers(find(cumcts>.999,1));
        
        % histogram
        subplot(2,2,1+2*(state-1))
        hold off
        hg = histogram(lt_state,'Normalization','pdf','BinMethod','scott');
        edges = hg.BinEdges + Tmin(state)-hg.BinEdges(1);
        hg = histogram(lt_state,edges,'Normalization','pdf');
        hg.FaceColor = colors{state};
        hg.EdgeColor = 'white';
        hg.FaceAlpha = 1;
        hold on
        ts = linspace(Tmin(state),centers(end),1001);
        plot(ts,exppdf_mod(ts,tauhats(j,state),Tmin(state),Tmax), 'k--', 'LineWidth',1)
        plot(centers,exppdf(centers,taus(j,state)), 'k:', 'LineWidth',1)
        xlim([0 tmax])
        ax = gca;
        ax.TickDir = 'out';
        xlabel('lifetime (s)', 'FontSize', 14)
        ylabel('Relative frequency / Probability density', 'FontSize', 14)
        box off

        % cdf
        subplot(2,2,2+2*(state-1))
        hold off
        semilogx(centers,cumcts,'Color',colors{state},'LineWidth',1)
        hold on
        ts = logspace(log10(Tmin(state)-.1),log10(centers(end)),2e3);
        semilogx(ts,expcdf_mod(ts,tauhats(j,state),Tmin(state),Tmax),'k--','LineWidth',1)
        semilogx(ts,expcdf(ts,taus(j,state)),'k:','LineWidth',1)
        xlim([Tmin(state)-.1 tmax])
        ylim([0 1])
        ax = gca;
        ax.TickDir = 'out';
        ax.YTick = [0 .5 1];
        %ax.YTickLabel = {};
        xlabel('lifetime (s)', 'FontSize', 14)
        ylabel('Cumulative frequency / Probability', 'FontSize', 14)
        box off

        subplot(2,2,1)
        title(['PDF for state 1 (bound), \tau_{final} = ' num2str(round(taus(j,1),2)) ' s, \tau_{MLE} = ' num2str(round(tauhats(j,1),2)) ' s, N = ' num2str(Ns{L}(j,1))], 'FontSize', 18, 'Interpreter', 'tex')
        subplot(2,2,2)
        title(['CDF for state 1 (bound), \tau_{final} = ' num2str(round(taus(j,1),2)) ' s, \tau_{MLE} = ' num2str(round(tauhats(j,1),2)) ' s, N = ' num2str(Ns{L}(j,1))], 'FontSize', 18, 'Interpreter', 'tex')
        subplot(2,2,3)
        title(['PDF for state 2 (unbound), \tau_{final} = ' num2str(round(taus(j,2),2)) ' s, \tau_{MLE} = ' num2str(round(tauhats(j,2),2)) ' s, N = ' num2str(Ns{L}(j,2))], 'FontSize', 18, 'Interpreter', 'tex')
        subplot(2,2,4)
        title(['CDF for state 2 (unbound), \tau_{final} = ' num2str(round(taus(j,2),2)) ' s, \tau_{MLE} = ' num2str(round(tauhats(j,2),2)) ' s, N = ' num2str(Ns{L}(j,2))], 'FontSize', 18, 'Interpreter', 'tex')
        
        [~,h] = suplabel(['Species ' num2str(j) ': ' Names{j}],'t');
        set(h,'FontSize',20)
    end
    print('-dpng','-r150',['LifetimeDistributions_' Names{j} '.png'])
end

%% Produce text files for main population scatter plotting
tmp = inputdlg({'Enter sample ID:'},'SID',1,{'M000'});
SID = tmp{1};
OUTdices = cell(size(INdices));
for j = 1:numel(INdices)
    OUTdices{j} = setdiff(1:size(lifetimes{j}.ALL,1),INdices{j});
    OUTdices{j} = reshape(OUTdices{j},[],1);
end

for j = 1:Nspecies    
    % write .txt file IN
    ID = [SID '_IN_' Names{j}];
    wave_names = {[ID '_mTb'],[ID '_mTu'],[ID '_SEMb'],[ID '_SEMu'], [ID '_sample'], [ID '_movie']};
    fileID=fopen(ID, 'w'); %open file to write
    for i=1:6            %write wavenames at each column header
        fprintf(fileID, [wave_names{i} '\t']);
    end
    fprintf(fileID,'\n');
    fclose(fileID);
    % append data
    dlmwrite(ID, lifetimes{j}.ScatDat(INdices{j},:), 'delimiter', '\t','-append')
    
    % write .txt file OUT
    ID = [SID '_OUT_' Names{j}];
    wave_names = {[ID '_mTb'],[ID '_mTu'],[ID '_SEMb'],[ID '_SEMu'], [ID '_sample'], [ID '_movie']};
    fileID=fopen(ID, 'w'); %open file to write
    for i=1:6            %write wavenames at each column header
        fprintf(fileID, [wave_names{i} '\t']);
    end
    fprintf(fileID,'\n');
    fclose(fileID);
    % append data
    dlmwrite(ID, lifetimes{j}.ScatDat(OUTdices{j},:), 'delimiter', '\t','-append')    
end
display('Done writing text file for scatter plot of main population')
%% Error estimate by bootstrapping: Get distribution
%
if ~exist('Nspecies','var')
    Nspecies = numel(lifetimes);
end
bootstat = cell(Nspecies,2);
bootsam = cell(Nspecies,2);
bootkhat = cell(Nspecies,1);
Nsam = 1e4;
ax = cell(2,1);
%
for j = 1:Nspecies
    LT = {vertcat(lifetimes{j}.ALL{INdices{j},1}) vertcat(lifetimes{j}.ALL{INdices{j},2})};
    Tmax = max(max(LT{1}),max(LT{2}))+0.1;
    for k = 2:-1:1
        [bootstat{j,k},bootsam{j,k}] = bootstrp(Nsam,@mean,LT{k});
    end
    display(['Species #' num2str(j) ' of ' num2str(Nspecies) ': Bootstrapping done. Getting corrected rates...'])
    bootkhat{j} = zeros(Nsam,2);
    for i = 1:Nsam
        bootkhat{j}(i,:) = get_corrected_rates({LT{1}(bootsam{j,1}(:,i)) LT{2}(bootsam{j,2}(:,i))},Tmin,Tmax);
    end
    display(['Species #' num2str(j) ' of ' num2str(Nspecies) ': Obtained corrected rates - moving on...'])
    figure
    for k = 1:2
        subplot(2,2,k)
        histogram(1./bootkhat{j}(:,k))
        title([Names{j} ', state ' num2str(k) ', corrected'], 'FontSize', 12)
        ax{1} = gca;
        subplot(2,2,k+2)
        histogram(bootstat{j,k})
        title([Names{j} ', state ' num2str(k) ', mean'], 'FontSize', 12)
        ax{2} = gca;
        for i = 1:2
            ax{i}.XLim = [min(ax{1}.XLim(1),ax{2}.XLim(1)) max(ax{1}.XLim(2),ax{2}.XLim(2))];
        end
    end
end
display('***************************************************')
display('All done.')
display('***************************************************')

%% Save bootstrap data
save rates_bootstrap.mat bootsam bootstat bootkhat

%% Get mean value and 3*sigma error
%{
pctlim = [0.005 0.995];
CIs = cell(1,2);
for k = 1:2
    CIs{k} = zeros(Nspecies,2);
    for j = 1:Nspecies
        tmp = sort(bootkhat{j}(:,k));
        CIs{k}(j,1) = tmp(ceil(numel(tmp)*pctlim(1)));
        CIs{k}(j,2) = tmp(floor(numel(tmp)*pctlim(2)));
    end
end
%}
if ~exist('SID','var')
    tmp = inputdlg({'Enter sample ID:'},'SID',1,{'M000'});
    SID = tmp{1};
end
if ~exist('ID','var')
    ID = inputdlg({'Enter tether ID'},'ID',1,{'L0'});
    ID = ID{1};
end
% make .txt file for Igor Pro plotting
waveNames = {['koffboot' ID], ...
            ['konboot' ID], ...
            ['ekoffboot' ID], ...
            ['ekonboot' ID], ...
            ['Keqboot' ID], ...
            ['eKeqboot' ID]};
datArray = zeros(Nspecies,numel(waveNames));
errMult = 3;
for j = 1:Nspecies
    for k = 1:2
        datArray(j,k) = mean(bootkhat{j}(:,k));
        datArray(j,k+2) = errMult*std(bootkhat{j}(:,k));
    end
end
datArray(:,5) = datArray(:,2)./datArray(:,1);
datArray(:,6) = datArray(:,5).*sqrt((datArray(:,3)./datArray(:,1)).^2+ ...
                                    (datArray(:,4)./datArray(:,2)).^2);

fileID = fopen([SID '_bootRates.txt'],'w');
for i = 1:numel(waveNames)
    fprintf(fileID, [waveNames{i} '\t']);
end
fprintf(fileID,'\n');
fclose(fileID);
dlmwrite([SID '_bootRates.txt'],datArray, 'delimiter', '\t', '-append')
%}
