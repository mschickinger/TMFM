%%
tau_err = cell(3,2);
k_err = cell(3,2);
foo = {L2;L6;L10};
for s = 1:2
    for i = 1:3
        tau_err{i,s} = [foo{i}.taus(:,s)';foo{i}.SEM(:,s)'];
        k_err{i,s} = [1./foo{i}.taus(:,s)'; zeros(1,length(foo{i}.taus))];
        e_rel = foo{i}.SEM(:,s)./foo{i}.taus(:,s);
        k_err{i,s}(2,:) = k_err{i,s}(1,:).*e_rel';
    end
end

K_err = k_err(:,1);
for i = 1:3
    K_err{i}(1,:) = (foo{i}.taus(:,1)./foo{i}.taus(:,2))';
    erel = sqrt((foo{i}.SEM(:,1)./foo{i}.taus(:,1)).^2 + (foo{i}.SEM(:,2)./foo{i}.taus(:,2)).^2);
    K_err{i}(2,:) = K_err{i}(1,:).*erel';
end

%% Lebensdauern
mult_err = 6;
figure('Units','normalized','Position',[0 0.5 1 0.5])
for s = 1:2
    subplot(1,3,s)
    hold off
    for i = 1:3
        C_plot = C((1+double(i==3)):3);      
        errorbar(C_plot,tau_err{i,s}(1,:),mult_err*tau_err{i,s}(2,:),'-x','MarkerSize',8,'LineWidth',0.75,'AlignVertexCenters','on')
        hold on
    end
    set(gca,'Xscale','log','XLim',[0.1 1.1],'YLim',[5 120])
end

%% Raten
mult_err = 6;
figure('Units','normalized','Position',[0 0.5 1 0.5])
for s = 1:2
    subplot(1,3,s)
    hold off
    for i = 1:3
        C_plot = C((1+double(i==3)):3);
        
        errorbar(C_plot,k_err{i,s}(1,:),mult_err*k_err{i,s}(2,:),'-x','MarkerSize',8,'LineWidth',0.75,'AlignVertexCenters','on')
        hold on
    end
    set(gca,'Xscale','log','XLim',[0.1 1.1],'YLim',[0 0.15])
end
%%
subplot(1,3,3)
hold off
for i = 1:3
    C_plot = C((1+double(i==3)):3);
    errorbar(log(C_plot),log(K_err{i}(1,:)),mult_err*K_err{i}(2,:),'-x','MarkerSize',8,'LineWidth',0.75,'AlignVertexCenters','on')
    hold on
end
%set(gca,'Xscale','log','Yscale','log','XLim',[0.1 1.1],'YLim',[0.1 3])
set(gca,'XLim',[-2 0.1],'YLim',[-1.75 1.25])