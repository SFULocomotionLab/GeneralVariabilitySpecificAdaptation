function pl3 = Fig3_MetCost_Adapt(yAdapt, metVAL, metVAL_zt, SFmeanVALnorm, AnkleVALnorm, AnkleVAL, AnkleVAL_zt, solTotalMeanVALnorm, gastrocTotalMeanVALnorm, orig, formatGraphs, SubjColors)

    % Figure 3: Met cost relationships and bar plot of adaptation

    % four variables
    % SF
    SFi = yAdapt{1,2};
    SF = SFi(:);
    % Ankle
    Anklei = yAdapt{2,2};
    Ankle = Anklei(:);
    % Soleus
    Soleusi = yAdapt{3,2};
    Soleus = Soleusi(:);
    % Gastroc
    Gastroci = yAdapt{4,2};
    Gastroc = Gastroci(:);
    % subjects for mixed effects model
    Subji = (1:10).*ones(12,10);
    Subj = Subji(:);

    % normalize met cost to first day of NW or ZT
    metNorm = nan(12,10);
    for k = 1:size(metVAL,2) % subjs
        if k == 6
            ind = 3; % subj 6 missing data on day 1
        else
            ind = 1;
        end
        metNorm(:,k) = (metVAL(:,k)-nanmean(metVAL_zt(ind:ind+1,k)))./(nanmean(metVAL_zt(ind:ind+1,k)));
    end
    Met = metNorm(:);

    % tables
    SF_tbl = table(Subj,Met,SF,'VariableNames',{'Group','MetCost','yVar'});
    Ankle_tbl = table(Subj,Met,Ankle,'VariableNames',{'Group','MetCost','yVar'});
    Soleus_tbl = table(Subj,Met,Soleus,'VariableNames',{'Group','MetCost','yVar'});
    Gastroc_tbl = table(Subj,Met,Gastroc,'VariableNames',{'Group','MetCost','yVar'});

    % four liner mixed effects models
    SF_lme = fitlme(SF_tbl,'MetCost ~ 1 + yVar + (1|Group)');
    Ankle_lme = fitlme(Ankle_tbl,'MetCost ~ 1 + yVar + (1|Group)');
    Soleus_lme = fitlme(Soleus_tbl,'MetCost ~ 1 + yVar + (1|Group)');
    Gastroc_lme = fitlme(Gastroc_tbl,'MetCost ~ 1 + yVar + (1|Group)');

    %     % multiple liner mixed effects model
    %     tbl = table(Subj, Met, SF, Ankle, Soleus, Gastroc, ...
    %         'VariableNames',{'Group','MetCost','SF','Ankle','Soleus','Gastroc'});
    %     multlme = fitlme(tbl,'MetCost ~ 1 + SF + Ankle + Soleus + Gastroc + (1|Group)');

    if orig == 1 % axis for original vars
        xlims = [0.8 1.2; 0 4; 0.4 1.6; 0.4 1.6];
    else % axis for additional vars
        xlims = [0 2; 0 5; 0 3; 0 3];
    end

    if formatGraphs == 1 % format graphs for manuscript
        dotVal = 6;
        dotAdapt = 2;
    else % regulat graphs
        dotVal = 30;
        dotAdapt = 10;
    end

    % loop through dimensions
    betaSave = nan(2,4);
    slopeCI = nan(2,4);
    Vars = {'SF','Ankle','Soleus','Gastroc'};
    VarsAdapt = {'SFmeanVALnorm','AnkleVALnorm','solTotalMeanVALnorm','gastrocTotalMeanVALnorm'};
    for i = 1:4
        pl3 = figure(360);
        ax(i) = subplot(1,4,i);
        hold on

        % variable range
        varlim = eval(string(cellstr(Vars(i)))); varlim = varlim(~isnan(Met));
        xfit = linspace(min(varlim),max(varlim),100);
        xtemp = eval(strcat(string(cellstr(Vars(i))),'i'));
        
        % mixed effects model
        beta = fixedEffects(eval(strcat(string(cellstr(Vars(i))),'_lme')));
        betaSave(:,i) = beta;
        [~,~,STATS] = randomEffects(eval(strcat(string(cellstr(Vars(i))),'_lme'))); % Compute the random-effects statistics (STATS)
        STATS.Level = nominal(STATS.Level);

        % met cost
        randEff = STATS.Estimate';
        mettemp = metNorm-randEff;

        % plotting individual fits
        for j = 1:10
            scatter(xtemp(:,j),mettemp(:,j),dotVal,'Filled','MarkerFaceColor',SubjColors(j,:))
        end

        % plotting fit
        y_hat = beta(1) + beta(2)*xfit;
        plot(xfit,y_hat,'-k','LineWidth',2)

        % plot 95 % CI
        tblnew = table();
        tblnew.Group = repmat(5,100,1);
        tblnew.yVar = xfit';
        tblnew.MetCost = linspace(min(metNorm(metNorm~=0)),max(metNorm(metNorm~=0)),100)';
        [~,yCI] = predict(eval(strcat(string(cellstr(Vars(i))),'_lme')),tblnew,'Conditional',0);
        plot(tblnew.yVar,yCI,'k-');
        
        % slope 95% CI
        lmei = eval(string(strcat(cellstr(Vars(i)),'_lme.Coefficients')));
        slopeCI(:,i) = [table2array(lmei(2,7)); table2array(lmei(2,8))];

        % R-squared and pval
        R2(i) = eval(string(strcat(cellstr(Vars(i)),'_lme.Rsquared.Ordinary')));
        pval(i) = table2array(lmei(2,6));

        % plot settings
        ylim([-0.8 0.2])
        yticks(-0.8:0.2:0.2);
        xlim(xlims(i,:))
        if orig == 1
            if i == 1
                xticks(xlims(i,1):0.1:xlims(i,2))
            elseif i == 2
                xticks(xlims(i,1):1:xlims(i,2))
            else
                xticks(xlims(i,1):0.4:xlims(i,2))
            end
        else
            xticks(xlims(i,1):1:xlims(i,2))
        end
            title({string(strcat(Vars(i),{' '},'Adaptation vs. Met Cost')),strcat('R-squared:',{' '},string(R2(i))),strcat('Pval:',{' '},string(pval(i)))})
        if formatGraphs == 1
            set(gcf,'units','inches','position',[0,0,6.25,2.25])
        else
            xlabel({string(strcat(Vars(i),{' '},'Adaptation')),'(normalized to NW/ZT)'})
            ylabel('Met Cost (W)')
            set(gca,'FontSize',15)
        end
        
        % calculating adaptation for each variable
        y = eval(string(VarsAdapt(i)));
        if i > 2 % muscle activity
            init = nanmean(y(1:2,1:10));
            fin = nanmean([y(9:10,1) y(end-1:end,2:10)]);
        else % other vars
            init = nanmean([y(1:2,1:5) y(3:4,6) y(1:2,7:10)]);
            fin = nanmean([y(9:10,1) y(end-1:end,2:10)]);
        end
        [~,pvaladapt(i)] = ttest(init,fin);
        percChange = ((fin-init)./init).*100;
        percChangeMean(i) = nanmean(percChange);
        percChangeStd(i) = nanstd(percChange);
        
        % plotting adaptation for each variable
        figure(23491)
        subplot(1,4,i)
        hold on
        b = bar(1:2,[nanmean(init); nanmean(fin)]);
        errorbar([nanmean(init); nanmean(fin)],[nanstd(init); nanstd(fin)],'k.','LineWidth',1)
        plot((1+0.1.*randn(length(init),1)),init,'k.','MarkerSize',10)
        plot((2+0.1.*randn(length(fin),1)),fin,'k.','MarkerSize',10)
        if orig == 1
            ylimi = [0.8 1.2; 0 4; 0.4 1.6; 0.4 1.6];
        else
            ylimi = [0 2; 0 5; 0 3; 0 3];
        end
        ylim(ylimi(i,:));
        xlim([0 3])
        xticks(1:2)
        xticklabels({'D1';'D6'})
        
        clearvars init fin percChange
    end
    
    % table for met cost relationships
    Rows = {'slope';'offset';'pval';'R-squared';'95% CI low';'95% CI high'};
    SF = [betaSave(2,1); betaSave(1,1); pval(1); R2(1); slopeCI(1,1); slopeCI(2,1)];
    Ankle = [betaSave(2,2); betaSave(1,2); pval(2); R2(2); slopeCI(1,2); slopeCI(2,2)];
    Soleus = [betaSave(2,3); betaSave(1,3); pval(3); R2(3); slopeCI(1,3); slopeCI(2,3)];
    Gastroc = [betaSave(2,4); betaSave(1,4); pval(4); R2(4); slopeCI(1,4); slopeCI(2,4)];
    metcostrelationships = table(Rows,SF,Ankle,Soleus,Gastroc);
    disp(metcostrelationships)
    
    % table for adaptation
    Rows = {'Initial vs. Final';'Percent Change';'Percent Change std'};
    SF = [pvaladapt(1); percChangeMean(1); percChangeStd(1)];
    Ankle = [pvaladapt(2); percChangeMean(2); percChangeStd(2)];
    Soleus = [pvaladapt(3); percChangeMean(3); percChangeStd(3)];
    Gastroc = [pvaladapt(4); percChangeMean(4); percChangeStd(4)];
    pvals = table(Rows,SF,Ankle,Soleus,Gastroc);
    disp(pvals)
    
    % met cost reduction
    metValinit = nanmean([metVAL(1:2,1:5) metVAL(3:4,6) metVAL(1:2,7:10)]);
    metValfinal = nanmean([metVAL(9:10,1) metVAL(end-1:end,2:10)]);
    percChange = ((metValfinal - metValinit)./metValinit)*100;
    costSavingsMean = nanmean(percChange);
    costSavingsStd = nanstd(percChange);
    % met cost perc ZT
    metValpercZT = ((metValfinal - nanmean(metVAL_zt))./nanmean(metVAL_zt))*100;
    costSavingsMean_percZT = nanmean(metValpercZT);
    costSavingsStd_percZT = nanstd(metValpercZT);
    [~,p_metCost] = ttest(metValfinal,nanmean(metVAL_zt),'tail','left');
    disp('met cost stats')
    disp([costSavingsMean_percZT, costSavingsStd_percZT])
    disp(p_metCost)
    
    if formatGraphs == 1
        set(gcf,'units','inches','position',[0,0,7,1.25])
    end
    
end