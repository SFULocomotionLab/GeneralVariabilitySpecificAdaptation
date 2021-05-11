function pl2 = Fig2_Exploration_BarGraph_RepSubj(SFstdVAL_nw, SFstdVAL_zt, SFstdVAL, SFmeanVAL_nw, AnkleStdVAL_zt, AnkleVAL_zt, AnkleStdVAL, AnkleVAL, solTotalStdVAL_nw, solTotalStdVAL, solTotalMeanVAL_nw, gastrocTotalStdVAL_nw, gastrocTotalStdVAL, gastrocTotalMeanVAL_nw, orig, formatGraphs, SubjColors, Fs_timescale)

    % Figure 2: Bar plot of 
    % without-assistancs day 1 vs. with-assistance day 1
    % with-assistance day 1 vs. with-assistance day 6

    % loop through if sensitivity analysis
    for s = 1:length(Fs_timescale)
        
        % loop through variables
        varsbaseline = {'SFstdVAL_nw';'AnkleStdVAL_zt';'solTotalStdVAL_nw';'gastrocTotalStdVAL_nw'};
        vars = {'SFstdVAL';'AnkleStdVAL';'solTotalStdVAL';'gastrocTotalStdVAL'};
        for i = 1:4
            
            % each variable
            ybaseline = eval(string(varsbaseline(i)));
            y = eval(string(vars(i)));
            
            % bpm to spm
            if i == 1 && orig == 1
                ybaseline = ybaseline./2;
                y = y./2;
            end
            
            % normal walking
            % average across this first day
            if i > 2 % muscle activity
                ybaseline_nw = nanmean(ybaseline(1:2,1:10,s));
            else % other vars
                ybaseline_nw = nanmean([ybaseline(1:2,1:5,s) ybaseline(3:4,6,s) ybaseline(1:2,7:10,s)]);
            end
            ybaselineMean = nanmean(ybaseline_nw);
            ybaselineStd = nanstd(ybaseline_nw);
            
            % check for changes in normal walking variability across experiments
            ybaselineFinal = nanmean([ybaseline(9:10,1,s) ybaseline(end-1:end,2:10,s)]);
            [~,pVar_nw(i)] = ttest(ybaseline_nw,ybaselineFinal);
            percChange_nw = ((ybaselineFinal-ybaseline_nw)./ybaselineFinal).*100;
            percChange_nwmean = nanmean(percChange_nw);
            percChange_nwstd = nanstd(percChange_nw);
            
            % zero torque for step width variable
            if i == 1 && orig == 0
                ybaseline_zt = nanmean([SFstdVAL_zt(1:2,1:5,s) SFstdVAL_zt(3:4,6,s) SFstdVAL_zt(1:2,7:10,s)]);
                ybaselineMean_zt = nanmean(ybaseline_zt);
                ybaselineStd_zt = nanstd(ybaseline_zt);
                
                % test between NW and ZT
                [~,pVar_nwzt] = ttest(ybaseline_nw,ybaseline_zt);
                percChange_nwzt = ((ybaseline_zt-ybaseline_nw)./ybaseline_nw).*100;
                percChange_nwzt_mean = nanmean(percChange_nwzt);
                percChange_nwzt_std = nanstd(percChange_nwzt);
            end

            % initial avg variability GA across two six-minute trials
            if i > 2 % muscle activity
                yInit = nanmean(y(1:2,1:10,s));
            else % other vars
                yInit = nanmean([y(1:2,1:5,s) y(3:4,6,s) y(1:2,7:10,s)]);
            end
            yMeanInit = nanmean(yInit);
            yStdInit = nanstd(yInit);

            % final avg variability GA across two six-minute trials
            yFinal = nanmean([y(9:10,1,s) y(end-1:end,2:10,s)]);
            yMeanFinal = nanmean(yFinal);
            yStdFinal = nanstd(yFinal);

            % nw vs day 1 variability stats
            [~,pval_nw_GA1(s,i)] = ttest(yInit,ybaseline_nw,'Tail','right');
            percChange_nw_GA1 = ((yInit-ybaseline_nw)./ybaseline_nw).*100;
            percChange_nw_GA1_mean(s,i) = nanmean(percChange_nw_GA1);
            percChange_nw_GA1_std(s,i) = nanstd(percChange_nw_GA1);

            % day 1 vs day 6 variability stats
            [~,pval_GA1_GA5(s,i)] = ttest(yFinal,yInit,'Tail','left');
            percChange_GA1_GA5 = ((yFinal-yInit)./yInit).*100;
            percChange_GA1_GA5_mean(s,i) = nanmean(percChange_GA1_GA5);
            percChange_GA1_GA5_std(s,i) = nanstd(percChange_GA1_GA5);

            % nw vs day 6 variability stats
            [~,pval_nw_GA5(s,i)] = ttest(yFinal,ybaseline_nw);
            percChange_nw_GA5 = ((yFinal-ybaseline_nw)./ybaseline_nw).*100;
            percChange_nw_GA5_mean(s,i) = nanmean(percChange_nw_GA5);
            percChange_nw_GA5_std(s,i) = nanstd(percChange_nw_GA5);
            
            % plotting
            pl2 = figure(23481);
            subplot(length(Fs_timescale),4,s+3*(s-1)+(i-1));
            hold on
            if i == 1 && orig == 0
                b = bar([1 2 3 4],[ybaselineMean; yMeanInit; yMeanFinal; ybaselineMean_zt]);
                plot((4+0.1.*randn(length(ybaseline_zt),1)),ybaseline_zt,'k.','MarkerSize',10)
                errorbar([ybaselineMean; yMeanInit; yMeanFinal; ybaselineMean_zt],[ybaselineStd; yStdInit; yStdFinal; ybaselineStd_zt],'k.','LineWidth',1)
                b.CData(4,:) = [0.655, 0.725, 0.80];
            else
                b = bar([1 2 3],[ybaselineMean; yMeanInit; yMeanFinal]);
                errorbar([ybaselineMean; yMeanInit; yMeanFinal],[ybaselineStd; yStdInit; yStdFinal],'k.','LineWidth',1)
            end
            plot((1+0.1.*randn(length(ybaseline_nw),1)),ybaseline_nw,'k.','MarkerSize',10)
            plot((2+0.1.*randn(length(yInit),1)),yInit,'k.','MarkerSize',10)
            plot((3+0.1.*randn(length(yFinal),1)),yFinal,'k.','MarkerSize',10)
            
            b.FaceColor = 'flat';
            if i == 2
                b.CData(1,:) = [0.655, 0.725, 0.80];
            else
                b.CData(1,:) = [0.776, 0.675, 0.592];
            end
            b.CData(2,:) = [0, 0.643, 0.89];
            b.CData(3,:) = [0, 0.643, 0.89];
            
            % plot details
            if orig == 1
                ylimi = [0 4; 0 6; 0 0.06; 0 0.06];
                titlei = {'step freq';'ankle angle range';'soleus';'medial gastroc'};
                xlim([0 4])
                xticks(1:3)
                xticklabels({'NW','D1','D6'})
            else
                ylimi = [0 0.04; 0 8; 0 0.2; 0 0.2];
                titlei = {'step width';'peak ankle angle';'rectus femoris';'biceps femoris'};
                if i == 1
                    xlim([0 5])
                    xticks(1:4)
                    xticklabels({'NW','D1','D6'})
                else
                    xlim([0 4])
                    xticks(1:3)
                    xticklabels({'NW','D1','D6','ZT'})
                end
            end
            ylim(ylimi(i,:));
            title(titlei(i));
            ylabel('variability')
            set(gca,'FontSize',12)
            pbaspect([1 1 1])
        end
        
        % table of statistics
        if orig == 1 % original four vars
            Rows = {'Initial vs. nw';'Percent Change';'Percent Change std';'Initial vs. Final';'Percent Change';'Percent Change std';'nw vs. Final';'Percent Change';'Percent Change std'};
            SWVar = [pval_nw_GA1(s,1); percChange_nw_GA1_mean(s,1); percChange_nw_GA1_std(s,1); pval_GA1_GA5(s,1); percChange_GA1_GA5_mean(s,1); percChange_GA1_GA5_std(s,1); pval_nw_GA5(s,1); percChange_nw_GA5_mean(s,1); percChange_nw_GA5_std(s,1)];
            AnkleVar = [pval_nw_GA1(s,2); percChange_nw_GA1_mean(s,2); percChange_nw_GA1_std(s,2); pval_GA1_GA5(s,2); percChange_GA1_GA5_mean(s,2); percChange_GA1_GA5_std(s,2); pval_nw_GA5(s,2); percChange_nw_GA5_mean(s,2); percChange_nw_GA5_std(s,2)];
            RFVar = [pval_nw_GA1(s,3); percChange_nw_GA1_mean(s,3); percChange_nw_GA1_std(s,3); pval_GA1_GA5(s,3); percChange_GA1_GA5_mean(s,3); percChange_GA1_GA5_std(s,3); pval_nw_GA5(s,3); percChange_nw_GA5_mean(s,3); percChange_nw_GA5_std(s,3)];
            BFVar = [pval_nw_GA1(s,4); percChange_nw_GA1_mean(s,4); percChange_nw_GA1_std(s,4); pval_GA1_GA5(s,4); percChange_GA1_GA5_mean(s,4); percChange_GA1_GA5_std(s,4); pval_nw_GA5(s,4); percChange_nw_GA5_mean(s,4); percChange_nw_GA5_std(s,4)];
        else % additional four vars includes ZT for step width
            Rows = {'nw vs. zt';'Percent Change';'Percent Change std';'Initial vs. nw';'Percent Change';'Percent Change std';'Initial vs. Final';'Percent Change';'Percent Change std';'nw vs. Final';'Percent Change';'Percent Change std'};
            SWVar = [pVar_nwzt; percChange_nwzt_mean; percChange_nwzt_std; pval_nw_GA1(s,1); percChange_nw_GA1_mean(s,1); percChange_nw_GA1_std(s,1); pval_GA1_GA5(s,1); percChange_GA1_GA5_mean(s,1); percChange_GA1_GA5_std(s,1); pval_nw_GA5(s,1); percChange_nw_GA5_mean(s,1); percChange_nw_GA5_std(s,1)];
            AnkleVar = [nan; nan; nan; pval_nw_GA1(s,2); percChange_nw_GA1_mean(s,2); percChange_nw_GA1_std(s,2); pval_GA1_GA5(s,2); percChange_GA1_GA5_mean(s,2); percChange_GA1_GA5_std(s,2); pval_nw_GA5(s,2); percChange_nw_GA5_mean(s,2); percChange_nw_GA5_std(s,2)];
            RFVar = [nan; nan; nan; pval_nw_GA1(s,3); percChange_nw_GA1_mean(s,3); percChange_nw_GA1_std(s,3); pval_GA1_GA5(s,3); percChange_GA1_GA5_mean(s,3); percChange_GA1_GA5_std(s,3); pval_nw_GA5(s,3); percChange_nw_GA5_mean(s,3); percChange_nw_GA5_std(s,3)];
            BFVar = [nan; nan; nan; pval_nw_GA1(s,4); percChange_nw_GA1_mean(s,4); percChange_nw_GA1_std(s,4); pval_GA1_GA5(s,4); percChange_GA1_GA5_mean(s,4); percChange_GA1_GA5_std(s,4); pval_nw_GA5(s,4); percChange_nw_GA5_mean(s,4); percChange_nw_GA5_std(s,4)];
        end
        pvalsVar = table(Rows,SWVar,AnkleVar,RFVar,BFVar);
        disp(pvalsVar)
    end
    set(gcf,'units','inches','position',[0,0,10.5,4])
    
    %% Representative subject
    if 0 %orig == 1
        % load data
        cd(strcat('/Volumes/GoogleDrive/My Drive/Locomotion Lab/Projects/Cost Optimization - exos/Katies Data/Learning/Saving'))

        files = {'TEBX-7_18-VAL-nw1SF.mat'; ...
                 'TEBX--7_18-VAL-zt1AnkleAngleandRange.mat'; ...
                 'TEBX-7_18-VAL-nw1EMG.mat'; ...
                 'TEBX-7_18-VAL-ga1SF.mat'; ...
                 'TEBX-7_18-VAL-ga1AnkleAngleandRange.mat'; ...
                 'TEBX-7_18-VAL-ga1EMG.mat'; ...
                 'TEBX-8_4-VAL-ga1SF.mat'; ...
                 'TEBX-8_4-VAL-ga1AnkleAngleandRange.mat'; ...
                 'TEBX-8_4-VAL-ga1EMG.mat';};
        names = {'SF';'Ankle';'EMG';'SF';'Ankle';'EMG';'SF';'Ankle';'EMG'};
        
        figure(203)
        mSF = 1; mAnkle = 2; mSol = 3; mGastroc = 4;
        for k = 1:length(files)
            clearvars y t
            
            load(string(files(k)))
            
            if strcmp(names(k),'SF')
                % bpm to spm
                y = y./2;
                % time
                t = t./60; t = t-t(2);
                if length(y) ~= length(t)
                    y = y(1:length(t));
                end
                % plot
                subplot(3,4,mSF)
                hold on
                plot(t,y)
                ylim([40 80])
                xticks(0:2:6)
                mSF = mSF+4;
            elseif strcmp(names(k),'Ankle')
                % plot
                subplot(3,4,mAnkle)
                hold on
                plot(t,y)
                ylim([-20 40])
                xticks(0:50:100)
                mAnkle = mAnkle+4;
            elseif strcmp(names(k),'EMG')
                % plot sol
                subplot(3,4,mSol)
                hold on
                plot(t,y{1})
                ylim([-1 2])
                xticks(0:50:100)
                mSol = mSol+4;
                
                % plot gastroc
                subplot(3,4,mGastroc)
                hold on
                plot(t,y{2})
                ylim([-1 2])
                xticks(0:50:100)
                mGastroc = mGastroc+4;
            end
            
        end
        
    end

end