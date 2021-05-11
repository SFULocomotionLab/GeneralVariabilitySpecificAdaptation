function pl5 = Fig5_Timescales(yStd, timepointsVAL, SubjColors, formatGraphs, n, orig)

    % plot specifics
    if formatGraphs == 1
        dotVal = 8;
        if n == 1
            ylims = [0 3; -2 6; 0 3; 0 3];
        elseif n == 2
            ylims = [0.8 1.2; 0 4; 0.6 1.6; 0.6 1.6];
        end
    elseif formatGraphs == 0 && n == 1
        titles = {'A. Step frequency var';'B. Ankle angle var';'C. Soleus var';'D. Gastroc var'};
        ylabels = {'SF var (normalized to NW)';'Ankle angle var (normalized to ZT)';'Soleus var (normalized to NW)';'Medial gastroc var (normalized to NW)'};
        dotVal = 15;
        ylims = [0 3; -2 6; 0 3; 0 3];
    elseif formatGraphs == 0 && n == 2
        titles = {'A. Step frequency adapt';'B. Ankle angle adapt';'C. Soleus adapt';'D. Gastroc adapt'};
        ylabels = {'SF adapt (normalized to NW)';'Ankle angle adapt (normalized to ZT)';'Soleus adapt (normalized to NW)';'Medial gastroc adapt (normalized to NW)'};
        dotVal = 15;
        ylims = [0.8 1.2; 0 4; 0.6 1.6; 0.6 1.6];
    end
    
    % exponential equations
    modelFun = @(b,t)b(1)*exp(-t/b(2)) + b(3); % + b(3)*exp(-t/b(4)) + b(5);
    
    % subjects for mixed effects model
    Subji = (1:10).*ones(12,10);
    
    % loop thru dimensions
    for j = 1:4 % for GA, OPT
        pl5 = figure(13948);
        if n == 1
            subplot(2,4,j+4)
        elseif n == 2
            subplot(2,4,j)
        end
        hold on
        
        % data
        clearvars xi yi x y
        yi = yStd{j,2};
        xi = timepointsVAL;
        
        % initial fit
        if j == 2 && n == 2 % skip ankle adapt
            merp = 1;
            y = yi(~isnan(yi));
        elseif orig == 0 && j == 1 && n == 1
            merp = 1;
            y = yi(~isnan(yi));
        else
            % guess
            if j <= 2
                start = [10 600 1];
            else
                start = [1 5 1];
            end

            % fit
            Subj = Subji(~isnan(yi));
            x = xi(~isnan(yi));
            y = yi(~isnan(yi));
            [betaInit(:,j),~,stats,B] = nlmefit(x,y,Subj,[],modelFun,start,'REParamsSelect',3);
            
            % plot
            t = linspace(min(x), max(x),100);
            plot(t,modelFun(betaInit(:,j),t),'k','LineWidth',2)
            
            % CI
            xrange = linspace(min(x),max(x),100);
            [ypred,delta] = nlpredci(modelFun,xrange,betaInit(:,j),stats.ires,'Covar',stats.covb); %,'MSE',stats.mse,'SimOpt','on');
            lower = ypred - delta;
            upper = ypred + delta;
            plot(xrange,lower,'-k')
            plot(xrange,upper,'-k')
            
            % R-squared
            R2(j) = sum((stats.ires).^2)/sum((y - mean(y)).^2);
            
            % fit residuals
            res = stats.ires;
            
            % subtract random effects
            yi = yi-B;
        end
        
        % plot data
        for i=1:size(yi,2)
            if i <=5
                colo = [0 0.5000 0.5000];
            else
                colo = [0.9570 0.5078 0.1875];
            end

            % plotting validation
            plot(xi(:,i),yi(:,i),'.','MarkerSize',dotVal,'Color',SubjColors(i,:))

        end
        
        % plot specifics
        ylim(ylims(j,:))
        if n == 1
            if j == 2
                yticks(ylims(j,1):2:ylims(j,2))
            else
                yticks(ylims(j,1):1:ylims(j,2))
            end
        elseif n == 2
            if j == 1
                yticks(ylims(j,1):0.1:ylims(j,2))
            elseif j == 2
                yticks(ylims(j,1):1:ylims(j,2))
            else
                yticks(ylims(j,1):0.2:ylims(j,2))
            end
        end
        xlim([0 600])
        xticks(0:200:600)
        if formatGraphs == 1
            set(gcf,'units','inches','position',[0,0,6,4.08])
            set(gca,'FontSize',5)
        else
            title(strcat(titles(j)))
            xlabel('Time in exoskeleton (mins)')
            ylabel(ylabels(j))
        end
        
    end
    
    % table for time constants
    Rows = {'offset';'multiplier';'time constant';'R-squared'};
    SF = [betaInit(3,1); betaInit(1,1); betaInit(2,1); R2(1)];
    Ankle = [betaInit(3,2); betaInit(1,2); betaInit(2,2); R2(2)];
    Soleus = [betaInit(3,3); betaInit(1,3); betaInit(2,3); R2(3)];
    Gastroc = [betaInit(3,4); betaInit(1,4); betaInit(2,4); R2(4)];
    timeconstants = table(Rows,SF,Ankle,Soleus,Gastroc);
    disp(timeconstants)
    
end