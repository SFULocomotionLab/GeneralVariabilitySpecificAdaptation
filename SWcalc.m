function [SW, tsave] = SWcalc(names, Fs, t, allDatai, jj, ii, tri)
    %% Calculate combined center of pressure in the fore-aft direction
    % Use GRF to calculate SW for now because I am not sure how foot
    % switched are dropping out
    % columns (variables) we care about (for now)
    % 56 = time
    % 41 = VO2, 42 = VCO2
    % 11 = hsL, 12 = hsR -- doesn't seem to have anything collected...
    % 13 = LFx, 14 = LFy, 15 = LFz
    % 16 = LMx, 17 = LMy, 18 = LMz
    % 19 = RFx, 20 = RFy, 21 = RFz
    % 22 = RMx, 23 = RMy, 24 = RMz
    
    % Define forces and moments
    % Left
    Fx1=allDatai(:,strcmp(names,'LFx'));
    Fy1=allDatai(:,strcmp(names,'LFy'));
    Fz1=allDatai(:,strcmp(names,'LFz'));
    Mx1=allDatai(:,strcmp(names,'LMx'));
    My1=allDatai(:,strcmp(names,'LMy'));
    Mz1=allDatai(:,strcmp(names,'LMz'));
    % Right
    Fx2=allDatai(:,strcmp(names,'RFx'));
    Fy2=allDatai(:,strcmp(names,'RFy'));
    Fz2=allDatai(:,strcmp(names,'RFz'));
    Mx2=allDatai(:,strcmp(names,'RMx'));
    My2=allDatai(:,strcmp(names,'RMy'));
    Mz2=allDatai(:,strcmp(names,'RMz'));

    %% FILTER

    % filter
%     n = 2;
%     Wn = 15/(Fs/2);
%     [b,a] = butter(n,Wn,'low');
% 
%     Mx1 = filtfilt(b,a,Mx1);
%     Fz1 = filtfilt(b,a,Fz1);
%     Mx2 = filtfilt(b,a,Mx2);
%     Fz2 = filtfilt(b,a,Fz2);
%     
%     Mx1 = Mx1(~isnan(Mx1));
%     Fz1 = Fz1(~isnan(Fz1));
%     Mx2 = Mx2(~isnan(Mx2));
%     Fz2 = Fz2(~isnan(Fz2));
    
    %% Step width
    
    % calibration
    C1 = 500; C2 = 500; C3 = 1000;
    C4 = 800; C5 = 400; C6 = 400;
    
    CM = [C1 C2 C3 C4 C5 C6];
    F1i = [Fx1 Fy1 Fz1 Mx1 My1 Mz1]; F1 = F1i.*CM;
    Fx1 = F1(:,1); Fy1 = F1(:,2); Fz1 = F1(:,3);
    Mx1 = F1(:,4); My1 = F1(:,5); Mz1 = F1(:,6);
    
    F2i = [Fx2 Fy2 Fz2 Mx2 My2 Mz2]; F2 = F2i.*CM;
    Fx2 = F2(:,1); Fy2 = F2(:,2); Fz2 = F2(:,3);
    Mx2 = F2(:,4); My2 = F2(:,5); Mz2 = F2(:,6);
            
    % CHECK TREADMILL DIMENSIONS
    % Calculate COPx for left side
    copx1 = -0.279*2 - (0.0025.*Fx1-My1)./Fz1;
    % Calculate COPx for right side
    copx2 = 0.279*2 - (0.0025.*Fx2-My2)./Fz2;

    % Deal with bad centers of pressure for low forces
    if isnan(copx1)
        copx1 = 0;
    end
    if isnan(copx2)
        copx2 = 0;
    end

    % Weighted average to deal with stepping on same belts for consecutive steps
    % Treadmill parameters - TO BE CHANGED TO COLLINS LAB DIMENSIONS
    copxi = ((copx1.*Fz1) + (copx2.*Fz2))./(Fz1+Fz2);
    
    % For now, I will estimate the center of treadmill to be:
    % using running average of 10s windows
    copx_movingmean = movmean(copxi,10*Fs);
    copx = copxi-copx_movingmean;

    %% Step frequency
    
    % Calculate COPy for left side
    copy1 = Mx1./Fz1;
    % Calculate COPy for right side
    copy2 = Mx2./Fz2;

    % Deal with bad centers of pressure for low forces
    if isnan(copy1)
        copy1 = 0;
    end
    if isnan(copy2)
        copy2 = 0;
    end

    % Weighted average to deal with stepping on same belts for consecutive steps
    % Treadmill parameters - TO BE CHANGED TO COLLINS LAB DIMENSIONS
    copyi = ((copy1.*Fz1) + (copy2.*Fz2))./(Fz1+Fz2);

    % For now, I will estimate the center of treadmill to be:
    % using running average of 10s windows
    copy_movingmean = movmean(copyi,10*Fs);
    copy = copyi-copy_movingmean;

    % refilter the COP data from 0.3-3 Hz
    n = 3;
    Wn = [0.3 3]/(Fs/2);
    [b,a] = butter(n,Wn,'bandpass');
    copy = filtfilt(b,a,copy);

    % differentiate COP to apply heelstrike conditions
    % initialize
    yp_dif = 0; yp_dift = 0; yp = copy(1); cop_step = 0;
    for j = 2:length(copy)
        yp_dif(j) = copy(j) - copy(j-1);
        yp_dift(j) = yp_dif(j-1);
        yp(j) = copy(j);

        if yp_dif(j) >= 0 && yp_dift(j) < 0 && yp(j) < 0
            cop_step(j) = 1;
        else
            cop_step(j) = 0;
        end

    end
    
    % check
    if 0
        figure
        hold on
        plot(copx)
        plot(cop_step)
    end
    
    %% Calculate SF
    
    k = 1; % initialize
    for l = 2:length(cop_step)
        if cop_step(l) == 1
            step_save(k) = l;
            k = k+1;
        end
    end

    SW = 0; tsave = 0; e = 2; % initialize
    for r = 2:length(copy)
        if r == step_save(e)
            % take difference between time in second, then multiply
            % by 60 to get bpm
            SW(e) = abs(copx(step_save(e-1))) + abs(copx(step_save(e)));
            tsave(e)  = t(step_save(e));
            if e < length(step_save)
                e = e+1;
            end
        end
    end
    
    I1 = SW>1;
    I2 = SW<=0;
    SW(I1) = nan;
    SW(I2) = nan;
    
    %% PLOT
    if 0
        SubjColors = [0         0         0;... % TA: 0	0	0
                  0.5977    0.5977    0.5977;... % TB: 153	153	153
                       0    0.5000    0.5000;... % TC: 0	128	128
                  0.7500    0.7500         0;... % TD: 192	192	0
                  0.5000         0         0;... % TE: 128	0	0
                  0.8984    0.0977    0.2930;... % PA: 230	25	75
                  0.9766    0.7422    0.7422;... % PB: 250	190	190
                  0.6445    0.6133    0.4922;... % PC: 165	157	126
                  0.9570    0.5078    0.1875;... % PD: 245	130	48
                  0.4766    0.8242    0.6680;... % PE: 122	211	171
                  0.2344    0.7031    0.2930;... % EA: 60	180	75
                       0    0.5078    0.7812;... % EB: 0	130	200
                  0.5664    0.1172    0.7031;... % EC: 145	30	180
                  0.8984    0.7422    0.9961;... % ED: 230	190	255
                  0.9961         0    0.5469]; % EE: 255	0	140
        if strcmp(tri,'ga1') || strcmp(tri,'ga2') || strcmp(tri,'nw1') || strcmp(tri,'nw2')
            figure(jj+10^1)
            if strcmp(tri,'ga1')
                subplot(4,8,ii)
                col = SubjColors(jj,:);
                triname = 'GA1: SW';
            elseif strcmp(tri,'ga2')
                subplot(4,8,ii+2*8)
                col = SubjColors(jj,:);
                triname = 'GA2: SW';
            elseif strcmp(tri,'nw1')
                subplot(4,8,ii+8)
                col = [0.776, 0.675, 0.592];
                triname = 'NW1: SW';
            elseif strcmp(tri,'nw2')
                subplot(4,8,ii+3*8)
                col = [0.776, 0.675, 0.592];
                triname = 'NW2: SW';
            end
            hold on
            plot(tsave,SW,'Color',[col 0.5])
            plot([tsave(1) tsave(end)],[nanmean(SW) nanmean(SW)],'-k')
            ylim([0 0.3])
            title(strcat('S',string(jj),{' '},triname))
            xlabel('Time in exo (mins)')
            ylabel('SW (bpm)')
        end
    end
end