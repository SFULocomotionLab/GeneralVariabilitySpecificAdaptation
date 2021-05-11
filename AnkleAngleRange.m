function [L_AnkleRange, R_AnkleRange, L_ankleAngle_timelock_resampled, R_ankleAngle_timelock_resampled, L_gaitcycle, R_gaitcycle] = AnkleAngleRange(names, Fs, allDatai, jj, ii, tri,tauL)

    % convert to degrees
    Lankleangle = allDatai(:,strcmp(names,'posaL')).*(180/pi);
    Rankleangle = allDatai(:,strcmp(names,'posaR')).*(180/pi);
    time = allDatai(:,end);
    
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
    n = 2;
    Wn = 15/(Fs/2);
    [b,a] = butter(n,Wn,'low');

    Mx1 = filtfilt(b,a,Mx1);
    Fz1 = filtfilt(b,a,Fz1);
    Mx2 = filtfilt(b,a,Mx2);
    Fz2 = filtfilt(b,a,Fz2);
    
    Mx1 = Mx1(~isnan(Mx1));
    Fz1 = Fz1(~isnan(Fz1));
    Mx2 = Mx2(~isnan(Mx2));
    Fz2 = Fz2(~isnan(Fz2));
    
    %% Time lock - use HS/TO to get 100% of gait cycle
    
    thresLow = mean(Fz1)/3; thresHigh = mean(Fz1);
    waitHS1 = 1; waitHS2 = 1; waitTO1 = 1; waitTO2 = 1;
    kHS1 = 1; kHS2 = 1; kTO1 = 1; kTO2 = 1;
    
    for i = 2:length(Fz1)
        % LEFT
        % wait 0.25 seconds between steps, steep positive delta
        delta1 = Fz1(i) - Fz1(i-1);
        if Fz1(i) > thresLow && Fz1(i) < thresHigh 
            if waitHS1 > 250 && waitTO1 > 150
                if delta1 > 0 && kHS1 == kTO1
                    HS1(kHS1) = i;
                    HS1_log(i) = 1;
                    kHS1 = kHS1+1;
                    waitHS1 = 1;
                end
            end
        end
        if Fz1(i) > thresLow && Fz1(i) < thresHigh 
            if waitTO1 > 250 && waitHS1 > 150
                if delta1 < 0 && kTO1 == kHS1-1
                    TO1(kTO1) = i+30;
                    TO1_log(i+30) = 1;
                    kTO1 = kTO1+1;
                    waitTO1 = 1;
                end
            end
        end
        waitHS1 = waitHS1+1; waitTO1 = waitTO1+1;
        
        % RIGHT
        % wait 0.25 seconds between steps
        delta2 = Fz2(i) - Fz2(i-1);
        if Fz2(i) > thresLow && Fz2(i) < thresHigh 
            if waitHS2 > 250 && waitTO2 > 150
                if delta2 > 0 && kHS2 == kTO2
                    HS2(kHS2) = i;
                    HS2_log(i) = 1;
                    kHS2 = kHS2+1;
                    waitHS2 = 1;
                end
            end
        end
        if Fz2(i) > thresLow && Fz2(i) < thresHigh 
            if waitTO2 > 250 && waitHS2 > 150
                if delta2 < 0 && kTO2 == kHS2-1
                    TO2(kTO2) = i+30;
                    TO2_log(i+30) = 1;
                    kTO2 = kTO2+1;
                    waitTO2 = 1;
                end
            end
        end
        waitHS2 = waitHS2+1; waitTO2 = waitTO2+1;
        
    end
    
    if 0
        figure
        hold on
        plot(HS1_log)
        plot(TO1_log)
        hold on
        plot(Fz1)
        
        figure
        hold on
        plot(HS2_log)
        plot(TO2_log)
        hold on
        plot(Fz2)
    end
    
    %% LEFT
    
    % time lock to stance (from left HS to left TO)
    L_ankleAngle_timelock = nan(1500,(length(HS1)-1));
    L_ankleTime_timelock = nan(1500,(length(HS1)-1));
    
    if length(HS1) > length(TO1)
        tLength = length(TO1);
    else
        tLength = length(HS1);
    end
    
    deltai = abs(TO1(1:tLength)-HS1(1:tLength));
    deltai = deltai(~isoutlier(deltai));
    deltai = deltai(deltai<500);
    upperLim = round(max(deltai))+1;
    lowerLim = round(min(deltai))-1;
    tLii = 1;
    for tLi = 1:tLength-1
        if ~isempty(TO1(TO1 > HS1(tLi) & TO1 < HS1(tLi+1))) % step
            L_ankleTOi = TO1(TO1 > HS1(tLi) & TO1 < HS1(tLi+1)); L_ankleTOi = L_ankleTOi(end);
            L_ankle_init = Lankleangle(HS1(tLi):L_ankleTOi);
            
            if length(L_ankle_init) > lowerLim && length(L_ankle_init) < upperLim && Fz1(HS1(tLi)+100) > mean(Fz1)
                L_length(tLii) = length(L_ankle_init);

                % range
                temp1 = max(L_ankle_init(1:round(length(L_ankle_init)/2)));
                temp2 = min(L_ankle_init(round(length(L_ankle_init)/2):end));
                temp3 = max(L_ankle_init(round(length(L_ankle_init)/2):end));
                
                L_AnkleRange(tLii) = range(L_ankle_init);
%                 L_AnkleRange(tLii) = range([temp1 temp2]);
%                 L_AnkleRange(tLii) = range([temp2 temp3]);
%                 L_AnkleRange(tLii) = temp2;

                % log angle
                L_ankleAngle_timelock(1:L_length(tLii),tLii) = L_ankle_init;
                L_ankleTime_timelock(1:L_length(tLii),tLii) = time(HS1(tLi):L_ankleTOi) - time(HS1(tLi));

                tLii = tLii+1;
                clearvars L_ankle_init L_ankleTOi
            end
        end
    end
    for tL = 1:(length(L_length))
        t = linspace(1,length(L_ankleAngle_timelock(1:L_length(tL),tL)),100);
        tinit = 1:length(L_ankleAngle_timelock(1:L_length(tL),tL));
        L_ankleAngle_timelock_resampled(:,tL) = interp1(tinit,L_ankleAngle_timelock(1:L_length(tL),tL),t);
        L_ankleTime_timelock_resampled(:,tL) = interp1(tinit,L_ankleTime_timelock(1:L_length(tL),tL),t);
        L_gaitcycle(:,tL) = linspace(0,100,100)';
    end
    L_ankle_mean = mean(L_ankleAngle_timelock_resampled');
    
    %% RIGHT
    
    % time lock to stance (from left HS to left TO)
    R_ankleAngle_timelock = nan(1500,(length(HS2)-1));
    R_ankleTime_timelock = nan(1500,(length(HS2)-1));
    
    if length(HS2) > length(TO2)
        tRength = length(TO2);
    else
        tRength = length(HS2);
    end
    
    deltai = abs(TO2(1:tRength)-HS2(1:tRength));
    deltai = deltai(~isoutlier(deltai));
    deltai = deltai(deltai<500);
    upperLim = round(max(deltai))+1;
    lowerLim = round(min(deltai))-1;
    tRii = 1;
    for tRi = 1:tRength-1
        if ~isempty(TO2(TO2 > HS2(tRi) & TO2 < HS2(tRi+1)))
            R_ankleTOi = TO2(TO2 > HS2(tRi) & TO2 < HS2(tRi+1)); R_ankleTOi = R_ankleTOi(end);
            R_ankle_init = Rankleangle(HS2(tRi):R_ankleTOi);
            
            if length(R_ankle_init) > lowerLim && length(R_ankle_init) < upperLim && Fz2(HS2(tRi)+100) > mean(Fz2)
            
                R_length(tRii) = length(R_ankle_init);

                % range
                temp1 = max(R_ankle_init(1:round(length(R_ankle_init)/2)));
                temp2 = min(R_ankle_init(round(length(R_ankle_init)/2):end));
                temp3 = max(R_ankle_init(round(length(R_ankle_init)/2):end));
                
                R_AnkleRange(tRii) = range(R_ankle_init);
%                 R_AnkleRange(tRii) = range([temp1 temp2]);
%                 R_AnkleRange(tRii) = range([temp2 temp3]);
%                 R_AnkleRange(tRii) = temp2;

                % log angle
                R_ankleAngle_timelock(1:R_length(tRii),tRii) = R_ankle_init;
                R_ankleTime_timelock(1:R_length(tRii),tRii) = time(HS2(tRi):R_ankleTOi) - time(HS2(tRi));

                tRii = tRii+1;
                clearvars R_ankleTOi R_ankle_init 
            end
        end
    end
    for tR = 1:(length(R_length))
        t = linspace(1,length(R_ankleAngle_timelock(1:R_length(tR),tR)),100);
        tinit = 1:length(R_ankleAngle_timelock(1:R_length(tR),tR));
        R_ankleAngle_timelock_resampled(:,tR) = interp1(tinit,R_ankleAngle_timelock(1:R_length(tR),tR),t);
        R_ankleTime_timelock_resampled(:,tR) = interp1(tinit,R_ankleTime_timelock(1:R_length(tR),tR),t);
        R_gaitcycle(:,tR) = linspace(0,100,100)';
        
    end
    R_ankle_mean = mean(R_ankleAngle_timelock_resampled');

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
        if strcmp(tri,'ga1') || strcmp(tri,'ga2') || strcmp(tri,'zt1') || strcmp(tri,'zt2')
            figure(jj+10^2)
            if strcmp(tri,'ga1')
                subplot(4,8,ii)
                col = SubjColors(jj,:);
                triname = 'GA1: L/R Ankle';
            elseif strcmp(tri,'ga2')
                subplot(4,8,ii+2*8)
                col = SubjColors(jj,:);
                triname = 'GA2: L/R Ankle';
            elseif strcmp(tri,'zt1')
                subplot(4,8,ii+8)
                col = [0.655, 0.725, 0.80];
                triname = 'ZT1: L/R Ankle';
            elseif strcmp(tri,'zt2')
                subplot(4,8,ii+3*8)
                col = [0.655, 0.725, 0.80];
                triname = 'ZT2: L/R Ankle';
            end
          
            hold on
            plot(L_gaitcycle,L_ankleAngle_timelock_resampled,'Color',col)
            plot([L_gaitcycle(1) L_gaitcycle(end)],[mean(L_AnkleRange) mean(L_AnkleRange)],'-k')
            
            plot(R_gaitcycle,R_ankleAngle_timelock_resampled,'Color',col)
            plot([R_gaitcycle(1) R_gaitcycle(end)],[mean(R_AnkleRange) mean(R_AnkleRange)],'-k')
            
            ylim([-40 60])
            title(strcat('S',string(jj),{' '},triname))
            xlabel('Time in exo (mins)')
            ylabel('Ankle (deg)')
        end
    end
end