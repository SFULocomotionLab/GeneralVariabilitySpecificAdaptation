function [R_gait_cycle,R_ankle_sol_timelock_resampled,R_ankle_gastroc_timelock_resampled, ...
            L_gait_cycle,L_ankle_sol_timelock_resampled,L_ankle_gastroc_timelock_resampled,...
            Lsolpeak,Rsolpeak,Lgastrocpeak,Rgastrocpeak, ...
            RtotalsolEMG,LtotalsolEMG,RtotalgastrocEMG,LtotalgastrocEMG] = EMGCalc_RFBF(Fs, EMGanalysisAll, LsolpeakNorm, RsolpeakNorm, LgastrocpeakNorm, RgastrocpeakNorm, jj, ii, tri)

    % EMG
    L_EMG_sol = EMGanalysisAll(:,1); R_EMG_sol = EMGanalysisAll(:,2);
    L_EMG_gastroc = EMGanalysisAll(:,3); R_EMG_gastroc = EMGanalysisAll(:,4);
    % forces
    Mx1 = EMGanalysisAll(:,5); Fz1 = EMGanalysisAll(:,6);
    Mx2 = EMGanalysisAll(:,7); Fz2 = EMGanalysisAll(:,8);
    % time
    time = 0:(1/Fs):length(L_EMG_sol)/Fs;
    
    % Steele 2017
    % EMG data were collected at 2000 Hz with an on-board bandpass filter 
    % applied with cut-offs at 20-450 Hz. The EMG data were then high-pass 
    % filtered at 40 Hz (3rd order Butterworth), rectified, and low-pass 
    % filtered at 10 Hz (3rd order Butterworth)
    
    L_EMG_sol = L_EMG_sol(~isnan(L_EMG_sol));
    L_EMG_gastroc = L_EMG_gastroc(~isnan(L_EMG_gastroc));
    R_EMG_sol = R_EMG_sol(~isnan(R_EMG_sol));
    R_EMG_gastroc = R_EMG_gastroc(~isnan(R_EMG_gastroc));
    
    % sol
    if length(R_EMG_sol) < length(L_EMG_sol)
        L_EMG_sol = L_EMG_sol(1:length(R_EMG_sol));
    elseif length(L_EMG_sol) < length(R_EMG_sol)
        R_EMG_sol = R_EMG_sol(1:length(L_EMG_sol));
    end
    % gastroc
    if length(R_EMG_gastroc) < length(L_EMG_gastroc)
        L_EMG_gastroc = L_EMG_gastroc(1:length(R_EMG_gastroc));
    elseif length(L_EMG_gastroc) < length(R_EMG_gastroc)
        R_EMG_gastroc = R_EMG_gastroc(1:length(L_EMG_gastroc));
    end

    %% High pass filter
    % 40 Hz (3rd order Butterworth)
    n = 3;
    C_off = 20;
    Wn = C_off/(Fs/2);
    [b,a] = butter(n,Wn,'high');

    % sol
    L_sol_Filt2 = filtfilt(b,a,L_EMG_sol);
    R_sol_Filt2 = filtfilt(b,a,R_EMG_sol);    
    % gastroc
    L_gastroc_Filt2 = filtfilt(b,a,L_EMG_gastroc);
    R_gastroc_Filt2 = filtfilt(b,a,R_EMG_gastroc);   
    
    %% Rectify
    % sol
    L_sol_rect = abs(L_sol_Filt2);
    R_sol_rect = abs(R_sol_Filt2);
    % gastroc
    L_gastroc_rect = abs(L_gastroc_Filt2);
    R_gastroc_rect = abs(R_gastroc_Filt2);

    %% Low pass filter
    % filtered at 10 Hz (3rd order Butterworth)
    n = 3;
    C_off = 6; %6
    Wn = C_off/(Fs/2);
    [b,a] = butter(n,Wn,'low');

    % sol
    L_sol_Filt3 = filtfilt(b,a,L_sol_rect);
    R_sol_Filt3 = filtfilt(b,a,R_sol_rect);    
    % gastroc
    L_gastroc_Filt3 = filtfilt(b,a,L_gastroc_rect);
    R_gastroc_Filt3 = filtfilt(b,a,R_gastroc_rect); 
    
    %% FILTER
    n = 2;
    Wn = 15/(Fs/2);
    [b,a] = butter(n,Wn,'low');
    
    Mx1 = filtfilt(b,a,Mx1(~isnan(Mx1)));
    Fz1 = filtfilt(b,a,Fz1(~isnan(Fz1)));
    Mx2 = filtfilt(b,a,Mx2(~isnan(Mx2)));
    Fz2 = filtfilt(b,a,Fz2(~isnan(Fz2)));
    
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
        if Fz1(i) < thresLow && Fz1(i) > thresHigh 
            if waitHS1 > 500 && waitTO1 > 300
                if delta1 < 0 && kHS1 == kTO1
                    HS1(kHS1) = i;
                    HS1_log(i) = 1;
                    kHS1 = kHS1+1;
                    waitHS1 = 1;
                end
            end
        end
        if Fz1(i) < thresLow && Fz1(i) > thresHigh 
            if waitTO1 > 500 && waitHS1 > 300
                if delta1 > 0 && kTO1 == kHS1-1
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
        if Fz2(i) < thresLow && Fz2(i) > thresHigh 
            if waitHS2 > 500 && waitTO2 > 300
                if delta2 < 0 && kHS2 == kTO2
                    HS2(kHS2) = i;
                    HS2_log(i) = 1;
                    kHS2 = kHS2+1;
                    waitHS2 = 1;
                end
            end
        end
        if Fz2(i) < thresLow && Fz2(i) > thresHigh 
            if waitTO2 > 500 && waitHS2 > 300
                if delta2 > 0 && kTO2 == kHS2-1
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
        plot(HS1_log.*-100)
        plot(TO1_log.*-100)
        hold on
        plot(Fz1)
        
        figure
        hold on
        plot(HS2_log.*-100)
        plot(TO2_log.*-100)
        hold on
        plot(Fz2)
    end
    
    %% LEFT
    
    if mean(Fz1) < 0
        
        % deal with lengths
        if length(HS1) > length(TO1)
            tLength = length(TO1);
        else
            tLength = length(HS1);
        end
        
        L_ankle_sol_timelock = nan(1500,tLength);
        L_ankle_gastroc_timelock = nan(1500,tLength);
        L_seconds = nan(1500,tLength);

        deltai = abs(HS1(2:tLength)-HS1(1:tLength-1));
        deltai = deltai(~isoutlier(deltai));
        deltai = deltai(deltai<1500);
        upperLim = round(max(deltai))+1;
        lowerLim = round(min(deltai))-1;
        tLii = 1;
        for tLi = 1:tLength-1
            if ~isempty(TO1(TO1 > HS1(tLi) & TO1 < HS1(tLi+1))) % step

                L_sol=L_sol_Filt3(HS1(tLi):HS1(tLi+1));
                L_gastroc=L_gastroc_Filt3(HS1(tLi):HS1(tLi+1));

                if length(L_sol) > lowerLim && length(L_sol) < upperLim && Fz1(HS1(tLi)+100) < mean(Fz1)
                    L_length_sol(tLii) = length(L_sol);
                    L_length_gastroc(tLii) = length(L_gastroc);

                    % norm
                    L_sol_norm = L_sol./LsolpeakNorm;
                    L_gastroc_norm = L_gastroc./LgastrocpeakNorm;
                    L_seconds_init = time(HS1(tLi):HS1(tLi+1)) - time(HS1(tLi));

                    % log
                    if max(L_sol_norm) < 8 && max(L_gastroc_norm) < 8
                        L_ankle_sol_timelock(1:L_length_sol(tLii),tLii) = L_sol_norm;
                        L_ankle_gastroc_timelock(1:L_length_gastroc(tLii),tLii) = L_gastroc_norm;
                        L_seconds(1:length(L_seconds_init),tLii) = L_seconds_init; L_seconds(L_seconds==0) = nan; L_seconds(1,:) = zeros(1,size(L_seconds,2));

                        tLii = tLii+1;
                    end
                end
                clearvars L_sol L_ankle_sol L_gastroc L_ankle_gastroc L_ankleTOi
            end
        end
        for tL = 1:tLii-1
            L_ankle_sol_timelock_resampled(:,tL) = interp1(L_seconds(1:L_length_sol(tL),tL),L_ankle_sol_timelock(1:L_length_sol(tL),tL),linspace(L_seconds(1,tL),L_seconds(L_length_sol(tL),tL),100));
            L_ankle_gastroc_timelock_resampled(:,tL) = interp1(L_seconds(1:L_length_gastroc(tL),tL),L_ankle_gastroc_timelock(1:L_length_gastroc(tL),tL),linspace(L_seconds(1,tL),L_seconds(L_length_gastroc(tL),tL),100));
            L_timelock_resampled(:,tL) = interp1(L_seconds(1:L_length_gastroc(tL),tL),L_seconds(1:L_length_gastroc(tL),tL),linspace(L_seconds(1,tL),L_seconds(L_length_gastroc(tL),tL),100));
            L_gait_cycle(:,tL) = linspace(0,100,100);
        end

        if isnan(L_EMG_sol)
            L_ankle_sol_timelock_resampled = nan(length(L_ankle_gastroc_timelock_resampled),1);
        end
    else
        L_gait_cycle = nan(1,length(R_gait_cycle));
        L_timelock_resampled = nan(1,length(R_timelock_resampled));
        L_ankle_sol_timelock_resampled = nan(1,length(R_ankle_sol_timelock_resampled));
        L_ankle_gastroc_timelock_resampled = nan(1,length(R_ankle_gastroc_timelock_resampled));
        tLii = length(HS2);
        HS1_log = nan;
    end
    
    %% RIGHT
    
    if mean(Fz2) < 0
        
        % deal with lengths
        if length(HS2) > length(TO2)
            tRength = length(TO2);
        else
            tRength = length(HS2);
        end
        
        R_ankle_sol_timelock = nan(1500,tRength);
        R_ankle_gastroc_timelock = nan(1500,tRength);
        R_seconds = nan(1500,tRength);

        deltai = abs(HS2(2:tRength)-HS2(1:tRength-1));
        deltai = deltai(~isoutlier(deltai));
        deltai = deltai(deltai<1500);
        upperLim = round(max(deltai))+1;
        lowerLim = round(min(deltai))-1;
        tRii = 1;
        for tRi = 1:tRength-1
            if ~isempty(TO2(TO2 > HS2(tRi) & TO2 < HS2(tRi+1))) % step

                R_sol=R_sol_Filt3(HS2(tRi):HS2(tRi+1));
                R_gastroc=R_gastroc_Filt3(HS2(tRi):HS2(tRi+1));

                if length(R_sol) > lowerLim && length(R_sol) < upperLim && Fz2(HS2(tRi)+100) < mean(Fz2)
                    R_length_sol(tRii) = length(R_sol);
                    R_length_gastroc(tRii) = length(R_gastroc);

                    % norm
                    R_sol_norm = R_sol./RsolpeakNorm;
                    R_gastroc_norm = R_gastroc./RgastrocpeakNorm;
                    R_seconds_init = time(HS2(tRi):HS2(tRi+1)) - time(HS2(tRi));

                    % log
                    if max(R_sol_norm) < 8 && max(R_gastroc_norm) < 8
                        R_ankle_sol_timelock(1:R_length_sol(tRii),tRii) = R_sol_norm;
                        R_ankle_gastroc_timelock(1:R_length_gastroc(tRii),tRii) = R_gastroc_norm;
                        R_seconds(1:length(R_seconds_init),tRii) = R_seconds_init; R_seconds(R_seconds==0) = nan; R_seconds(1,:) = zeros(1,size(R_seconds,2));

                        tRii = tRii+1;
                    end
                end
                clearvars R_sol R_ankle_sol R_gastroc R_ankle_gastroc R_ankleTOi
            end
        end
        for tR = 1:tRii-1
            R_ankle_sol_timelock_resampled(:,tR) = interp1(R_seconds(1:R_length_sol(tR),tR),R_ankle_sol_timelock(1:R_length_sol(tR),tR),linspace(R_seconds(1,tR),R_seconds(R_length_sol(tR),tR),100));
            R_ankle_gastroc_timelock_resampled(:,tR) = interp1(R_seconds(1:R_length_gastroc(tR),tR),R_ankle_gastroc_timelock(1:R_length_gastroc(tR),tR),linspace(R_seconds(1,tR),R_seconds(R_length_gastroc(tR),tR),100));
            R_timelock_resampled(:,tR) = interp1(R_seconds(1:R_length_gastroc(tR),tR),R_seconds(1:R_length_gastroc(tR),tR),linspace(R_seconds(1,tR),R_seconds(R_length_gastroc(tR),tR),100));
            R_gait_cycle(:,tR) = linspace(0,100,100);
        end

        if isnan(R_EMG_sol)
            R_ankle_sol_timelock_resampled = nan(length(R_ankle_gastroc_timelock_resampled),1);
        end
    else
        R_gait_cycle = nan(1,length(L_gait_cycle));
        R_timelock_resampled = nan(1,length(L_timelock_resampled));
        R_ankle_sol_timelock_resampled = nan(1,length(L_ankle_sol_timelock_resampled));
        R_ankle_gastroc_timelock_resampled = nan(1,length(L_ankle_gastroc_timelock_resampled));
        tRii = tLii;
        HS2_log = nan;
    end
    
    %% Calculate peak
    
    if mean(Fz1) < 0
        Lsolpeak = max(L_ankle_sol_timelock_resampled);
        Lgastrocpeak = max(L_ankle_gastroc_timelock_resampled);
    else
        Lsolpeak = nan;
        Lgastrocpeak = nan;
    end
    
    if mean(Fz2) < 0
        Rsolpeak = max(R_ankle_sol_timelock_resampled);
        Rgastrocpeak = max(R_ankle_gastroc_timelock_resampled);
    else
        Rsolpeak = nan;
        Rgastrocpeak = nan;
    end
    
    %% Calculate total
    if mean(Fz1) < 0
        for m = 1:tLii-1
            L_xi = L_timelock_resampled(:,m); L_x = L_xi(~isnan(L_xi)); 
            Lsol_yi = L_ankle_sol_timelock_resampled(:,m); Lsol_y = Lsol_yi(~isnan(Lsol_yi));
            Lgastroc_yi = L_ankle_gastroc_timelock_resampled(:,m); Lgastroc_y = Lgastroc_yi(~isnan(Lgastroc_yi));
            LtotalsolEMG(m) = trapz(L_x(~isnan(Lsol_y)),Lsol_y(~isnan(Lsol_y)));
            LtotalgastrocEMG(m) = trapz(L_x(~isnan(Lgastroc_y)),Lgastroc_y(~isnan(Lgastroc_y)));
        end
    else
        LtotalsolEMG = nan(1,length(RtotalsolEMG));
        LtotalgastrocEMG = nan(1,length(RtotalgastrocEMG));
    end

    if mean(Fz2) < 0
        for m = 1:tRii-1
            R_xi = R_timelock_resampled(:,m); R_x = R_xi(~isnan(R_xi)); 
            Rsol_yi = R_ankle_sol_timelock_resampled(:,m); Rsol_y = Rsol_yi(~isnan(Rsol_yi));
            Rgastroc_yi = R_ankle_gastroc_timelock_resampled(:,m); Rgastroc_y = Rgastroc_yi(~isnan(Rgastroc_yi));
            RtotalsolEMG(m) = trapz(R_x(~isnan(Rsol_y)),Rsol_y(~isnan(Rsol_y)));
            RtotalgastrocEMG(m) = trapz(R_x(~isnan(Rgastroc_y)),Rgastroc_y(~isnan(Rgastroc_y)));
        end
    else
        RtotalsolEMG = nan(1,length(LtotalsolEMG));
        RtotalgastrocEMG = nan(1,length(LtotalgastrocEMG));
    end
    
    %% plotting
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
            for k = 1:2
                if k == 1
                    figure(jj+10^3)
                elseif k == 2
                    figure(jj+10^4)
                end
                
                if strcmp(tri,'ga1')
                    if k == 1
                        figure(jj+10^3)
                        triname = 'GA1: L/R Sol';
                    elseif k == 2
                        figure(jj+10^4)
                        triname = 'GA1: L/R Gastroc';
                    end
                    subplot(4,8,ii)
                    col = SubjColors(jj,:);
                elseif strcmp(tri,'ga2')
                    if k == 1
                        figure(jj+10^3)
                        triname = 'GA2: L/R Sol';
                    elseif k == 2
                        figure(jj+10^4)
                        triname = 'GA2: L/R Gastroc';
                    end
                    subplot(4,8,ii+2*8)
                    col = SubjColors(jj,:);
                elseif strcmp(tri,'nw1')
                    if k == 1
                        figure(jj+10^3)
                        triname = 'NW1: L/R Sol';
                    elseif k == 2
                        figure(jj+10^4)
                        triname = 'NW1: L/R Gastroc';
                    end
                    subplot(4,8,ii+8)
                    col = [0.776, 0.675, 0.592];
                elseif strcmp(tri,'nw2')
                    if k == 1
                        figure(jj+10^3)
                        triname = 'NW2: L/R Sol';
                    elseif k == 2
                        figure(jj+10^4)
                        triname = 'NW2: L/R Gastroc';
                    end
                    subplot(4,8,ii+3*8)
                    col = [0.776, 0.675, 0.592];
                end

                % PLOT
                if k == 1
                    hold on
                    plot(L_gait_cycle,L_ankle_sol_timelock_resampled,'Color',col)
                    plot([L_gait_cycle(1) L_gait_cycle(end)],[mean(LtotalsolEMG) mean(LtotalsolEMG)],'-k')

                    plot(R_gait_cycle,R_ankle_sol_timelock_resampled,'Color',col)
                    plot([R_gait_cycle(1) R_gait_cycle(end)],[mean(RtotalsolEMG) mean(RtotalsolEMG)],'-k')
                    ylim([-0.5 10])
                    title(strcat('S',string(jj),{' '},triname))
                    xlabel('Time in exo (mins)')
                    ylabel('Sol (norm to NW)')
                elseif k == 2
                    hold on
                    plot(L_gait_cycle,L_ankle_gastroc_timelock_resampled,'Color',col)
                    plot([L_gait_cycle(1) L_gait_cycle(end)],[mean(LtotalgastrocEMG) mean(LtotalgastrocEMG)],'-k')

                    plot(R_gait_cycle,R_ankle_gastroc_timelock_resampled,'Color',col)
                    plot([R_gait_cycle(1) R_gait_cycle(end)],[mean(RtotalgastrocEMG) mean(RtotalgastrocEMG)],'-k')
                    ylim([-0.5 10])
                    title(strcat('S',string(jj),{' '},triname))
                    xlabel('Time in exo (mins)')
                    ylabel('Gastroc (norm to NW)')
                end
            end
        end
    end
    
end