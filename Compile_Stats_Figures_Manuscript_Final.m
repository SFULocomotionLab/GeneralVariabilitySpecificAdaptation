clear; close all

% determine if load raw data or load individual files
loadSavedData = 0; % 0 = load raw data, 1 = load individual files
% determine if saving
saving = 0; % 0 = do not save data, 1 = do save data
% determine if plotting manuscript figures
formatGraphs = 0; % 0 = regular plots, 1 = format plots for manuscript

% collection frequency for all measures
Fs = 500; % Hz
% except EMG
Fs_EMG = 1000; % Hz

% variability/adaptation filtering
Fs_timescale = 1/30;
% for full sensitivity analysis
% Fs_timescale = 1./[10 30 50];

% determine variables for analyses
orig = 0; % 1 = original 4 variables, 0 = other 4 variables

% individual subject colours
% T = static, P = continued optimization, E = re-optimization
SubjColors = [0         0         0;...      % TA: 0	0	0
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

% subjects and their weight
subjT = subjTnaming;
subjTopt = subjToptnaming;
allSubj = {'TA','TB','TC','TD','TE',...
           'PA','PB','PC','PD','PE',...
           'EA','EB','EC','ED','EE'};
allSubjWeights = [58 66 58 64 92,...
                  55 68 75 75 73,...
                  70 93 59 83 73];

% paths for individual participant files and other functions
addpath('/Volumes/GoogleDrive/My Drive/Locomotion Lab/Projects/Cost Optimization - exos/Katies Data/Learning')
addpath('/Volumes/GoogleDrive/My Drive/Locomotion Lab/Projects/Cost Optimization - exos/Katies Data/Analysis - manuscript')

%% initializing variables

% Representative subject
repSubj = 5;

% SF adapt
SFmeanVAL = nan(12,10); SFmeanVALnorm = nan(12,10);
SFmeanVAL_nw = nan(12,10); SFmeanVAL_zt = nan(12,10);

% SF variability
SFstdVAL = nan(12,10,length(Fs_timescale)); SFstdVALnorm = nan(12,10,length(Fs_timescale));
SFstdVAL_nw = nan(12,10,length(Fs_timescale)); SFstdVAL_zt = nan(12,10,length(Fs_timescale));

% Ankle adapt
AnkleVAL = nan(12,10); AnkleVALnorm = nan(12,10);
AnkleVAL_zt = nan(12,10);

% Ankle variability
AnkleStdVAL = nan(12,10,length(Fs_timescale)); AnkleStdVALnorm = nan(12,10,length(Fs_timescale));
AnkleStdVAL_zt = nan(12,10,length(Fs_timescale));

% Soleus adapt
solTotalMeanVAL = nan(12,10); solTotalMeanVALnorm = nan(12,10);
solTotalMeanVAL_nw = nan(12,10); solTotalMeanVAL_zt = nan(12,10);

% Soleus variability
solTotalStdVAL = nan(12,10,length(Fs_timescale)); solTotalStdVALnorm = nan(12,10,length(Fs_timescale));
solTotalStdVAL_nw = nan(12,10,length(Fs_timescale)); solTotalStdVAL_zt = nan(12,10,length(Fs_timescale));

% Gastroc adapt
gastrocTotalMeanVAL = nan(12,10); gastrocTotalMeanVALnorm = nan(12,10);
gastrocTotalMeanVAL_nw = nan(12,10); gastrocTotalMeanVAL_zt = nan(12,10);

% Gastroc variability
gastrocTotalStdVAL = nan(12,10,length(Fs_timescale)); gastrocTotalStdVALnorm = nan(12,10,length(Fs_timescale));
gastrocTotalStdVAL_nw = nan(12,10,length(Fs_timescale)); gastrocTotalStdVAL_zt = nan(12,10,length(Fs_timescale));

% Metabolic cost
metVAL = nan(12,10);
metVAL_zt = nan(12,10);

% trials with poor signal quality (EMG)
if orig == 1
    NoEMG = [4 1 4; 6 3 6; 7 2 7];
else
    NoEMG = [4 1 4; 2 5 6; 5 4 7; 6 4 6; 6 4 7];
end

% trials with poor signal quality (forces for detecting steps)
NoForce = [6 1 2; 6 1 3; 6 1 4; 6 1 5; 6 1 6; 6 1 7;...
           7 2 7; 7 6 7; 6 4 7; 6 5 6];

% time
timepointsVAL = nan(12,10);
timepointsADAPT = nan(21,10);

%% loop through subjects and days (and filter cut-off if sensitivity analysis)

for ss = 1:length(Fs_timescale) % loop through sensitivity analysis
    
    % high pass filter for calculating variability
    n = 3;
    WnHP = Fs_timescale(ss)/(1/2);
    [bHP,aHP] = butter(n,WnHP,'high');
    
    for jj = 1:length(subjT)-5 % loop through all but 5 re-opt subjects
        
        % counting
        kVal_zt = 0; % zero torque
        kVal_nw = 0; % normal walking
        kVal_ga = 0; % generic assistance

        % initialize days and walking time
        cd('/Volumes/GoogleDrive/My Drive/Locomotion Lab/Projects/Cost Optimization - exos/Katies Data/Learning')
        all_days = subjT{jj,2}.day;
        timeExo = 0;

        for ii = 1:length(subjT{jj,2}.day) % loop through days
            
            % navigate to participant's folder on this day
            cd(strcat(string(subjT(jj)),'/',string(all_days(ii))))

            % calculate resting metabolics during quiet standing on this day
            allfiles = subjT{jj,2}.files{ii};
            qsfile = strcmp(allfiles(:,1),'qs');
            load(char(allfiles(qsfile)))
            qsMean = mean(met((end-3*60*Fs):end));
            
            %% additional trials
            % amount of time with exoskeleton assistance (including HILO)
            if ii > 1 % but only second day onwards
                for g = 1:4
                    timeExo = timeExo + 9*2; % 8 CL for 2 mins each + 2 mins GA = 9*2
                end
            end
            
            %% double-reversal trials

            % files
            validFilesCell = subjT{jj,2}.files{ii};
            validFilesMat = cell2mat(validFilesCell(:,2));

            % quiet standing and normal walking validation trials
            qsi = strcmp(validFilesCell(:,1),'qs');
            nw1i = strcmp(validFilesCell(:,1),'nw1');
            nw2i = strcmp(validFilesCell(:,1),'nw2');
            nw = [cell2mat(validFilesCell(nw1i,2)); cell2mat(validFilesCell(nw2i,2))];
            nwSimple = [cell2mat(validFilesCell(nw1i,1)); cell2mat(validFilesCell(nw2i,1))];

            % zero torque and exoskeleton assistance trials
            allDelete = qsi+nw1i+nw2i;
            exoTrialsi = find(allDelete~=1);
            exoTrials = validFilesCell(exoTrialsi,:);
            validFilesSorti = validFilesMat(exoTrialsi,:);
            [validFilesSort,sorti] = sort(fliplr(validFilesSorti(:,end)));
            exoTrialsSort = exoTrials(sorti);

            %% BASELINE normal walking conditions

            for nwtri = 1:size(nw,1)
                if loadSavedData == 0
                    % trial num for checking for bad trials
                    triNum = find(strcmp(validFilesCell(:,2),nw(nwtri,:))~=0);
                    % load data
                    cd(strcat('/Volumes/Learning Da/allData/',string(subjT(jj)),'/',string(all_days(ii))))
                    load(strcat(string(nwSimple(nwtri,:)),'.mat'));
                    % only consider last 3 minutes
                    ts = allData.data(end-3*60*Fs+1:end,end)-allData.data(end-3*60*Fs+1,end);
                    allDatai = allData.data(end-3*60*Fs+1:end,:);

                    %% step frequency
                    
                    % no force data for subj 6 day 1
                    check = ismember(NoForce,[jj ii triNum],'rows');
                    if ismember(1,check) % jj == 6 && ii == 1
                        SF = nan(1,length(ts));
                        SF_HP = nan(1,length(ts));
                        SF_LP = nan(1,length(ts));
                    else
                        if orig == 1
                            [SF,~] = SFcalc(names, Fs, ts, allDatai, jj, ii, nwSimple(nwtri,:));
                        else
                            [SF,~] = SWcalc(names, Fs, ts, allDatai, jj, ii, nwSimple(nwtri,:));
                        end
                        
                        % high-pass filter for calculating variability
                        SF_HP = filtfilt(bHP,aHP,SF(~isnan(SF)));
                        
                        % representative subject
                        if jj == repSubj && nwtri == 1  && ii == 1 && orig == 1
                            allDatarep = allData.data(end-6*60*Fs+1:end,:);
                            [y,t] = SFcalc(names, Fs, ts, allDatai, jj, ii, nan);
                        end
                    end
                    
                    %% muscle activity
                    
                    % bad EMG signal
                    check = ismember(NoEMG,[jj ii triNum],'rows');
                    if ismember(1,check)
                        LsolpeakNorm(nwtri) = nan; RsolpeakNorm(nwtri) = nan;
                        LgastrocpeakNorm(nwtri) = nan; RgastrocpeakNorm(nwtri) = nan;
                    else

                        % determine normalizing constants
                        if orig == 1 % original vars (soleus and gastroc)
                            load(strcat(nw(nwtri,:),'EMGanalysis_sol_gastroc','.mat'));
                            if size(EMGanalysisAll,1) < 180000
                                EMGanalysisAlli = EMGanalysisAll;
                            else
                                EMGanalysisAlli = EMGanalysisAll(end-3*60*Fs_EMG+1:end,:);
                            end
                            [~,~,~, ... % timeseries
                            ~,~,~,... % timeseries
                            Lsolpeak,Rsolpeak,Lgastrocpeak,Rgastrocpeak, ... % peak
                            ~,~,~,~] = EMGCalc_SolGastroc(Fs_EMG, EMGanalysisAlli, 1, 1, 1, 1, jj, ii, nan);
                        
                        else % additional vars (rectus femoris and biceps femoris)
                            load(strcat(nw(nwtri,:),'EMGanalysis_RFBF','.mat'));
                            if size(EMGanalysisAll,1) < 180000
                                EMGanalysisAlli = EMGanalysisAll;
                            else
                                EMGanalysisAlli = EMGanalysisAll(end-3*60*Fs_EMG+1:end,:);
                            end
                            [~,~,~, ... % timeseries
                            ~,~,~,... % timeseries
                            Lsolpeak,Rsolpeak,Lgastrocpeak,Rgastrocpeak, ... % peak
                            ~,~,~,~] = EMGCalc_RFBF(Fs_EMG, EMGanalysisAlli, 1, 1, 1, 1, jj, ii, nan);
                        end
                        
                        % soleus normalizing constants
                        LsolpeakNorm(nwtri) = mean(Lsolpeak);
                        RsolpeakNorm(nwtri) = mean(Rsolpeak);
                        % gastroc normalizing constants
                        LgastrocpeakNorm(nwtri) = mean(Lgastrocpeak);
                        RgastrocpeakNorm(nwtri) = mean(Rgastrocpeak);

                        % determine nw avg and std values, normalized to previously calculated average peak
                        if orig == 1
                            [~,~,~, ... % timeseries
                            ~,~,~,... % timeseries
                            ~,~,~,~, ... % peak
                            RtotalsolEMG,LtotalsolEMG,RtotalgastrocEMG,LtotalgastrocEMG] = EMGCalc_SolGastroc(Fs_EMG, EMGanalysisAlli, LsolpeakNorm(nwtri), RsolpeakNorm(nwtri), LgastrocpeakNorm(nwtri), RgastrocpeakNorm(nwtri), jj, ii, nwSimple(nwtri,:));
                        else
                            [~,~,~, ... % timeseries
                            ~,~,~,... % timeseries
                            ~,~,~,~, ... % peak
                            RtotalsolEMG,LtotalsolEMG,RtotalgastrocEMG,LtotalgastrocEMG] = EMGCalc_RFBF(Fs_EMG, EMGanalysisAlli, LsolpeakNorm(nwtri), RsolpeakNorm(nwtri), LgastrocpeakNorm(nwtri), RgastrocpeakNorm(nwtri), jj, ii, nwSimple(nwtri,:));
                        end
                        
                        % filter
                        % if no right foot data
                        if ~isnan(RtotalsolEMG)
                            if length(RtotalsolEMG) > 9
                                RtotalsolEMG_HP = filtfilt(bHP,aHP,RtotalsolEMG);
                                RtotalgastrocEMG_HP = filtfilt(bHP,aHP,RtotalgastrocEMG);
                            end
                        end
                        % if no left foot data
                        if ~isnan(LtotalsolEMG)
                            if length(LtotalsolEMG) > 9
                                LtotalsolEMG_HP = filtfilt(bHP,aHP,LtotalsolEMG);
                                LtotalgastrocEMG_HP = filtfilt(bHP,aHP,LtotalgastrocEMG);
                            end
                        end

                        % average
                        solTotalMEANVAL_nwi = mean([mean(LtotalsolEMG) mean(RtotalsolEMG)]);
                        gastrocTotalMEANVAL_nwi = mean([mean(LtotalgastrocEMG) mean(RtotalgastrocEMG)]);
                        % variability
                        solTotalSTDVAL_nwi = mean([std(LtotalsolEMG_HP) std(RtotalsolEMG_HP)]);
                        gastrocTotalSTDVAL_nwi = mean([std(LtotalgastrocEMG_HP) std(RtotalgastrocEMG_HP)]);

                        % representative subject
                        if jj == repSubj && nwtri == 1  && ii == 1 && orig == 1
                            clearvars t y yLP
                            y = [];
                            
                            [RtEMGplot_nw,RSolplot_nw,RGastrocplot_nw, ... % timeseries
                            LtEMGplot_nw,LSolplot_nw,LGastrocplot_nw,... % timeseries
                            ~,~,~,~, ... % peak
                            ~,~,~,~] = EMGCalc_SolGastroc(Fs_EMG, EMGanalysisAlli, LsolpeakNorm(nwtri), RsolpeakNorm(nwtri), LgastrocpeakNorm(nwtri), RgastrocpeakNorm(nwtri), jj, ii, 'n');
                            
                            t = LtEMGplot_nw(:,1);
                            y{1} = LSolplot_nw;
                            y{2} = LGastrocplot_nw;
                        end
                            
                        %% saving
                        
                        % SF
                        cd(strcat('/Volumes/GoogleDrive/My Drive/Locomotion Lab/Projects/Cost Optimization - exos/Katies Data/Learning/Saving'))
                        if saving == 1 && jj == repSubj && nwtri == 1  && ii == 1 && orig == 1
                            save(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-nw1SF.mat'),'SF_HP','SF','t','y')
                        elseif saving == 1
                            save(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',nw(nwtri,:),'SF.mat'),'SF_HP','SF')
                        end
                        
                        % EMG
                        if saving == 1 && jj == repSubj && nwtri == 1 && ii == 1 && orig == 1
                            save(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-nw1EMG.mat'),'LsolpeakNorm','RsolpeakNorm','LgastrocpeakNorm','RgastrocpeakNorm',...
                                        'solTotalMEANVAL_nwi','gastrocTotalMEANVAL_nwi','solTotalSTDVAL_nwi','gastrocTotalSTDVAL_nwi', ...
                                        't','y');
                        elseif saving == 1
                            save(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',nw(nwtri,:),'EMG.mat'),'LsolpeakNorm','RsolpeakNorm','LgastrocpeakNorm','RgastrocpeakNorm',...
                                        'solTotalMEANVAL_nwi','gastrocTotalMEANVAL_nwi','solTotalSTDVAL_nwi','gastrocTotalSTDVAL_nwi')
                        end
                    end
                else
                    cd(strcat('/Volumes/GoogleDrive/My Drive/Locomotion Lab/Projects/Cost Optimization - exos/Katies Data/Learning/Saving'))
                    load(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',nw(nwtri,:),'EMG.mat'))
                    load(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',nw(nwtri,:),'SF.mat'))
                end

                kVal_nw = kVal_nw+1;

                % soleus
                solTotalMeanVAL_nw(kVal_nw,jj) = solTotalMEANVAL_nwi;
                solTotalStdVAL_nw(kVal_nw,jj,ss) = solTotalSTDVAL_nwi;

                % gastroc
                gastrocTotalMeanVAL_nw(kVal_nw,jj) = gastrocTotalMEANVAL_nwi;
                gastrocTotalStdVAL_nw(kVal_nw,jj,ss) = gastrocTotalSTDVAL_nwi;

                % step frequency
                SFmeanVAL_nw(kVal_nw,jj) = nanmean(SF);
                SFstdVAL_nw(kVal_nw,jj,ss) = nanstd(SF_HP);
                
            end

            %% ASSISTANCE conditions

            for tri = 1:length(exoTrialsSort)
                exoTrialsClassifyi = cell2mat(exoTrialsSort(tri,:));
                exoTrialsClassify = exoTrialsClassifyi(1);
                disp(strcat(string(subjT(jj)),'--',string(ss),'--',string(all_days(ii)),'--',exoTrialsClassifyi))

                if loadSavedData == 0
                    % load data
                    cd(strcat('/Volumes/Learning Da/allData/',string(subjT(jj)),'/',string(all_days(ii))))
                    load(strcat(string(exoTrialsClassifyi),'.mat'));
                    % only consider last 3 minutes
                    ts = allData.data(end-3*60*Fs+1:end,end)-allData.data(end-3*60*Fs+1,end);
                    allDatai = allData.data(end-3*60*Fs+1:end,:);

                    if exoTrialsClassifyi(1) == 'g' || exoTrialsClassifyi(1) == 'z'
                        % trials without force data
                        triNum = find(strcmp(validFilesCell(:,1),exoTrialsSort(tri,1))~=0);
                        check = ismember(NoForce,[jj ii triNum],'rows');
                        if ismember(1,check) % bad trials
                            SF_HP = nan(1,length(ts));
                            SF_LP = nan(1,length(ts));
                            SFtimeseries_x = nan; SFtimeseries_y = nan;

                            L_avgAnkleRange = nan;
                            R_avgAnkleRange = nan;

                            L_AvgAnkleStd = nan;
                            R_AvgAnkleStd = nan;
                            AnkleAngletimeseries_y = nan;

                            netMetPow = nan;
                        else % good trials
                            
                            %% step frequency
                            if orig == 1
                                [SF,~] = SFcalc(names, Fs, ts, allDatai, jj, ii, exoTrialsClassifyi);
                            else
                                [SF,~] = SWcalc(names, Fs, ts, allDatai, jj, ii, exoTrialsClassifyi);
                            end
                            
                            % high-pass filter for calculating variability
                            SF_HP = filtfilt(bHP,aHP,SF(~isnan(SF)));
                            
                            %% ankle angle
                            if orig == 1
                                [L_AnkleRange, R_AnkleRange, ~, ~, ~, ~] = AnkleAngleRange(names, Fs, allDatai, jj, ii, exoTrialsClassifyi, tauL(end-3*60*Fs+1:end,:));
                            else
                                [L_AnkleRange, R_AnkleRange, ~, ~, ~, ~] = AnkleAnglePeak(names, Fs, allDatai, jj, ii, exoTrialsClassifyi, tauL(end-3*60*Fs+1:end,:));
                            end
                            
                            % range
                            L_avgAnkleRange = nanmean(L_AnkleRange);
                            R_avgAnkleRange = nanmean(R_AnkleRange);

                            % variability
                            L_AnkleRange_HP = filtfilt(bHP,aHP,L_AnkleRange(~isnan(L_AnkleRange)));
                            L_AvgAnkleStd = nanstd(L_AnkleRange_HP);
                            R_AnkleRange_HP = filtfilt(bHP,aHP,R_AnkleRange(~isnan(R_AnkleRange)));
                            R_AvgAnkleStd = nanstd(R_AnkleRange_HP);

                            %% metabolics
                            meti = (met - qsMean);
                            netMetPow = mean(meti((end-3*60*Fs+1):end));
                            
                            %% representative subject
                            if jj == repSubj && orig == 1
                                clearvars t y
                                
                                % SF
                                [SFploti,tSFploti] = SFcalc(names, Fs, ts, allDatai, jj, ii, nan);
                                
                                % ankle
                                [~, ~, LAnkleploti, RAnkleploti, LtAnkleploti, RtAnkleploti] = AnkleAngleRange(names, Fs, allDatai, jj, ii, nan);
                                
                                cd(strcat('/Volumes/GoogleDrive/My Drive/Locomotion Lab/Projects/Cost Optimization - exos/Katies Data/Learning/Saving'))
                                % first day
                                if strcmp(exoTrialsClassifyi,'ga1') && ii == 1
                                    t = tSFploti(~isnan(tSFploti)); t = t(t~=0); t = t-t(1);
                                    y = SFploti;
                                    save(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'SF.mat'),'SF_HP','SF','t','y')
                                    
                                    t = LtAnkleploti(:,1);
                                    y = LAnkleploti;
                                    save(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'AnkleAngleandRange.mat'),'L_avgAnkleRange','R_avgAnkleRange','t','y')
                                % last day
                                elseif strcmp(exoTrialsClassifyi,'ga1') && ii == length(subjT{jj,2}.day)
                                    t = tSFploti(~isnan(tSFploti)); t = t(t~=0); t = t-t(1);
                                    y = SFploti;
                                    save(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'SF.mat'),'SF_HP','SF','t','y')
                                    
                                    t = LtAnkleploti(:,1);
                                    y = LAnkleploti;
                                    save(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'AnkleAngleandRange.mat'),'L_avgAnkleRange','R_avgAnkleRange','t','y')
                                % zero torque
                                elseif strcmp(exoTrialsClassifyi,'zt1') && ii == 1
                                    t = LtAnkleploti(:,1);
                                    y = LAnkleploti;
                                    save(strcat(string(subjT(jj)),'--',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'AnkleAngleandRange.mat'),'L_avgAnkleRange','R_avgAnkleRange','t','y')
                                end
                            end
                        end
                    
                        %% muscle activity
                        triNum = find(strcmp(validFilesCell(:,1),exoTrialsSort(tri,1))~=0);
                        check = ismember(NoEMG,[jj ii triNum],'rows');
                        if ismember(1,check)
                            LtotalMEANsolEMG = nan; RtotalMEANsolEMG = nan;
                            LtotalMEANgastrocEMG = nan; RtotalMEANgastrocEMG = nan;
                            LtotalSTDsolEMG = nan; RtotalSTDsolEMG = nan;
                            LtotalSTDgastrocEMG = nan; RtotalSTDgastrocEMG = nan;
                        else
                            cd(strcat('/Volumes/Learning Da/allData/',string(subjT(jj)),'/',string(all_days(ii))))
                            exoTrialsSortEMG = exoTrials(sorti,:);
                            exoTrialsClassifyEMG = cell2mat(exoTrialsSortEMG(tri,2));

                            % normalizing constants (average peak NW for same day)
                            LsolpeakNorm = nanmean(LsolpeakNorm);
                            RsolpeakNorm = nanmean(RsolpeakNorm);
                            LgastrocpeakNorm = nanmean(LgastrocpeakNorm);
                            RgastrocpeakNorm = nanmean(RgastrocpeakNorm);

                            if orig == 1
                                load(strcat(exoTrialsClassifyEMG,'EMGanalysis_sol_gastroc','.mat'));
                                if length(EMGanalysisAll) < 180000
                                    EMGanalysisAlli = EMGanalysisAll;
                                else
                                    EMGanalysisAlli = EMGanalysisAll(end-3*60*Fs_EMG+1:end,:);
                                end
                                [~,~,~, ... % timeseries
                                ~,~,~,... % timeseries
                                ~,~,~,~, ... % peak
                                RtotalsolEMG,LtotalsolEMG,RtotalgastrocEMG,LtotalgastrocEMG] = EMGCalc_SolGastroc(Fs_EMG, EMGanalysisAlli, LsolpeakNorm, RsolpeakNorm, LgastrocpeakNorm, RgastrocpeakNorm, jj, ii, exoTrialsClassifyi);
                            else
                                load(strcat(exoTrialsClassifyEMG,'EMGanalysis_RFBF','.mat'));
                                if length(EMGanalysisAll) < 180000
                                    EMGanalysisAlli = EMGanalysisAll;
                                else
                                    EMGanalysisAlli = EMGanalysisAll(end-3*60*Fs_EMG+1:end,:);
                                end
                                [~,~,~, ... % timeseries
                                ~,~,~,... % timeseries
                                ~,~,~,~, ... % peak
                                RtotalsolEMG,LtotalsolEMG,RtotalgastrocEMG,LtotalgastrocEMG] = EMGCalc_RFBF(Fs_EMG, EMGanalysisAlli, LsolpeakNorm, RsolpeakNorm, LgastrocpeakNorm, RgastrocpeakNorm, jj, ii, exoTrialsClassifyi);
                            end

                            % filter
                            % if no right foot data
                            if ~isnan(RtotalsolEMG)
                                if length(RtotalsolEMG) > 9
                                    RtotalsolEMG_HP = filtfilt(bHP,aHP,RtotalsolEMG);
                                    RtotalgastrocEMG_HP = filtfilt(bHP,aHP,RtotalgastrocEMG);
                                end
                            end
                            % if no left foot data
                            if ~isnan(LtotalsolEMG)
                                if length(LtotalsolEMG) > 9
                                    LtotalsolEMG_HP = filtfilt(bHP,aHP,LtotalsolEMG);
                                    LtotalgastrocEMG_HP = filtfilt(bHP,aHP,LtotalgastrocEMG);
                                end
                            end

                            % soleus average
                            LtotalMEANsolEMG = nanmean(LtotalsolEMG);
                            RtotalMEANsolEMG = nanmean(RtotalsolEMG);
                            % soleus variability
                            LtotalSTDsolEMG = nanstd(LtotalsolEMG_HP);
                            RtotalSTDsolEMG = nanstd(RtotalsolEMG_HP);
                            % gastroc average
                            LtotalMEANgastrocEMG = nanmean(LtotalgastrocEMG);
                            RtotalMEANgastrocEMG = nanmean(RtotalgastrocEMG);
                            % gastroc variability
                            LtotalSTDgastrocEMG = nanstd(LtotalgastrocEMG_HP);
                            RtotalSTDgastrocEMG = nanstd(RtotalgastrocEMG_HP);

                            % representative subject
                            if jj == repSubj && orig == 1
                                clearvars t y yLP

                                [RtEMGploti,RSolploti,RGastrocploti, ... % timeseries
                                LtEMGploti,LSolploti,LGastrocploti,... % timeseries
                                ~,~,~,~, ... % peak
                                ~,~,~,~] = EMGCalc_SolGastroc(Fs_EMG, EMGanalysisAlli, LsolpeakNorm, RsolpeakNorm, LgastrocpeakNorm, RgastrocpeakNorm, jj, ii, 'n');

                                if strcmp(exoTrialsClassifyi,'ga1') && ii == 1
                                    t = LtEMGploti(:,1);
                                    y{1} = LSolploti;
                                    y{2} = LGastrocploti;
                                % last day
                                elseif strcmp(exoTrialsClassifyi,'ga1') && ii == length(subjT{jj,2}.day)
                                    t = LtEMGploti(:,1);
                                    y{1} = LSolploti;
                                    y{2} = LGastrocploti;
                                end
                            end
                        end
                        cd(strcat('/Volumes/GoogleDrive/My Drive/Locomotion Lab/Projects/Cost Optimization - exos/Katies Data/Learning/Saving'))
                        if saving == 1 && jj == repSubj && strcmp(exoTrialsClassifyi,'ga1') && ii == 1 && orig == 1
                            save(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'EMG.mat'),'LtotalMEANsolEMG','RtotalMEANsolEMG','LtotalMEANgastrocEMG','RtotalMEANgastrocEMG',...
                                'LtotalSTDsolEMG','RtotalSTDsolEMG','LtotalSTDgastrocEMG','RtotalSTDgastrocEMG', ...
                                't','y');
                        elseif saving == 1 && jj == repSubj && strcmp(exoTrialsClassifyi,'ga1') && ii == length(subjT{jj,2}.day) && orig == 1
                            save(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'EMG.mat'),'LtotalMEANsolEMG','RtotalMEANsolEMG','LtotalMEANgastrocEMG','RtotalMEANgastrocEMG',...
                                'LtotalSTDsolEMG','RtotalSTDsolEMG','LtotalSTDgastrocEMG','RtotalSTDgastrocEMG', ...
                                't','y');
                        elseif saving == 1
                            save(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'EMG.mat'),'LtotalMEANsolEMG','RtotalMEANsolEMG','LtotalMEANgastrocEMG','RtotalMEANgastrocEMG',...
                                'LtotalSTDsolEMG','RtotalSTDsolEMG','LtotalSTDgastrocEMG','RtotalSTDgastrocEMG');
                        end
                    end
                end

                if loadSavedData == 0
                    cd(strcat('/Volumes/GoogleDrive/My Drive/Locomotion Lab/Projects/Cost Optimization - exos/Katies Data/Learning/Saving'))
                    if saving == 1 && jj == repSubj && strcmp(exoTrialsClassifyi,'ga1') && ii == 1 && orig == 1
                    elseif saving == 1 && jj == repSubj && strcmp(exoTrialsClassifyi,'ga1') && ii == length(subjT{jj,2}.day) && orig == 1
                    elseif saving == 1 && jj == repSubj && strcmp(exoTrialsClassifyi,'zt1') && ii == 1 && orig == 1
                    elseif saving == 1
                        save(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'.mat'),'SF_HP','SF')
                        save(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'AnkleAngleandRange.mat'),'L_avgAnkleRange','R_avgAnkleRange')
                        save(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'AnkleAngleStd.mat'),'L_AvgAnkleStd','R_AvgAnkleStd')
                        save(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'Metabolics.mat'),'netMetPow')
                    end
                elseif exoTrialsClassifyi(1) == 'g' || exoTrialsClassifyi(1) == 'z'
                    cd(strcat('/Volumes/GoogleDrive/My Drive/Locomotion Lab/Projects/Cost Optimization - exos/Katies Data/Learning/Saving'))
                    load(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'.mat'))
                    if jj == repSubj && orig == 1 && strcmp(exoTrialsClassifyi,'zt1') && ii == 1
                        load(strcat(string(subjT(jj)),'--',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'AnkleAngleandRange.mat'))
                    else
                        load(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'AnkleAngleandRange.mat'))
                    end
                    load(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'AnkleAngleStd.mat'))
                    load(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'Metabolics.mat'))
                    load(strcat(string(subjT(jj)),'-',string(all_days(ii)),'-VAL-',exoTrialsClassifyi,'EMG.mat'))
                end

                % generic assistance
                if exoTrialsClassify == 'g'

                    kVal_ga = kVal_ga + 1;
                    
                    % determine time walking with exoskeleton assistance
                    timeExo = timeExo + 6;
                    timepointsVAL(kVal_ga,jj) = timeExo;

                    %% step frequency
                    % SF mean
                    SFmeanVAL(kVal_ga,jj) = nanmean(SF);

                    % SF variability
                    SFstdVAL(kVal_ga,jj,ss) = nanstd(SF_HP);
                    
                    %% ankle
                    % mean
                    AnkleVAL(kVal_ga,jj) = nanmean([L_avgAnkleRange R_avgAnkleRange]);

                    % variability
                    AnkleStdVAL(kVal_ga,jj,ss) = nanmean([L_AvgAnkleStd R_AvgAnkleStd]);

                    %% metabolics
                    metVAL(kVal_ga,jj) = netMetPow./allSubjWeights(jj);

                    %% muscle activity
                    triNum = find(strcmp(validFilesCell(:,1),exoTrialsSort(tri,1))~=0);
                    check = ismember(NoEMG,[jj ii triNum],'rows');
                    if ismember(1,check)
                        solTotalMeanVAL(kVal_ga,jj) = nan;
                        solTotalStdVAL(kVal_ga,jj,ss) = nan;
                        gastrocTotalMeanVAL(kVal_ga,jj) = nan;
                        gastrocTotalStdVAL(kVal_ga,jj,ss) = nan;
                    elseif jj == 3 && ii == 1 && triNum == 6
                        % Soleus
                        solTotalMeanVAL(kVal_ga,jj) = nanmean([LtotalMEANsolEMG RtotalMEANsolEMG]);
                        solTotalStdVAL(kVal_ga,jj,ss) = nanmean([LtotalSTDsolEMG RtotalSTDsolEMG]);

                        % Gastroc
                        gastrocTotalMeanVAL(kVal_ga,jj) = RtotalMEANgastrocEMG;
                        gastrocTotalStdVAL(kVal_ga,jj,ss) = RtotalSTDgastrocEMG;
                    else
                        % Soleus
                        solTotalMeanVAL(kVal_ga,jj) = nanmean([LtotalMEANsolEMG RtotalMEANsolEMG]);
                        solTotalStdVAL(kVal_ga,jj,ss) = nanmean([LtotalSTDsolEMG RtotalSTDsolEMG]);

                        % Gastroc
                        gastrocTotalMeanVAL(kVal_ga,jj) = nanmean([LtotalMEANgastrocEMG RtotalMEANgastrocEMG]);
                        gastrocTotalStdVAL(kVal_ga,jj,ss) = nanmean([LtotalSTDgastrocEMG RtotalSTDgastrocEMG]);
                    end

                % optimized assistance
                elseif exoTrialsClassify == 'o'
                    timeExo = timeExo + 6;

                % zero torque
                elseif exoTrialsClassify == 'z'
                    kVal_zt = kVal_zt + 1;
                    
                    % SF mean
                    SFmeanVAL_zt(kVal_zt,jj) = nanmean(SF);
                    % SF variability
                    SFstdVAL_zt(kVal_zt,jj,ss) = nanstd(SF_HP);
                    
                    % Ankle angle
                    AnkleVAL_zt(kVal_zt,jj) = nanmean([L_avgAnkleRange R_avgAnkleRange]);
                    % Ankle angle variability
                    AnkleStdVAL_zt(kVal_zt,jj,ss) = nanmean([L_AvgAnkleStd R_AvgAnkleStd]);

                    % Soleus
                    solTotalMeanVAL_zt(kVal_zt,jj) = nanmean([LtotalMEANsolEMG RtotalMEANsolEMG]);
                    solTotalStdVAL_zt(kVal_zt,jj,ss) = nanmean([LtotalSTDsolEMG RtotalSTDsolEMG]);
                    % Gastroc
                    gastrocTotalMeanVAL_zt(kVal_zt,jj) = nanmean([LtotalMEANgastrocEMG RtotalMEANgastrocEMG]);
                    gastrocTotalStdVAL_zt(kVal_zt,jj,ss) = nanmean([LtotalSTDgastrocEMG RtotalSTDgastrocEMG]);
                    
                    % metabolics
                    metVAL_zt(kVal_zt,jj) = netMetPow./allSubjWeights(jj);
                end
            end
            
            cd('../')
        end
        
        %% first day index
        if jj == 6
            ind = 3; % subj 6 missing data on day 1
        else
            ind = 1;
        end

        %% Time

        % ga validation
        timepointsVALi = timepointsVAL(~isnan(timepointsVAL(:,jj)),jj);

        %% SF

        % adaptation
        SFmeanVALnorm(:,jj) = SFmeanVAL(:,jj)./nanmean(SFmeanVAL_nw(ind:ind+1,jj));
        yAdapt{1,1} = 'SF'; yAdapt{1,2} = SFmeanVALnorm;

        % variability
        SFstdVALnorm(:,jj,ss) = SFstdVAL(:,jj,ss)./nanmean(SFstdVAL_nw(ind:ind+1,jj,ss));
        yStd{1,1} = 'SF'; yStd{1,2} = SFstdVALnorm;

        %% ankle angle range

        % adaptation
        AnkleVALnorm(:,jj) = AnkleVAL(:,jj)./nanmean(AnkleVAL_zt(ind:ind+1,jj));
        yAdapt{2,1} = 'Ankle'; yAdapt{2,2} = AnkleVALnorm;

        % variability
        AnkleStdVALnorm(:,jj,ss) = AnkleStdVAL(:,jj,ss)./nanmean(AnkleStdVAL_zt(ind:ind+1,jj,ss));
        yStd{2,1} = 'Ankle'; yStd{2,2} = AnkleStdVALnorm;

        %% Soleus

        % adaptation
        solTotalMeanVALnorm(:,jj) = solTotalMeanVAL(:,jj)./nanmean(solTotalMeanVAL_nw(ind:ind+1,jj));
        yAdapt{3,1} = 'Soleus'; yAdapt{3,2} = solTotalMeanVALnorm;
        
        % variability
        solTotalStdVALnorm(:,jj,ss) = solTotalStdVAL(:,jj,ss)./nanmean(solTotalStdVAL_nw(ind:ind+1,jj,ss));
        yStd{3,1} = 'Soleus'; yStd{3,2} = solTotalStdVALnorm;

        %% Gastroc

        % adaptation
        gastrocTotalMeanVALnorm(:,jj) = gastrocTotalMeanVAL(:,jj)./nanmean(gastrocTotalMeanVAL_nw(ind:ind+1,jj));
        yAdapt{4,1} = 'Gastroc'; yAdapt{4,2} = gastrocTotalMeanVALnorm;
        
        % variability
        gastrocTotalStdVALnorm(:,jj,ss) = gastrocTotalStdVAL(:,jj,ss)./nanmean(gastrocTotalStdVAL_nw(ind:ind+1,jj,ss));
        yStd{4,1} = 'Gastroc'; yStd{4,2} = gastrocTotalStdVALnorm;

        %% plotting each participant's adaptation

        if 0 %loadSavedData == 0
            titleAdapt = {'SF';'Ankle';'Soleus';'Medial Gastroc'};
            
            varsMeanVAL = {'SFmeanVAL';'AnkleVAL';...
                'solTotalMeanVAL';'gastrocTotalMeanVAL'};
            varsMeanVALnorm = {'SFmeanVALnorm';'AnkleVALnorm';...
                'solTotalMeanVALnorm';'gastrocTotalMeanVALnorm'};
            ylabelAdapt = {'SF Adapt (bpm)';'Ankle Adapt (deg)';...
                'Soleus Adapt (norm to NW)';'Gastroc Adapt (norm to NW)'};
            ylabelAdaptnorm = {'SF Adapt (norm to NW)';'Ankle Adapt (norm to ZT)';...
                'Soleus Adapt (norm to NW)';'Gastroc Adapt (norm to NW)'};
            ylimAdapt = [80 140; 0 60; 0 1; 0 1];
            ylimAdaptnorm = [0.8 1.2; 0 4; 0.6 1.6; 0.6 1.6; 0 5];
            
            for dimAdapt = 1:8
                if dimAdapt < 5
                    figureMean(dimAdapt) = figure(jj+10^dimAdapt);
                    axAdapt(jj,dimAdapt) = subplot(4,8,[7+2*8 7+3*8]);
                    yVAL = eval(string(varsMeanVAL(dimAdapt)));
                else
                    figureMean(dimAdapt) = figure(jj+10^(dimAdapt-4));
                    axAdapt(jj,dimAdapt) = subplot(4,8,[8+2*8 8+3*8]);
                    yVAL = eval(string(varsMeanVALnorm(dimAdapt-4)));
                end
                hold on

                % fitting exp
                modelFun = @(b,daySave)b(1)*exp(-daySave/b(2)) + b(3);

                % plotting data
                plot(timepointsVAL(~isnan(yVAL(:,jj)),jj),yVAL(~isnan(yVAL(:,jj)),jj),'.','MarkerSize',15,'Color',SubjColors(jj,:));

                if dimAdapt < 5
                    ylabel(ylabelAdapt(dimAdapt))
                    ylim(ylimAdapt(dimAdapt,:))
                else
                    if dimAdapt == 2 || dimAdapt == 6
                        TCMean(dimAdapt-4,jj) = nan;
                    else
                        try
                        % fit
                            start = [1 100 1];
                            y = yVAL(~isnan(yVAL(:,jj)),jj);
                            t = timepointsVAL(~isnan(yVAL(:,jj)),jj);
                            xunc = fitnlm(t,y,modelFun,start);
                            tt = linspace(min(t), max(t),100)';
                            line(tt,predict(xunc,tt),'LineWidth', 1.5, 'color', SubjColors(jj,:))
                            TCMean(dimAdapt-4,jj) = table2array(xunc.Coefficients(2,1));
                        catch
                            TCMean(dimAdapt-4,jj) = nan;
                        end
                    end
                    
                    % plot settings
                    title(strcat('TC=',string(round(TCMean(dimAdapt-4,jj)))))
                    plot([0 600],[1 1],'--k')
                    ylabel(ylabelAdaptnorm(dimAdapt-4))
                    ylim(ylimAdaptnorm(dimAdapt-4,:))
                end
                
                % plot settings
                xlim([0 600])
                xticks(0:200:600)
                xlabel('Time in exo (mins)')
            end

            %% plotting each participant's exploration

            titleExplore = {'SF Var';'Ankle Var';'Soleus Var';'Medial Gastroc Var'};
            varsStdVAL = {'SFstdVAL';'AnkleStdVAL';...
                'solTotalStdVAL';'gastrocTotalStdVAL'};
            varsStdVALnorm = {'SFstdVALnorm';'AnkleStdVALnorm';...
                'solTotalStdVALnorm';'gastrocTotalStdVALnorm'};
            ylabelExplore = {'SF Explore (bpm)';'Ankle Explore (deg)';...
                'Soleus Explore (norm to NW)';'Gastroc Explore (norm to NW)'};
            ylabelExplorenorm = {'SF Explore (norm to NW)';'Ankle Explore (norm to ZT)';...
                'Soleus Explore (norm to NW)';'Gastroc Explore (norm to NW)'};
            ylimExplore = [0 6; 0 6; 0 0.08; 0 0.08];
            ylimExplorenorm = [0 3; -2.5 8; 0 3; 0 3];
            
            for dimExplore = 1:8
                if dimExplore < 5
                    figureVar(dimExplore) = figure(jj+10^dimExplore);
                    axExplore(jj,dimExplore) = subplot(4,8,[7 7+8]);
                    yVAL = eval(string(varsStdVAL(dimExplore)));
                else
                    figureVar(dimExplore) = figure(jj+10^(dimExplore-4));
                    axExplore(jj,dimExplore) = subplot(4,8,[8 8+8]);
                    yVAL = eval(string(varsStdVALnorm(dimExplore-4)));
                end
                hold on

                % fitting exp
                modelFun = @(b,daySave)b(1)*exp(-daySave/b(2)) + b(3);

                % plotting data
                plot(timepointsVAL(~isnan(yVAL(:,jj)),jj),yVAL(~isnan(yVAL(:,jj)),jj),'.','MarkerSize',15,'Color',SubjColors(jj,:));

                if dimExplore < 5
                    ylabel(ylabelExplore(dimExplore))
                    ylim(ylimExplore(dimExplore,:))
                else
                    % fit
                    try
                        start = [10 400 1];
                        y = yVAL(~isnan(yVAL(:,jj)),jj);
                        t = timepointsVAL(~isnan(yVAL(:,jj)),jj);
                        xunc = fitnlm(t,y,modelFun,start);
                        tt = linspace(min(t), max(t),100)';
                        line(tt,predict(xunc,tt),'LineWidth', 1.5, 'color', SubjColors(jj,:))
                        TCVar(dimExplore-4,jj) = table2array(xunc.Coefficients(2,1));
                    catch
                        TCVar(dimExplore-4,jj) = nan;
                    end
                    
                    % plot settings
                    title(strcat('TC=',string(round(TCVar(dimExplore-4,jj)))))
                    plot([0 600],[1 1],'--k')
                    ylabel(ylabelExplorenorm(dimExplore-4))
                    ylim(ylimExplorenorm(dimExplore-4,:))
                    if saving == 1 && orig == 1
                        cd('/Volumes/GoogleDrive/My Drive/Locomotion Lab/Projects/Cost Optimization - exos/Katies Data/Analysis - manuscript/Check Figures - main 4 dims')
                        figureName = strcat(titleAdapt(dimExplore-4),'-',allSubj(jj),'.png');
                        saveas(figureVar(dimExplore-4),string(figureName))
                        save(strcat(string(subjT(jj)),'-TC.mat'),'TCVar','TCMean')
                    elseif saving == 0 && orig == 0
                        cd('/Volumes/GoogleDrive/My Drive/Locomotion Lab/Projects/Cost Optimization - exos/Katies Data/Analysis - manuscript/Check Figures - other 4 dims')
                        figureName = strcat(titleAdapt(dimExplore-4),'-',allSubj(jj),'.png');
                        saveas(figureVar(dimExplore-4),string(figureName))
                        save(strcat(string(subjT(jj)),'-TC.mat'),'TCVar','TCMean')
                    end
                end
                
                % plot settings
                xlim([0 600])
                xticks(0:200:600)
                xlabel('Time in exo (mins)')
                set(gcf,'units','inches','position',[0,0,15,10])

            end

            cd('/Volumes/GoogleDrive/My Drive/Locomotion Lab/Projects/Cost Optimization - exos/Katies Data/Learning')
            
            %% plotting each participant's metabolics

            figure(9)
            subplot(2,5,jj);
            hold on

            % data
            plot(timepointsVAL(~isnan(metVAL(:,jj)),jj),metVAL(~isnan(metVAL(:,jj)),jj),'.','MarkerSize',15,'Color',SubjColors(jj,:));
            y = metVAL(~isnan(metVAL(:,jj)),jj);
            t = timepointsVAL(~isnan(metVAL(:,jj)),jj);

            % fit
            modelFun = @(b,daySave)b(1)*exp(-daySave/b(2)) + b(3);
            start = [1 550 1];            
            xunc = fitnlm(t,y,modelFun,start);
            tt = linspace(min(t), max(t),100)';
            line(tt,predict(xunc,tt),'LineWidth', 1.5, 'color', SubjColors(jj,:))
            TCMet(jj,1) = table2array(xunc.Coefficients(2,1));

            % plot settings
            title(allSubj(jj))
            ylabel('Net metabolic power (W)')
            ylim([0 5])
            xlim([0 600])
            xlabel('Time in exo (mins)')
            xticks(0:200:600)

        end
    end
end

%% group analysis

% save sensitivity analysis
if length(Fs_timescale) > 1
    cd(strcat('/Volumes/GoogleDrive/My Drive/Locomotion Lab/Projects/Cost Optimization - exos/Katies Data/Learning/Saving'))
    save('sensitivityAnalysis.mat','SFstdVAL','SFstdVAL_nw','SFstdVALnorm', ...
        'solTotalStdVAL','solTotalStdVAL_nw','solTotalStdVALnorm', ...
        'gastrocTotalStdVAL','gastrocTotalStdVAL_nw','gastrocTotalStdVALnorm')
end

close all

%% Figure 2. Changes in variability as participants gain experience with exoskeleton

disp('---------------------STATS FOR FIGURE 2---------------------')
pl2 = Fig2_Exploration_BarGraph_RepSubj(SFstdVAL_nw, SFstdVAL_zt, SFstdVAL, SFmeanVAL_nw, AnkleStdVAL_zt, AnkleVAL_zt, AnkleStdVAL, AnkleVAL, solTotalStdVAL_nw, solTotalStdVAL, solTotalMeanVAL_nw, gastrocTotalStdVAL_nw, gastrocTotalStdVAL, gastrocTotalMeanVAL_nw, orig, formatGraphs, SubjColors, Fs_timescale);

%% Figure 3. Adaptation along variables that affect energetic cost

disp('---------------------STATS FOR FIGURE 3---------------------')
pl3 = Fig3_MetCost_Adapt(yAdapt, metVAL, metVAL_zt, SFmeanVALnorm, AnkleVALnorm, AnkleVAL, AnkleVAL_zt, solTotalMeanVALnorm, gastrocTotalMeanVALnorm, orig, formatGraphs, SubjColors);

%% Figure 4. Timescales of variability and adaptation

disp('---------------------STATS FOR FIGURE 4 (ADAPTATION)---------------------')
% adaptation
pl5 = Fig4_Timescales(yAdapt, timepointsVAL, SubjColors, formatGraphs, 2, orig);

disp('---------------------STATS FOR FIGURE 4 (VARIABILITY)---------------------')
% variability
pl5 = Fig4_Timescales(yStd, timepointsVAL, SubjColors, formatGraphs, 1, orig);

%% Figure 5. Differences in timescales of variability and adaptation between variables

disp('---------------------STATS FOR FIGURE 5 (ADAPTATION)---------------------')

loadi = 1; % 0 = load previous simulations, 1 = run 10000 simulations
if loadi == 0
    timeConstantsAdapt = resamplingMethod(yAdapt, timepointsVAL, SubjColors, formatGraphs, 2, orig);
else
    cd('/Volumes/GoogleDrive/My Drive/Locomotion Lab/Projects/Cost Optimization - exos/Katies Data')
    if orig == 1
        load('adaptMC.mat')
        dim = [1 3 4];
    else
        load('adaptMC_otherdim.mat')
        dim = [3 4];
    end
end

% plotting adaptation TCs
figure(4501)
ax2=subplot(1,2,1);
hold on
boxplot(timeConstantsAdapt,'Symbol','')
ylim([0 400])
yticks(0:100:400)
title('Time constant of adaptation')
xticklabels({'SF','Ankle','Sol','Gastroc'})
ylabel('TC (minutes)')

% stats - adaptation
medianAdapt = nanmedian(timeConstantsAdapt);
iqrAdapt = quantile(timeConstantsAdapt,[0.25 0.75]);
[p,~,statsAdapt] = kruskalwallis(timeConstantsAdapt(:,dim));
pvalsAdapt = multcompare(statsAdapt,'CType','dunn-sidak');
pvalsAdaptDisp = pvalsAdapt(:,end);

% table
Rows = {'median';'iqr low';'iqr high'};
SF = [medianAdapt(1); iqrAdapt(1,1); iqrAdapt(2,1)];
Ankle = [medianAdapt(2); iqrAdapt(1,2); iqrAdapt(2,2)];
Sol = [medianAdapt(3); iqrAdapt(1,3); iqrAdapt(2,3)];
Gastroc = [medianAdapt(4); iqrAdapt(1,4); iqrAdapt(2,4)];
adaptSimTC = table(Rows,SF,Sol,Gastroc);
disp(adaptSimTC)
disp(pvalsAdaptDisp)

%% additional analyses
if orig == 0
    load('adaptMC.mat')
    % only consider +ve time constants
    timeConstantsAdapt_new = nan(10000,4);
    for i = 1:4
        temp = timeConstantsAdapt(:,i);
        temptemp = temp(temp>0);
        timeConstantsAdapt_new(1:length(temptemp),i) = temptemp;
    end
    x1 = timeConstantsAdapt_new;
    load('adaptMC_otherdim.mat')
    timeConstantsAdapt_new = nan(10000,4);
    for i = 1:4
        temp = timeConstantsAdapt(:,i);
        temptemp = temp(temp>0);
        timeConstantsAdapt_new(1:length(temptemp),i) = temptemp;
    end
    x2 = timeConstantsAdapt_new;
    timeConstantsAdapt_newnew = [x1(:,3) x2(:,3:4)];
    [pAdapt,tblAdapt,statsAdapt] = anova1(timeConstantsAdapt_newnew);
    pvalsAdapt=multcompare(statsAdapt);
    
    % table
    Rows = {'pval'};
    p1 = pvalsAdapt(1,end);
    p2 = pvalsAdapt(2,end);
    p3 = pvalsAdapt(3,end);
    compOridAddit = table(Rows,p1,p2,p3);
    disp(compOridAddit)
end

%% variability

disp('---------------------STATS FOR FIGURE 5 (VARIABILITY)---------------------')


if loadi == 0
    timeConstantsExplore = resamplingMethod(yStd, timepointsVAL, SubjColors, formatGraphs, 1, orig);
else
    cd('/Volumes/GoogleDrive/My Drive/Locomotion Lab/Projects/Cost Optimization - exos/Katies Data')
    if orig == 1
        load('exploreMC.mat')
        dim = [1 3 4];
    else
        load('exploreMC_otherdim.mat')
        dim = [3 4];
    end
end

% plotting variability TCs
figure(4501)
ax2=subplot(1,2,1);
hold on
boxplot(timeConstantsExplore,'Symbol','')
ylim([0 400])
yticks(0:100:400)
title('Time constant of variability')
xticklabels({'SF','Ankle','Sol','Gastroc'})
ylabel('TC (minutes)')

% stats - variability
medianExplore = nanmedian(timeConstantsExplore);
iqrExplore = quantile(timeConstantsExplore,[0.25 0.75]);
[p,~,statsExplore] = kruskalwallis(timeConstantsExplore(:,dim));
pvalsExplore = multcompare(statsExplore,'CType','dunn-sidak');
pvalsExploreDisp = pvalsExplore(:,end);

% table
Rows = {'median';'iqr low';'iqr high'};
SF = [medianExplore(1); iqrExplore(1,1); iqrExplore(2,1)];
Ankle = [medianExplore(2); iqrExplore(1,2); iqrExplore(2,2)];
Sol = [medianExplore(3); iqrExplore(1,3); iqrExplore(2,3)];
Gastroc = [medianExplore(4); iqrExplore(1,4); iqrExplore(2,4)];
ExploreSimTC = table(Rows,SF,Ankle,Sol,Gastroc);
disp(ExploreSimTC)
disp(pvalsExploreDisp)

%% additional analyses
if orig == 0
    load('exploreMC.mat')
    % only consider +ve time constants
    timeConstantsExplore_new = nan(10000,4);
    for i = 1:4
        temp = timeConstantsExplore(:,i);
        temptemp = temp(temp>0);
        timeConstantsExplore_new(1:length(temptemp),i) = temptemp;
    end
    x1 = timeConstantsExplore_new;
    load('exploreMC_otherdim.mat')
    timeConstantsExplore_new = nan(10000,4);
    for i = 1:4
        temp = timeConstantsExplore(:,i);
        temptemp = temp(temp>0);
        timeConstantsExplore_new(1:length(temptemp),i) = temptemp;
    end
    x2 = timeConstantsExplore_new;
    timeConstantsExplore_newnew = [x1(:,3) x2(:,3:4)];
    [pExplore,tblExplore,statsExplore] = anova1(timeConstantsExplore_newnew);
    pvalsExplore=multcompare(statsExplore);
    
    % table
    Rows = {'pval'};
    p1 = pvalsExplore(1,end);
    p2 = pvalsExplore(2,end);
    p3 = pvalsExplore(3,end);
    compOridAddit = table(Rows,p1,p2,p3);
    disp(compOridAddit)
end

%% Checking for differences within NW

vars = {'SFmeanVAL_nw';'AnkleVAL_zt';'solTotalMeanVAL_nw';'gastrocTotalMeanVAL_nw'};
for i = 1:4
    y = eval(string(vars(i)));
    if i > 2 % muscle activity
        init = nanmean(y(1:2,1:10));
        fin = nanmean([y(9:10,1) y(end-1:end,2:10)]);
    else % other vars
        init = nanmean([y(1:2,1:5) y(3:4,6) y(1:2,7:10)]);
        fin = nanmean([y(9:10,1) y(end-1:end,2:10)]);
    end
    [~,padapt_nw(i)] = ttest(init,fin);
    percChange_nw = ((fin-init)./init).*100;
    percChange_nwMean(i) = nanmean(percChange_nw);
    percChange_nwStd(i) = nanstd(percChange_nw);
    
    clearvars init fin percChange_nw
end

% table
Rows = {'Baseline';'Percent Change mean';'Percent Change std'};
SF = [padapt_nw(1); percChange_nwMean(1); percChange_nwStd(1)];
Ankle = [padapt_nw(2); percChange_nwMean(2); percChange_nwStd(2)];
Sol = [padapt_nw(3); percChange_nwMean(3); percChange_nwStd(3)];
Gastroc = [padapt_nw(4); percChange_nwMean(4); percChange_nwStd(4)];
pvalsEstablishedPolicy = table(Rows,SF,Ankle,Sol,Gastroc);
disp(pvalsEstablishedPolicy)

%% Checking for differences between ZT and NW

% adapt
varsAdapt_nw = {'SFmeanVAL_nw';'solTotalMeanVAL_nw';'gastrocTotalMeanVAL_nw'};
varsAdapt_zt = {'SFmeanVAL_zt';'solTotalMeanVAL_zt';'gastrocTotalMeanVAL_zt'};
% explore
varsExplore_nw = {'SFstdVAL_nw';'solTotalStdVAL_nw';'gastrocTotalStdVAL_nw'};
varsExplore_zt = {'SFstdVAL_zt';'solTotalStdVAL_zt';'gastrocTotalStdVAL_zt'};
for i = 1:6
    if i > 3
        y_nw = eval(string(varsExplore_nw(i-3)));
        y_zt = eval(string(varsExplore_zt(i-3)));
    else
        y_nw = eval(string(varsAdapt_nw(i)));
        y_zt = eval(string(varsAdapt_zt(i)));
    end
    nw = nanmean(y_nw);
    zt = nanmean(y_zt);
    
    [~,p_nw_zt(i)] = ttest(nw,zt);
    percChange_nw_zt = ((zt-nw)./nw).*100;
    percChange_nw_zt_Mean(i) = nanmean(percChange_nw_zt);
    percChange_nw_zt_Std(i) = nanstd(percChange_nw_zt);
    
    clearvars y_nw nw y_zt zt percChange_nw_zt
end

% adapt
Rows = {'MAGNITUDE nw vs zt';'Percent Change mean';'Percent Change std'};
SF = [p_nw_zt(1); percChange_nw_zt_Mean(1); percChange_nw_zt_Std(1)];
Sol = [p_nw_zt(2); percChange_nw_zt_Mean(2); percChange_nw_zt_Std(2)];
Gastroc = [p_nw_zt(3); percChange_nw_zt_Mean(3); percChange_nw_zt_Std(3)];
NWvsZT = table(Rows,SF,Sol,Gastroc);
disp(NWvsZT)
% explore
Rows = {'VARIABILITY nw vs zt';'Percent Change mean';'Percent Change std'};
SF = [p_nw_zt(4); percChange_nw_zt_Mean(4); percChange_nw_zt_Std(4)];
Sol = [p_nw_zt(5); percChange_nw_zt_Mean(5); percChange_nw_zt_Std(5)];
Gastroc = [p_nw_zt(6); percChange_nw_zt_Mean(6); percChange_nw_zt_Std(6)];
NWvsZT = table(Rows,SF,Sol,Gastroc);
disp(NWvsZT)