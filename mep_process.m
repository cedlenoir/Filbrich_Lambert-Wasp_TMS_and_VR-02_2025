%% 1) LOAD DATA AND SET A FEW PARAMETERS
% channelEMG = [1 2];
% channelACC = [4 5]; % 4 -> Y; 5 -> Z; (3 -> X)

% load data from displayed list of .mat files contained in data folder
OS = ispc;
if true(OS) % PC
    dataPath = '.\data';
    savePath = '.\results';
else % mac
    dataPath = './data';
    savePath = './results';
end
file_list = dir(fullfile(dataPath,'*.mat*'));
clc
for ifile = 1:size(file_list,1)
    disp(strcat(['#',num2str(ifile),' ',file_list(ifile).name]))
    disp(' ')
end

% dialog box
prompt = {'\fontsize{15} Subject ID? :','\fontsize{15} which file number #? : ','\fontsize{15} Condition? : PRE -> 1 OR POST -> 2',...
    '\fontsize{15} Hand side? : left -> 1 OR right -> 2','\fontsize{15} EMG trial plots? : YES -> 1 OR NO -> 0','\fontsize{15} ACCELERO trial plots? : YES -> 1 OR NO -> 0'};
dlgtitle = 'LOAD DATA';
opts.Interpreter = 'tex';
dims = repmat([1 60],6,1);
definput = {'000','','','2','0','0'};
info = inputdlg(prompt,dlgtitle,dims,definput,opts);
subject_id = char(info(1));
chosen_file = str2double(info(2));
cond = str2double(info(3));
if cond == 1
    condition = 'PRE';
elseif cond == 2
    condition = 'PST';
end
hand_side = str2double(info(4));
if hand_side == 1
    tested_hand = 'left';
elseif hand_side == 2
    tested_hand = 'right';
end
plot_opt_EMG = str2double(info(5));
plot_opt_ACC = str2double(info(6));
filename = file_list(chosen_file).name(1:end-4);
savename = strcat('sub-',subject_id,'-',condition,'-',tested_hand,'-hand_',filename);

% store data in structure "data"
try
    datastruct = load(file_list(chosen_file).name);
    storedvars = fieldnames(datastruct) ;
    FirstVarName = storedvars{2};
    data = datastruct.(FirstVarName);
    clc
    disp(strcat(['file loaded : ',file_list(chosen_file).name]))
catch
    disp('wave data not found!')
end

% extract numbers of trials and channels
trials = double(data.frames);

% design bandpass filter Butterworth 20-450 Hz; order 4
filtOrder = 4;
low_cutoff = 20;
high_cutoff = 450;
% sampling rate
sr = 2000;
fnyquist = sr/2;

% bandpass-filter the EMG data
[b,a] = butter(filtOrder,[low_cutoff,high_cutoff]./fnyquist,'bandpass');
data.emg_filt = filtfilt(b,a,data.values(:,1:2,:));

% store acc data
data.acc = data.values(:,4:5,:);

% define baseline window as the first 400 ms of the frame
baselineWind = 0.4;
% trigger position is the middle of the frame in Signal
trigger_idx = 0.5*sr;
% define response time window for ACCELERO, from trigger index  to +200 ms
responseACCWind = [trigger_idx trigger_idx+0.2*sr];
% define response time window for EMG, from trigger index  to +200 ms
responseEMGWind = [trigger_idx trigger_idx+0.04*sr];
% define MEP peak-to-peak amplitude criterion >50 uV
emg_threshold = 0.05;

pause(0.1)

%% 2) PREPROCESS EMG DATA

avg_emg_bl = zeros(trials,1);

for itrial = 1:trials
    for ichan = 1:2
        % baseline correction
        data.emg_filt(:,ichan,itrial) = data.emg_filt(:,ichan,itrial) - mean(data.emg_filt(1:baselineWind*sr,ichan,itrial));

        %%%%%%%%%%% check if there is backgroun noise / pre-stimulus muscular activity ?

        avg_emg_bl(ichan,itrial) = mean(data.emg_filt(1:baselineWind*sr,ichan,itrial));

        % find max and min peaks in the response time window, and 10 ms after TMS pulse to avoid artefact and get its latency
        [data.peak_max_amp_emg(itrial,ichan), data.peak_max_time_emg(itrial,ichan)] = max(squeeze(data.emg_filt((responseEMGWind(1)+(0.01*sr)):responseEMGWind(2),ichan,itrial)));
        [data.peak_min_amp_emg(itrial,ichan), data.peak_min_time_emg(itrial,ichan)] = min(squeeze(data.emg_filt((responseEMGWind(1)+(0.01*sr)):responseEMGWind(2),ichan,itrial)));
        data.peak_max_time_emg(itrial,ichan) = data.peak_max_time_emg(itrial,ichan) + responseEMGWind(1)+(0.01*sr)-1;
        data.peak_min_time_emg(itrial,ichan) = data.peak_min_time_emg(itrial,ichan) + responseEMGWind(1)+(0.01*sr)-1;

        data.p2p_amp_emg(itrial,ichan) = abs(data.peak_max_amp_emg(itrial,ichan)) + abs(data.peak_min_amp_emg(itrial,ichan));

        % check the duration between min and max, and warning if greater
        % than 15 ms
        data.p2p_time_intv(itrial,ichan) = abs(data.peak_max_time_emg(itrial,ichan) - data.peak_min_time_emg(itrial,ichan));
        if data.p2p_time_intv(itrial,ichan) > 0.015*sr
            data.p2p_amp_emg(itrial,ichan) = NaN;
            data.MEP_latency(itrial,ichan) = NaN;
        else
        end

        % check if latency of MEPs is valid (< 40 ms)
        if ~isnan(data.p2p_amp_emg(itrial,ichan))

            if data.peak_max_time_emg(itrial,ichan) < responseEMGWind(2) && data.peak_max_time_emg(itrial,ichan) < data.peak_min_time_emg(itrial,ichan)
                data.p2p_amp_emg(itrial,ichan) = data.p2p_amp_emg(itrial,ichan);
                data.MEP_latency(itrial,ichan) = data.peak_max_time_emg(itrial,ichan);
            elseif data.peak_min_time_emg(itrial,ichan) < responseEMGWind(2) && data.peak_min_time_emg(itrial,ichan) < data.peak_max_time_emg(itrial,ichan)
                data.p2p_amp_emg(itrial,ichan) = data.p2p_amp_emg(itrial,ichan);
                data.MEP_latency(itrial,ichan) = data.peak_min_time_emg(itrial,ichan);
            else
                data.p2p_amp_emg(itrial,ichan) = NaN;
                data.MEP_latency(itrial,ichan) = NaN;
            end

        elseif isnan(data.p2p_amp_emg(itrial,ichan))
            data.MEP_latency(itrial,ichan) = NaN;
        end

        % check if the MEP p2p amplitude are > than threshold
        if ~isnan(data.p2p_amp_emg(itrial,ichan)) && data.p2p_amp_emg(itrial,ichan) > emg_threshold
            data.p2p_amp_emg(itrial,ichan) = data.p2p_amp_emg(itrial,ichan);
        else
            data.p2p_amp_emg(itrial,ichan) = NaN;
            data.MEP_latency(itrial,ichan) = NaN;
        end

    end
end

% list the indices of valid MEPs among all trials
for ichan= 1:2
    for itrial = 1:trials
        idxMEP(itrial,ichan) = ~isnan(data.p2p_amp_emg(itrial,ichan));
    end
    valid_MEP_trials{:,ichan} = find(idxMEP(:,ichan));
end

%% Plot EMG and peaks

if plot_opt_EMG == 1
    % remove trials with no MEP

    for ichan = 1:2
        for itrial = 1:length(valid_MEP_trials{1,ichan})

            figure('Position',[600,0,600,500],'Color','w');
            %     subplot(2,1,1)
            pe1 = plot(data.emg_filt(:,ichan,valid_MEP_trials{1,ichan}(itrial)),'b');
            title(strcat('EMG-',file_list(chosen_file).name(1:end-4),'-chan-"',data.chaninfo(1).title,'"-trial#',num2str(valid_MEP_trials{1,ichan}(itrial))))
            hold on
            pe2 = plot(data.peak_max_time_emg(valid_MEP_trials{1,ichan}(itrial),ichan),data.peak_max_amp_emg(valid_MEP_trials{1,ichan}(itrial),ichan),'xr','MarkerSize',20);
            pe3 = plot(data.peak_min_time_emg(valid_MEP_trials{1,ichan}(itrial),ichan),data.peak_min_amp_emg(valid_MEP_trials{1,ichan}(itrial),ichan),'xr','MarkerSize',20);
            ax= gca;
            ax.Box = 'off';
            ax.XLabel.String = 'time(ms)';
            ax.YLabel.String = strcat('amplitude (mV)');
            ax.XTickLabel = {'-500','-400','-300','-200','-100','0','100','200','300','400','500'};
            xline(responseEMGWind(2),'k--','response window','LabelHorizontalAlignment','right')
            xline(trigger_idx,'r--','TMS','LabelHorizontalAlignment','left')
            legend(pe3,'peak-to-peak','box','off')
        end
    end

    pause(1)
    disp('To continue and close all figures press any key !')
    pause()
    close all; clc

else
end


%% plot all EMG traces
for ichan = 1:2
    for itrial = 1:trials
        figure('Position',[0,0,600,500],'Color','w');
        pe4 = plot(data.emg_filt(:,ichan,itrial),'r');
        title(strcat('EMG-',file_list(choose_file).name(6:end-4),'-chan#',num2str(ichan),'-',data.chaninfo(2).title,'-trial#',num2str(itrial)))
        hold on
        pe5 = plot(data.peak_max_time_emg(itrial,ichan),data.peak_max_amp_emg(itrial,ichan),'xr','MarkerSize',20);
        pe6 = plot(data.peak_min_time_emg(itrial,ichan),data.peak_min_amp_emg(itrial,ichan),'xr','MarkerSize',20);
        ax= gca;
        ax.Box = 'off';
        ax.XLabel.String = 'time(ms)';
        ax.YLabel.String = strcat('amplitude (mV)');
        xline(responseEMGWind(1),'k--')
        xline(responseEMGWind(2),'k--')
        xline(trigger_idx,'r-','TMS','LabelHorizontalAlignment','left')
        legend(pe4,'ALL','box','off')
    end
end