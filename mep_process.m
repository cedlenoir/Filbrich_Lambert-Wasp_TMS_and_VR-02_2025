% 
% 
% decide on folder and subject level naming of the acquisitions
% sub-00$
% main_folder contains the raw_data folder, pre-processed data folder,
% results folder, (scripts folder?)
% 
% (sub-001 is test_giulia)
% (sub-002 is Laurie)
%
% % Cédric Lenoir, NeTMeD, IoNS, UCLouvain, January 2026
% 
%% 1) LOAD DATA AND SET A FEW PARAMETERS

% set the paths and initialize letswave
% Select the general data folder which contains all subjects folders
main_folder = uigetdir('C:\','Select main folder');
cd(main_folder)
raw_data_folder = fullfile(main_folder,'raw_data');
results_folder = fullfile(main_folder,'results');
pre_proc_folder = fullfile(main_folder,'pre-processed');

if ~isfolder(raw_data_folder); mkdir(raw_data_folder); end
if ~isfolder(results_folder); mkdir(results_folder); end
if ~isfolder(pre_proc_folder); mkdir(pre_proc_folder); end

addpath(genpath(main_folder))

% initialize letswave 6
% check if lw is already on the path if not select the folder of the toolobx
if contains(path, 'C:\Users\cedlenoir\Documents\MATLAB\letswave6-master')
    letswave();
    clc
else
    lw_path = uigetdir('C:\Users\cedlenoir\Documents\MATLAB','Select Letswave 6 folder');
    addpath(genpath(lw_path))
    letswave();
    clc
end

% dialog box
prompt = {'\fontsize{12} Subject ID? :','\fontsize{12} Sensitized arm (L/R): ','\fontsize{12} Stimulated hemisphere (L/R):',...
    '\fontsize{12} TMS PRE-CAPS start index? : ','\fontsize{12} TMS PRE-CAPS stop index? : ',...
    '\fontsize{12} TMS POST-CAPS start index? : ','\fontsize{12} TMS POST-CAPS stop index? : ',...
    '\fontsize{12} TMS POST-CAPS start index? : ','\fontsize{12} TMS POST-CAPS stop index? : ',...
    '\fontsize{12} TMS POST-CAPS start index? : ','\fontsize{12} TMS POST-CAPS stop index? : ',...
    '\fontsize{12} TMS POST-CAPS start index? : ','\fontsize{12} TMS POST-CAPS stop index? : ',...
    '\fontsize{12} EMG trial plots? : YES -> 1 OR NO -> 0','\fontsize{12} Comments: ',};
dlgtitle = 'LOAD SUBJECT DATA';
opts.Interpreter = 'tex';
dims = repmat([1 80],15,1);
definput = {'002','R','L','114','137','138','164','165','192','193','220','221','246','1',''};
info = inputdlg(prompt,dlgtitle,dims,definput,opts);
subject_id = char(info(1));
sensi_arm = char(info(2));
stim_hemi = char(info(3));
ixd_tms_pre_start = str2double(info(4));
ixd_tms_pre_stop = str2double(info(5));
ixd_tms_pst_start = str2double([info(6) info(8) info(10) info(12)]);
ixd_tms_pst_stop = str2double([info(7) info(9) info(11) info(13)]);
plot_opt_EMG = str2double(info(14));
notes = str2double(info(15));

% store info in structure
sub_info = struct;
sub_info.sub_ID = subject_id;
sub_info.sensitized_arm = sensi_arm;
sub_info.stimulated_hemisph = stim_hemi;
sub_info.TMS_triggers.pre = [ixd_tms_pre_start;ixd_tms_pre_stop];
sub_info.TMS_triggers.pst = [ixd_tms_pst_start;ixd_tms_pst_stop];
sub_info.comments = notes;

% folder name of Visor EMG data
subject_folder = fullfile(raw_data_folder,strcat('sub-',subject_id),'Sessions');
sub_pre_proc_folder = fullfile(pre_proc_folder,sprintf('sub-%s',subject_id));
if ~isfolder(sub_pre_proc_folder); mkdir(sub_pre_proc_folder); end

% list sessions
session_list = dir(subject_folder);
session_list = session_list(~ismember({session_list.name}, {'.', '..'}));

% list EMG CNT files
file_list = dir(fullfile(session_list.folder,session_list.name,'*emg.cnt*'));

% import CNT file
[out_data,~] = RLW_import_CNT(fullfile(session_list.folder,session_list.name,file_list.name));
filename = file_list.name(1:end-8);
savename = strcat('sub-',subject_id,'_',filename);
out_data.header.name = savename;
CLW_save(sub_pre_proc_folder,out_data.header,out_data.data);
clear out_data
[header, data] = CLW_load(fullfile(sub_pre_proc_folder,savename));
sr = 1/header.xstep;

% pre-processign steps in lw
% high pass filter Butterworth 4 Hz; order 4
low_cutoff = 4;
order = 4;
[filt_header, filt_data] = RLW_butterworth_filter(header,data,'filter_type','highpass','low_cutoff',low_cutoff,'filter_order',order);

% segmentation
xstart = -0.2;
xduration = 0.7;
[seg_header, seg_data] = RLW_segmentation(filt_header, filt_data, {'1'},'x_start',xstart,'x_duration',xduration);
seg_header.chanlocs.labels = 'EMG1';

% DC & Linear detrend
[dc_header, dc_data] = RLW_dc_removal(seg_header,seg_data,'linear_detrend',1);
dc_header.name = strcat(header.name,' HPfilt ep DC');

% save dataset
CLW_save(sub_pre_proc_folder,dc_header, dc_data);

% keep TMS events after rMT
% arrange epochs of valid TMS triggers for each TMS block
valid_tms_idx = cell(5,1);
for iblock = 1:5
    if iblock == 1
        valid_tms_idx{iblock,1} = (ixd_tms_pre_start:ixd_tms_pre_stop);
    else
        valid_tms_idx{iblock,1} = (ixd_tms_pst_start(iblock-1):ixd_tms_pst_stop(iblock-1));
    end
end
for iblock = 1:5
    [block_header{iblock,1}, block_data{iblock,1}] = RLW_arrange_epochs(dc_header, dc_data,valid_tms_idx{iblock,1});
    block_header{iblock,1}.name = strcat(dc_header.name,' BLK ',num2str(iblock));
    CLW_save(sub_pre_proc_folder,block_header{iblock,1}, block_data{iblock,1});
end

% define baseline window to check for baseline activity
bsln_time_window = [-0.2 0];
bsln_sample_window = [1 round(0.2*sr+1)];

% First check: for baseline single trial EMG activity (RMS) per block
% threshold is mean(RMS) + 2.5*std(RMS)
% (as in Sulcova D. et al. bioRxiv 2022; Grandjean and Duque NIMG 2020)
idx_exclud = cell(5,1);
for iblock = 1:5
    temp_data = squeeze(block_data{iblock,1});

    % compute RMS for each trial
    for itrial = 1:size(temp_data,1)
        bsln_rms{iblock,1}(itrial,1) = rms(temp_data(itrial,bsln_sample_window(1):bsln_sample_window(2)));
    end
    % threshold mean RMS +/- 3*SD (permissive)
    avg_bsln_rms(iblock,1) = mean(bsln_rms{iblock,1});
    sd_bsln_rms(iblock,1) = std(bsln_rms{iblock,1});
    thrshld(iblock,1) = avg_bsln_rms(iblock,1)+2.5*sd_bsln_rms(iblock,1);
    % loop until no outlier is found (exclude MEP + warning if baseline trial
    % RMS > mean RMS +/- 2.5*SD, until no outliers)
    idx_val = 1;
    idx_exc = 1;
    for itrial = 1:size(temp_data,1)
        if bsln_rms{iblock,1}(itrial,1) < thrshld(iblock,1)
            valid_data{iblock,1}(idx_val,:) = temp_data(itrial,:);
            idx_val = idx_val+1;
        else
            idx_exclud{iblock,1}(idx_exc,1) = itrial;
            idx_exc = idx_exc+1;
        end
    end
    % if at least one MEP has been excluded redo the RMS check
    if size(valid_data{iblock,1},1) < block_header{iblock,1}.datasize(1)

        for itrial = 1:size(valid_data{iblock,1},1)
            bsln_rms2{iblock,1}(itrial,1) = rms(valid_data{iblock,1}(itrial,bsln_sample_window(1):bsln_sample_window(2)));
        end
    % threshold mean RMS +/- 3*SD (permissive)
    avg_bsln_rms2(iblock,1) = mean(bsln_rms2{iblock,1});
    sd_bsln_rms2(iblock,1) = std(bsln_rms2{iblock,1});
    thrshld2(iblock,1) = avg_bsln_rms2(iblock,1)+2.5*sd_bsln_rms2(iblock,1);
    % loop until no outlier is found (exclude MEP + warning if baseline trial
    % RMS > mean RMS +/- 3*SD, until no outliers)
    idx_val = 1;
    idx_exc = 1;
    for itrial = 1:size(valid_data{iblock,1},1)
        if bsln_rms2{iblock,1}(itrial,1) < thrshld2(iblock,1)
            valid_data{iblock,1}(idx_val,:) = valid_data{iblock,1}(itrial,:);
            idx_val = idx_val+1;
        else
            idx_exclud{iblock,1}(idx_exc,1) = itrial;
            idx_exc = idx_exc+1;
        end
    end
    else
    end
    clear temp_data
end

% display the number of MEP excluded for each block
clc
for iblock = 1:5
    if isempty(idx_exclud{iblock,1})
    disp(strcat(['No MEPs exclude in block ',num2str(iblock)]))
    else
        disp(strcat([num2str(length(idx_exclud{iblock,1})),' MEPs excluded in block ',num2str(iblock)]))
    end
end

% Second check: threshold on baseline RMS exceeding +/-15 µV
% (as in Sulcova et al. BioRxiv 2022; Morozova et al. Sci Reports 2024)
abs_thrshld = 15;
idx_exclud = cell(5,1);
for iblock = 1:5
    for itrial = 1:size(valid_data{iblock,1},1)
        bsln_rms2{iblock,1}(itrial,1) = rms(valid_data{iblock,1}(itrial,bsln_sample_window(1):bsln_sample_window(2)));
    end

    idx_val = 1;
    idx_exc = 1;
    for itrial = 1:size(valid_data{iblock,1},1)
        if bsln_rms2{iblock,1}(itrial,1) < abs_thrshld
            valid_data2{iblock,1}(idx_val,:) = valid_data{iblock,1}(itrial,:);
            idx_val = idx_val+1;
        else
            idx_exclud{iblock,1}(idx_exc,1) = itrial;
            idx_exc = idx_exc+1;
        end
    end
end
% plot MEP and criteria
xval = -0.2:1/sr:0.5;
for iblock = 1:5
    for itrial = 1:size(idx_exclud{iblock,1},1)
        f = figure('Color','w','Position',[0 0 1500 700]);
        plot(xval,valid_data{iblock,1}(idx_exclud{iblock,1}(itrial,1),:),'k','LineWidth',1)
        xline(0,'--r','TMS','LabelOrientation','horizontal')
        yline(0,'-k')
        yline(bsln_rms2{iblock,1}(idx_exclud{iblock,1}(itrial,1),1),'--k','RMS','LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','left')
        yline(thrshld(iblock,1),'--g','threshold','LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','center')
        yline(abs_thrshld,'--m','absolute threshold','LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','right')
        ax = gca;
        ax.TickDir = 'out';
        ax.XLabel.String = 'time (s)';
        ax.YLabel.String = 'amplitude (µV)';
        ax.YLim = [-200 200];
        ax.Box = 'off';
        title(strcat(['sub-',subject_id,' MEP#',num2str(idx_exclud{iblock,1}(itrial,1)),' in block ',num2str(iblock)]))
        pause()
    end
end

% define MEP peak-to-peak amplitude criterion > 50 µV in specific window
emg_threshold = 50;
% extract MEP latencies




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