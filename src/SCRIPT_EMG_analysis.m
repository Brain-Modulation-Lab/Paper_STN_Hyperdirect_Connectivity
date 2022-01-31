% SCRIPT EMG analysis

% compare EMG data with ECoG data
%
% Ahmed Jorge Feb 2019


clear all; close all; clc







%% DBS3025
% overall difficult to sync times from neuroomega to ecog. Witek will look
% at it later.


%% Read emg table info for DBS3025

% read table
clc
cfg             = [];
cfg.chantype    = 'emg';
cfg.path        = '/Volumes/Nexus/DBS/DBS3025/Raw/NeuroOmega/Matlab';
cfg.mpx_path    = '/Volumes/Nexus/DBS/DBS3025/Raw/NeuroOmega';
no_info = bml_neuroomega_info_raw(cfg);

% save table (to avoid long loading times)

%% get actual EMG signal for DBS3025

% manually load table (EMG_no_info_DBS3025.mat)
% !!! you have to manually load it!!!

%load related tables
SUBJECT = 'DBS3025';
cfg             = [];
cfg.chantype    = 'emg';
cfg.path        = '/Volumes/Nexus/DBS/DBS3025/Raw/NeuroOmega/Matlab';
cfg.mpx_path    = '/Volumes/Nexus/DBS/DBS3025/Raw/NeuroOmega';
sync      = bml_annot_read(['annot/' SUBJECT '_sync.txt']);
session   = bml_annot_read(['annot/' SUBJECT '_session.txt']);
electrode = bml_annot_read(['annot/' SUBJECT '_electrode.txt']);

% set up parameters
cfg=[];
cfg.roi=no_info(no_info.depth == -1.998,:); %-4.001 depth is indicated in the experiment log
cfg.chantype = 'emg';
cfg.filetype = 'neuroomega.mat';
cfg.allow_missing = true;

%not sure what this is for
electrode_emg_stim = electrode(electrode.type=="emg",:); 
cfg.electrode = electrode_emg_stim;
%cfg.electrode.starts = repmat(cfg.roi.starts(1),size(cfg.electrode.starts));
%cfg.electrode.ends = repmat(cfg.roi.ends(end),size(cfg.electrode.starts));
cfg.timetol = 0.001; %increase if necessary
no_emg = bml_load_continuous(cfg);

%% read ECoG data

% load raw data at 1kHz to identify sets
clear raw_all address
pt= 3025;
cfg4 = [];
cfg4.pt = pt;
cfg4.desired_f_sample = 1;
cfg4.plot_databrowser = 0;
[raw_all, ~] = ecog_getraw(cfg4);  

% mark stim trials (rough approximation of stim time done at 1kHz)
clear stim_marker_1khz stim_marker_time_1khz
cfg4 =[];
cfg4.pt = pt;
cfg4.sett = 999;
cfg4.flag = 1;
[stim_marker_1khz, stim_marker_time_1khz] = define_trials_ecog(cfg4,raw_all);



%% 

clc

timeoi = [ -0.100 1.100];
setoi  = 4;
setoia = find(stim_marker_1khz(:,2) == setoi);

pc=1;
for stim = setoia(1):setoia(end)
        
    indx = find(no_emg.time{1}>stim_marker_time_1khz(stim,1)-0.0001 & no_emg.time{1}<stim_marker_time_1khz(stim,1)+0.0001);
    boi = indx(1)+timeoi(1)*no_emg.fsample:indx(1)+timeoi(2)*no_emg.fsample;
    x1 = no_emg.time{1}(boi);
    y1 = no_emg.trial{1}(:,boi);
    
    clear indx
    indx = find(raw_all.time{1}>stim_marker_time_1khz(stim,1)-0.0001 & raw_all.time{1}<stim_marker_time_1khz(stim,1)+0.0001);
    toi = indx(1)+timeoi(1)*1000:indx(1)+timeoi(2)*1000;
    x2 = raw_all.time{1}(toi);
    y2 = raw_all.trial{1}(1,toi);
    
    figure(51)
    subplot(6,5,pc); pc=pc+1;
    for i = 1:4
        plot(x1,smooth(y1(i,:),40),'k'); hold on
    end
    plot(x2,y2,'r','LineWidth',3); hold on
    ylim([-10000 10000])
    
    
end












%% ------------------------------------------------------------------------
% coded with Alan (deprecrate?)

cd('/Volumes/Nexus/DBS/DBS3025/Preprocessed Data/Sync/')

SUBJECT = 'DBS3025';
sync      = bml_annot_read(['annot/' SUBJECT '_sync.txt']);
session   = bml_annot_read(['annot/' SUBJECT '_session.txt']);
electrode = bml_annot_read(['annot/' SUBJECT '_electrode.txt']);

NEW_PATH = '/Volumes/Nexus/DBS';
sync.folder = strrep(sync.folder, '\\136.142.16.9\Nexus\DBS', NEW_PATH);
sync.folder = strrep(sync.folder, '\', '/');

fprintf('passed 91\n')

%loading neuroomega macro
cfg=[];
cfg.chantype = 'emg';
cfg.roi = sync(sync.filetype=="neuroomega.mat" & sync.chantype=="analog",:);
sync_no_macro=bml_sync_transfer_neuroomega_chantype(cfg);

fprintf('passed 99\n')

cfg=[];
cfg.roi=sync_no_macro(end-1,:);
cfg.chantype = 'emg';
cfg.allow_missing = true;
cfg.electrode = electrode(electrode.type=="emg",:);
cfg.timetol = 0.001; %increase if necessary
no_emg = bml_load_continuous(cfg);

fprintf('passed 109\n')

%% convert to table
info_ns2 = bml_roi_table(no_emg,'trellis');

fprintf('passed 115\n')

 
% load actual raw
[raw, ~] = bml_load_continuous(cfg);

%% ------------------------------------------------------------------------









































%% ------------------------------------------------------------------------
% DBS 3031 (ECoG and EMG on the same file)


SUBJECT='DBS3031';
DATE=datestr(now,'yyyymmdd');
PATH_DATA='/Volumes/Nexus/DBS';
PATH_TRELLIS=[PATH_DATA filesep SUBJECT filesep 'Raw' filesep 'Trellis'];
%sync = bml_annot_read([PATH_DATA '/' SUBJECT '/Preprocessed Data/Sync/annot/' SUBJECT '_sync.txt']);
%sync(sync.filetype=="trellis.ns5",:)
cfg=[];
cfg.path=PATH_TRELLIS;
cfg.pattern='*.ns2';                %could be either .ns2=1kHz or .ns5=30kHz
cfg.filetype='trellis.ns2';         %could be either .ns2=1kHz or .ns5=30KHz
info_ns2 = bml_info_raw(cfg);
info_ns2 = bml_roi_table(info_ns2,'trellis');
cfg.chantype = 'lfp';
cfg.roi = info_ns2(7,:);
[raw, ~] = bml_load_continuous(cfg);

clear stim_marker_1khz stim_marker_time_1khz
cfg4 =[];
cfg4.pt = pt;
cfg4.sett = 999;
cfg4.flag = 1;
[stim_marker_1khz, stim_marker_time_1khz] = define_trials_ecog(cfg4,raw_all);

 
%% plot trial by trial EMG
clc

timeoi = [ -0.050 0.150];
setoi  = 1;
setoia = find(stim_marker_1khz(:,2) == setoi);


% emg
indx_1 = find(raw.label=="lfp 265");
indx_2 = find(raw.label=="lfp 267");

pc=1;
for stim = 101:130%setoia(1):setoia(1)+29    %setoia(end)
        
    clear indx
    indx = find(raw_all.time{1}==stim_marker_time_1khz(stim,1));
    toi = indx(1)+timeoi(1)*1000:indx(1)+timeoi(2)*1000;
    x1 = raw_all.time{1}(toi);
    y1 = raw_all.trial{1}(indx_1,toi);
    y2 = raw_all.trial{1}(indx_2,toi);
    y3 = y2-y1;
    
    %N.B. 149 shows the stim really nicely
    clear indx
    indx = find(raw_all.time{1}==stim_marker_time_1khz(stim,1));
    toi = indx(1)+timeoi(1)*1000:indx(1)+timeoi(2)*1000;
    x9 = raw_all.time{1}(toi);
    y9 = raw_all.trial{1}(149,toi);
    
    figure(51)
    subplot(6,5,pc); pc=pc+1;
    plot(x1,y1,'b'); hold on
    plot(x1,y2,'g'); hold on
    plot(x1,y3,'k'); hold on
    
    plot(x9,y9,'r'); hold on

    ylim([-500 500])
  
end



%% 
pc=1;
for stim = setoia(1):setoia(end)
        
    clear indx
    indx = find(raw_all.time{1}==stim_marker_time_1khz(stim,1));
    toi = indx(1)+timeoi(1)*1000:indx(1)+timeoi(2)*1000;
    x1 = linspace(timeoi(1),timeoi(2),(abs(timeoi(2))+abs(timeoi(1)))*1000+1);
    y1 = raw_all.trial{1}(indx_1,toi);
    y2 = raw_all.trial{1}(indx_2,toi);
    y3 = y2-y1;
    
    
    if abs(max(y3))<150
        y4(stim,:) = y3;
    else
        y4(stim,:) = nan*length(y3);
    end
    
    %{
    %N.B. 149 shows the stim really nicely
    clear indx
    indx = find(raw_all.time{1}>stim_marker_time_1khz(stim,1)-0.0001 & raw_all.time{1}<stim_marker_time_1khz(stim,1)+0.0001);
    toi = indx(1)+timeoi(1)*1000:indx(1)+timeoi(2)*1000;
    x2 = raw_all.time{1}(toi);
    y2 = raw_all.trial{1}(149,toi);
    %}
    
    %figure(51)
    %plot(x1,y3,'k'); hold on
    %plot(x2,y2,'r','LineWidth',1); hold on
    %ylim([-150 150])
  
end


figure(55)
plot(x1,nanmedian(y4,1)); hold on
y41 = quantile(y4,0.25);
y43 = quantile(y4,0.75);
plot(x1,y41)
plot(x1,y43)


figure(56)
plot(x1,smooth(nanmean(y4,1),1)); hold on
%y41 = smooth(mean(y4,1),10)+std(y4)';
%y43 = std(y4);
%plot(x1,y41)
%plot(x1,y43)



















