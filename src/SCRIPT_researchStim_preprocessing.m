% ECoG Research Stim Preprocessing
% 
% this script will convert raw .ns5 ecog data into a processed trial by
% trial dataset (pro)
%
% 
% raw data --> get set times --> get pro data (trial by trial)
% 1 set = 30 trials (each trial is one stimulation done, each set corresponds to a particular tract (central, posterior, medial, etc)
%
% ajorge 2018


%% loading packages

%fieldtrip
addpath('/Users/brainmodulationlab/git/fieldtrip')
addpath('/Users/brainmodulationlab/git/NPMK')
addpath('/Users/brainmodulationlab/git/bml')
ft_defaults
bml_defaults



%% save a few files of each patient locally

% clc
% %save electrode table
% for pt = 3023
%     clear filenamee foldernamee1 foldernamee2
%     foldernamee1 = sprintf('/Volumes/Nexus/DBS/DBS%4.0f/Preprocessed Data/Sync/annot/DBS%4.0f_electrode.txt',pt,pt);
%     
%     %foldernamee2 = sprintf('/Users/ajorge/ajorge/data/ecog_stnpathways/electrode_annot/');
%     
%     foldernamee2 = sprintf('/Users/ajorge/Desktop/w.txt')
%     cd(foldernamee2)
%     copyfile('foldernamee1.txt', 'foldernamee2.txt')
%     
% end
% 




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define pt

clear all
close all
clc

pt = 3030 %[3016 3023 3024] %[ 3019 3020 3022 3023 3024 3025 3026 3027 3028 3029 3030 3031]   %[3005 3006 3008 3010 3011 3012 3016 3017 3018 3019 3020 3022 3023 3024 3025 3026 3027 3028 3029 3030 3031 3032]                %3010 is 1000 seconds long. crashes


% load raw data at 1kHz to identify sets
clear raw_all address
cfg4 = [];
cfg4.pt = pt;
cfg4.desired_f_sample = 1;
cfg4.plot_databrowser = 0;
[raw_all, output] = ecog_getraw(cfg4);   

infoo_ns2 = output.infoo_ns2;

% mark stim trials (rough approximation of stim time done at 1kHz)
clear stim_marker_1khz stim_marker_time_1khz
cfg4 =[];
cfg4.pt = pt;
cfg4.sett = 999;
cfg4.flag = 1;
cfg4.fileidd = output.fileidd;
[stim_marker_1khz, stim_marker_time_1khz] = define_trials_ecog(cfg4,raw_all);


%% load raw data at 30kHz set by set

for sett =  [1 2 3 4 5 6 7 8 9 10 11 12 25]% [3 4]%[1:6]%[1:max( stim_marker_1khz(find(stim_marker_1khz(:,2)<25),2)) 25]
    
    
    
    %% load raw data at 30kHz set by set (30kHz files are big, load one set at a time)
    clear raw address
    cfg4 = [];
    cfg4.pt = pt;
    cfg4.desired_f_sample = 30;
    indx = find(stim_marker_1khz(:,2)==sett);
    cfg4.stim_marker_time_1khz = stim_marker_time_1khz(indx,1);
    cfg4.infoo_ns2 = infoo_ns2;
    [raw, output] = ecog_getraw(cfg4);   
    figure; plot(raw.time{1}, raw.trial{1}(26,:),'k'); hold on
    
    %% mark stim trials (finer approximation of stim time done at 30kHz)
    clear stim_marker2 stim_marker_time
    close all
    cfg4 = [];
    cfg4.pt = pt;
    cfg4.sett = sett;
    cfg4.flag = 30;
    [stim_marker_30khz, stim_marker_time_30khz] = define_trials_ecog(cfg4,raw);
    stim_marker_30khz(:,1) = stim_marker_30khz(:,1) + output.s1;
    stim_marker_30khz(:,2) = sett;
    stim_marker_time_30khz(:,2) = sett; 
    
    %% save values
    cd('/Users/ajorge/ajorge/data/ecog_stnpathways/stim_markers/')
    load('otherdata','otherdata');    
    otherdata{pt}.stim_marker_30khz{sett}       = stim_marker_30khz;
    otherdata{pt}.stim_marker_time_30khz{sett}  = stim_marker_time_30khz;
    otherdata{pt}.ptaddress{sett}               = output.path2pt;
    otherdata{pt}.datasamplerate{sett}          = raw.fsample;
    save('otherdata','otherdata');
    clear raw
    
    %% obtain processed data (trial based using stim_marker2 info)
    clear pro
    cfg = [];
    cfg.pt = pt;
    cfg.stim_marker_time_30khz = otherdata{pt}.stim_marker_time_30khz{sett};
    cfg.stim_marker_30khz = otherdata{pt}.stim_marker_30khz{sett} ;
    cfg.s1 = output.s1;
    cfg.address = otherdata{pt}.ptaddress{sett};
    cfg.trialflag = 10;
    cfg.fsample = 30000;
    cfg.sett = sett;

    pro = ecog_raw2pro(cfg);
    
    %deal with fragmented ns5 files for 3016, 3028, 3029
    if cfg4.pt == 3016 && output.fileidd == 4
        pro.trialinfo = pro.trialinfo+3;
        sett_save = sett+3;
    elseif cfg4.pt == 3028 && output.fileidd == 3
        pro.trialinfo = pro.trialinfo+5;
        sett_save = 6;
    elseif cfg4.pt == 3028 && output.fileidd == 4
        pro.trialinfo = pro.trialinfo+6;
        sett_save = sett+6;
    elseif cfg4.pt == 3029 && output.fileidd == 3
        pro.trialinfo = pro.trialinfo+6;
        sett_save = sett+6;    
    else
        sett_save = sett;
    end
    
    

    cd /Users/ajorge/ajorge/data/ecog_stnpathways/pro/
    %filenamee = sprintf('pro_DBS%4.0f_sett%02.0f',pt,sett);
    filenamee = sprintf('pro_DBS%4.0f_sett%02.0f',pt,sett_save);
    save(filenamee,'pro')

    

    

end









%% take out "bad" trials (trials with high variance)
cfg          = []; 
%cfg.viewmode = 'vertical';
cfg.layout    	= 'bmlecogarray.mat';
pro = ft_rejectvisual(cfg, pro);
















