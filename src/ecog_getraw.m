function [raw, output] = ecog_getraw(cfg4)
%
% function to get raw data for the stim experiments from ECOG files
%
% either obtain the 1 kHz data
% or 30 kHz data
%
% ajorge jun 2018






%% which file id?

if cfg4.pt == 3004
    fileidd = 1;
    
elseif cfg4.pt == 3006
    fileidd = 1;
    
elseif cfg4.pt == 3008
    fileidd = 3;
    
elseif cfg4.pt == 3010
    fileidd = 3;
elseif cfg4.pt == 3011
    fileidd = 3;
    
elseif cfg4.pt == 3012 && cfg4.desired_f_sample == 1
    fileidd = 2; 
elseif cfg4.pt == 3012 && cfg4.desired_f_sample == 30
    fileidd = 1;
    
elseif cfg4.pt == 3016
    fileidd = 4;  % 3 and 4, but selected 4 only (set 4 5 6)
    
    
elseif cfg4.pt == 3025
    fileidd = 3;
elseif cfg4.pt == 3026
    fileidd = 3;
  
   
    
elseif cfg4.pt == 3028
    warning('pt 3028 is special, make sure pointing to the correct filedd')    
    sprintf('set 1-5  fileidd = 2\nset 6    fileidd = 3\nset 7-12 fileidd = 4\n')
    prompt = 'Desired fileidd =  ';
    fileidd = input(prompt);        
    
elseif cfg4.pt == 3029
    %{
    warning('pt 3029 is special, make sure pointing to the correct filedd')    
    sprintf('set 1-6  fileidd = 2\nset 7-12 fileidd = 3\n')
    prompt = 'Desired fileidd =  ';
    fileidd = input(prompt);     
    %}
    fileidd=3
    
elseif cfg4.pt == 3031
    fileidd = 7;
elseif cfg4.pt == 3032
    fileidd = 2;
    
    
elseif cfg4.pt == 4068
    fileidd = 4;   
    
elseif cfg4.pt == 4079
    fileidd = 3; 
    
elseif cfg4.pt == 4081
    fileidd = 3;
    
elseif cfg4.pt == 4088
    fileidd = 4;
    
else
    %defaults to:
    fileidd = 2;
    warning('need to check that this fileidd is correct at some point \n')
end


output.fileidd = fileidd;


%% loading 1 kHz data

if cfg4.desired_f_sample == 1
    

    SUBJECT=sprintf('DBS%4.0f',cfg4.pt);
    %SUBJECT='DBS3010';

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
    cfg.roi = info_ns2(fileidd,:);
    %[raw, file_raw_map] = bml_load_continuous(cfg);
    [raw, ~] = bml_load_continuous(cfg);

    %path2pt
    output.path2pt = sprintf('%s/%s',cfg.path, string(cfg.roi.name));

    %give time (seconds) when ns2 started
    output.infoo_ns2 = cfg.roi;
    
    %view raw data
    try
    if cfg4.plot_databrowser == 1
        cfg2=[];
        cfg2.viewmode = 'vertical';
        %cfg2.blocksize = 'all';
        %cfg2.channel = 'all'
        ft_databrowser(cfg2,raw)
    end
    catch
    end


    %{
    %fieldtrip automatic artifact detection
    cfg.continuous='yes';
    cfg.artfctdef.zvalue.channel='all';
    cfg.artfctdef.zvalue.cutoff=20;
    cfg.artfctdef.zvalue.interactive = 'yes';
    cfg.artfctdef.zvalue.bpfreq      = [0.9 1.1];
    [cfg, artifact_muscle] = ft_artifact_zvalue(cfg,raw);
    ft_rejectartifact(cfg,raw);
    %}

    
    %we'll load the 30kHz data for these particular patients (they were
    %only stimulated at 1mA and is better appreciated at 30khz). So why not
    %load the 30kHz for all patients? it takes a long time and takes lots
    %of RAM (computer migth simply crash, and we have 32GB of ram memory!!!)
    if cfg4.pt==3004 || cfg4.pt==3005 || cfg4.pt==3006 || cfg4.pt==3008 || cfg4.pt==3011 || cfg4.pt == 4079
        SUBJECT=sprintf('DBS%4.0f',cfg4.pt);
        PATH_DATA='/Volumes/Nexus/DBS';
        PATH_TRELLIS=[PATH_DATA filesep SUBJECT filesep 'Raw' filesep 'Trellis'];
        cfg=[];
        cfg.path=PATH_TRELLIS;
        cfg.pattern='*.ns5';                %could be either .ns2=1kHz or .ns5=30kHz
        cfg.filetype='trellis.ns5';         %could be either .ns2=1kHz or .ns5=30KHz
        info_ns5    = bml_info_raw(cfg);
        info_ns5    = bml_roi_table(info_ns5,'trellis');
        cfg.roi     = info_ns5(fileidd,:);
        cfg.chantype = 'raw';
        
        %cfg.channel = {'raw 1','raw 5'};
        
        [raw, ~] = bml_load_continuous(cfg);
    end






    
    
    
   

elseif cfg4.desired_f_sample == 30
%% loading 30 kHz data
    % 30 kHz data is too big to load at once.
    % load 1kHz data first, identify ROIs, then load 30 kHz data (set by set)


    SUBJECT=sprintf('DBS%4.0f',cfg4.pt);
    PATH_DATA='/Volumes/Nexus/DBS';
    PATH_TRELLIS=[PATH_DATA filesep SUBJECT filesep 'Raw' filesep 'Trellis'];

    cfg=[];
    cfg.path=PATH_TRELLIS;
    cfg.pattern='*.ns5';                %could be either .ns2=1kHz or .ns5=30kHz
    cfg.filetype='trellis.ns5';         %could be either .ns2=1kHz or .ns5=30KHz
    
    %patch this error        
    if cfg4.pt == 3012 || cfg4.pt == 4079
        
    else
        cfg.chantype = 'raw';
    end
    
    %for i = [1:63]
    %    mat4{i,1} = sprintf('raw %1.0f',i);
    %end
    %cfg.channel = mat4;
    
    
    clear info_ns5
    info_ns5    = bml_info_raw(cfg);
    info_ns5    = bml_roi_table(info_ns5,'trellis');
    cfg.roi     = info_ns5(fileidd,:);
    %cfg.chantype = 'raw';
    
    
    %check to see that ns2_info start time is the same as ns5_info start time
    if cfg4.pt == 3010  || cfg4.pt == 3012
        
        delta_ns2ns5 = cfg4.infoo_ns2.starts - cfg.roi.starts;
        
        cfg.roi.starts = cfg.roi.starts + delta_ns2ns5;
        cfg.roi.ends   = cfg.roi.ends   + delta_ns2ns5;
        cfg.roi.t1     = cfg.roi.t1     + delta_ns2ns5;
        cfg.roi.t2     = cfg.roi.t2     + delta_ns2ns5;
        
    elseif cfg4.pt == 3005 || cfg4.pt == 3006 || cfg4.pt == 3008 || cfg4.pt == 3011
        %ignore since we are reading the 30kHz file to obtain initial
        %markers anyway
        
    elseif cfg4.infoo_ns2.starts ~= cfg.roi.starts
        cfg4.infoo_ns2.starts
        cfg.roi.starts
        warning('ns2 ns5 starting times are not equal')
        pause
    end

      
    
    %redefine ROI based on sett values
    filestart2setstart = cfg4.stim_marker_time_1khz(1)   - cfg.roi.starts         - 1;
    filestart2setend   = cfg4.stim_marker_time_1khz(end) - cfg.roi.starts         + 2;



    cfg.roi.starts  = cfg4.stim_marker_time_1khz(1)     - 1;    %substract 1 second for padding
    cfg.roi.ends    = cfg4.stim_marker_time_1khz(end)   + 2;    %add 1 second for padding 
    cfg.roi.duration    = cfg.roi.ends - cfg.roi.starts;

    cfg.roi.t1      = cfg.roi.starts;
    cfg.roi.t2      = cfg.roi.ends;

     
    
    
    cfg.roi.s1      = filestart2setstart*30000;
    cfg.roi.s2      = filestart2setend*30000;
    
    
    
    %[raw, file_raw_map] = bml_load_continuous(cfg);
    cfg.chantype = 'raw';
    [raw, ~] = bml_load_continuous(cfg);

    %path2pt
    output.path2pt = sprintf('%s/%s',cfg.path, string(cfg.roi.name));
    
    output.s1       = cfg.roi.s1;
    output.s2       = cfg.roi.s2;

end


