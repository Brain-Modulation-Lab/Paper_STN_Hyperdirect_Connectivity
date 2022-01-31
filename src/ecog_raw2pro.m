function pro = ecog_raw2pro(cfg2)


% example:
% cfg = [];
% cfg.pt = pt;
% cfg.stim_marker2 = stim_marker2;
% cfg.address = address;
% cfg.trialflag = 1;



% trialflag == 1
% process ecog data into trials
% it marks as a trial one single stim (usually 1 second)
% it marks as a set one complete set of stim (usually 30 seconds)

% trialflag==2
% process data into sets
% set0 = baseline period
% set1 = first stim set (30 seconds worth of data)
% ...

% trialflag == 5
% # divide into trials (eg. -250 to 1250 m) 
% # preprocess into trials
% # blank the stimulation period (eg. estimated between 0.005 to 0.010 ms)
% # splines this blanking distance
% # filters line noise, high pass and low pass filters
% # rereferences
% # re preprocess with new filters and re-references



% ajorge 2018




%% translate basic values

flag            = cfg2.trialflag;
if flag == 1 || flag == 2 || flag == 3 || flag == 5
    stim_marker2    = cfg2.stim_marker2;
elseif flag == 10
    stim_marker2    = cfg2.stim_marker_30khz;
end

pt              = cfg2.pt;
fsample         = cfg2.fsample;
%--------------------------------------------------------------------------
 




%%

%start building the configuration cfg
cfg = [];
cfg.datafile = cfg2.address;



%--------------------------------------------------------------------------
% channel selection

if fsample == 1000
    %select anything with the word lfp
    cfg.channel       = 'lfp*';
elseif fsample == 30000
    %select anything with the word raw
    cfg.channel       = 'raw*';
else
    warning('unable to select channels')
end

% another way to select channels individually
%{
for i = [1:63]
    mat4{i,1} = sprintf('raw %1.0f',i);
end
cfg.channel = mat4;
%}
%--------------------------------------------------------------------------




%--------------------------------------------------------------------------
% patient specific baseline time periods
% baselines are different for every patient (lenght and quality of data). manually select these settings  

if     pt == 3008
    basei = 2;
    basef = 30;
    
elseif pt == 3011
    basei = 2;
    basef = 20;%floor(stim_marker2(1,1)/1000);
    
elseif pt==3017
    % i looked into the baseline freq-time plot, it looks like there is
    % some weird noise between 1 and 12, and 35:45. For now, eliminate
    % those time periods from baseline. 12:35 looks homogeneous.
    basei = 12;
    basef = 35; %floor(stim_marker2(1,1)/1000);  %the first stim bin   
    
elseif pt==3018
    basei = 2;
    basef = floor(stim_marker2(1,1)/fsample);  %the first stim bin  
    
else
    warning('need to add patient-specific baseline time period')
    basei = 2;  %start at the second "trial" to be able to look backwards
    basef = 30;
end

%--------------------------------------------------------------------------












%--------------------------------------------------------------------------
% apply filters (60Hz, high pass, low pass, etc)
if flag == 1 || flag == 2 || flag == 3 || flag == 10

    %filter out line noise 
    %cfg.bsfilter = 'yes';
    %cfg.bsfreq = [59 61; 119 121; 179 181];
    
    %filter out low frequencies
    %cfg.hpfilter = 'yes';
    %cfg.hpfreq = 0.01; 
    
    %cfg.lpfilter = 'yes';
    %cfg.lpfreq = 500;
    

   
    %{
    %high pass filter
    try
        cfg.hpfreq = cfg.hpfreq; 
        cfg.hpfilter = 'yes';
       % cfg.hpfreq = 8; 
    catch
    end
    
    %low pass filter
    %if exist('cfg.lpfreq') == 1
    try
        cfg.lpfreq = cfg.lpfreq; 
        cfg.lpfilter = 'yes';
        %cfg.lpfreq = 12; 
    catch
    end
    %}
    
    
    %this patient has a big noise with this particular frequency
    %{
    if pt == 3020
        cfg.dftfilter = 'yes';
        cfg.dftfreq   = [182.7037]*(1:60);
        cfg.dftreplace = 'neighbour';
        
        cfg.dftbandwidth = [8]*ones(1,60);
        cfg.dftneighbourwidth = [16]*ones(1,60);
        
        warning('this setting is on')
        pause

    end
    %}
    
    
    
end
%--------------------------------------------------------------------------






%--------------------------------------------------------------------------
% re-referecing
if flag == 1 || flag == 2 || flag == 3 
    %[cfg.channel]     = channelselection_preprocessing(dummyraw);
    
    clear mat4 cfg.channel
    j=1;
    
    
    %---
    %rereferncing using all channels except the last few rows on the array
    %{
    for i = [1:15 22:36 43:57]
        mat4{j,1} = sprintf('lfp %1.0f',i);
        j=j+1;
    end
    %}
    %---
    %rereferncing using all channels
    %{
    for i = [1:63]
        mat4{j,1} = sprintf('lfp %1.0f',i);
        j=j+1;
    end
    %}
    
    %---
    
    %rereferncing using selected channels (a group, say, on the motor cortex)
    
    %{
    for i = [8 9 10     29 30 31    50 51 52]
        mat4{j,1} = sprintf('raw %1.0f',i);
        j=j+1;
    end
    %}
    
    %{
    for i = [1:63 129:192]
        mat4{j,1} = sprintf('raw %1.0f',i);
        j=j+1;
    end
    %}
    
    %{
    cfg.reref         = 'yes';      %'no' or 'yes' (default = 'no')
    cfg.refchannel    = mat4;      %cell-array with new EEG reference channel(s), this can be 'all' for a common average reference
    %cfg.refchannel    = 'all';      %cell-array with new EEG reference channel(s), this can be 'all' for a common average reference
    cfg.refmethod     = 'median';   %'avg', 'median', or 'bipolar' for bipolar derivation of sequential channels (default = 'avg')
    %}
    
    
    
    %---
    
    %bipolar montage 
    
    %{
    % horizontal pairs (1 and 22, 22 and 43 ...)
    for p = 1:42
        mont(p,p)       = 1;
        mont(p,p+21)    =-1;
    end
    for i = [1:63]
        oldchans{i,1} = sprintf('raw %1.0f',i);
    end
    for i = [1:42]
        newchans{i,1} = sprintf('bipol %1.0f',i);
    end
    %}
    
    %{
    % vertical pairs (1 and 2, 3 and 4 ...)
    for p = 1:20
        mont(p,p) = 1;
        mont(p,p+1)=-1;
    end
    for i = [1:63]
        oldchans{i,1} = sprintf('raw %1.0f',i);
    end
    for i = [1:20]
        newchans{i,1} = sprintf('bipol %1.0f',i);
    end
    
    
    cfg.montage.tra       = mont;
    cfg.montage.labelold  = oldchans;
    cfg.montage.labelnew  = newchans;
    
    %cfg.reref         = 'yes';      %'no' or 'yes' (default = 'no')
    %cfg.refchannel    = 'all';      %cell-array with new EEG reference channel(s), this can be 'all' for a common average reference
    cfg.refmethod     = 'bipolar';   %'avg', 'median', or 'bipolar' for bipolar derivation of sequential channels (default = 'avg')
    %}
    
    
    
    
    
    
    
    %cfg.demean = 'yes';
    %cfg.baselinewindow = [-0.0015 -0.0005];        %seconds
    %cfg.detrend = 'yes';
    
    
   
    
    
    
end
%--------------------------------------------------------------------------











%--------------------------------------------------------------------------
% divide into trials
if flag==1 || flag==4 || flag==5
    
    offset = 250*fsample/1000; %how many ms/bins to collect before stimulation
    %offset = -3*fsample/1000;
    
    %--------------------------------------------------------------------------
    % add baseline trials
    % cut at x seconds before the first stim trial occurs and split that into 1000 ms trials
    clear matt
    i56=1;
    for i55 = basei:basef
        matt(i56,1) = fsample*(i55-1) + 1 - offset;
        matt(i56,2) = fsample*(i55)+offset;
        matt(i56,3) = -offset;              %offset;
        matt(i56,4) = 0;                    %mark it as set0
        i56=i56+1;
    end

    %--------------------------------------------------------------------------
    % add stim trials based on stim_markers
    clear matt2
    
        %---standard fieldtrip definitios for trl---
        %{
        % define the trials
        trl(:,1) = stimulus_sample + pretrig;  % start of segment
        trl(:,2) = stimulus_sample + posttrig; % end of segment
        trl(:,3) = pretrig;                    % how many samples prestimulus

        % add the other information
        % these columns will be represented after ft_preprocessing in "data.trialinfo"
        % the last column in this example contains a "correct" boolean flag for each trial
        trl(:,4) = stimulus_value;
        trl(:,5) = response_value;
        trl(:,6) = reaction_time;
        %}
        %---
    
    matt2 = [stim_marker2(:,1)-offset   stim_marker2(:,1)+fsample*ones(size(stim_marker2,1),1)+offset    0*ones(size(stim_marker2,1),1)-offset     stim_marker2(:,2) ];
    %matt2 = [stim_marker2(:,1)-offset   stim_marker2(:,1)+fsample*ones(size(stim_marker2,1),1)/2+offset 0*ones(size(stim_marker2,1),1) stim_marker2(:,2)];
    %--------------------------------------------------------------------------

    cfg.trl = [matt; matt2];
    
    
    %{
    % THIS just didn't work. could try again later
    % cfg - configuraton structure
    % cfg.epoch - ANNOT table with new epoching
    % cfg.t0 - reference time for each epoch. If not given the time is kept in
    %          global time coordinates. Can be string or char is given that matches a 
    %          variable of cfg.epoch or a numeric vector of the same length than cfg.epoch
    % cfg.regularize - if true, resulting times are forced to be equal. 
    % cfg.warn - logical indicating if warnings should be issued. Defaults to true
    %
    % returns a raw with new trials. The epoch ANNOT is added as a new field in the 
    % raw, changing the id to match the index of the corresponding trials if necessary.
    cfg2 = [];
    cfg2.epoch = [(1:length(cfg.trl))' cfg.trl(:,1) cfg.trl(:,2) cfg.trl(:,2)-cfg.trl(:,1)];
    %pro2 = bml_redefinetrial(cfg2,raw) 
    %pro4 = bml_load_epoched(cfg)
    %}
    
end
%--------------------------------------------------------------------------   
    
    




    
%--------------------------------------------------------------------------   
% divide into sets    
if flag == 2
    
    padding = 100;
    
    %define baseline set
    clear matt
    matt(1,1)=basei*fsample;
    matt(1,2)=basef*fsample+padding;
    matt(1,3)=0;
    matt(1,4)=0;
    
    %define stim sets
    clear matt2
    for i6 = 1:stim_marker2(end,2)
        
        inx = find(stim_marker2(:,2)==i6);
        
        matt2(i6,1)=stim_marker2(inx(1),1)-padding;
        matt2(i6,2)=stim_marker2(inx(end),1)+fsample+padding;
        matt2(i6,3)=0;
        matt2(i6,4)=i6;
    end
    
    cfg.trl = [matt; matt2];
      
end 
%--------------------------------------------------------------------------
    









%--------------------------------------------------------------------------
% divide into trials (centered at 0, from -??? ms to ??? ms)
if flag==10
    
    %select how many seconds to collect before and after stimulation
    secs_i = 0.010;
    secs_f = 0.100;
    
    bins_i = secs_i*fsample;
    bins_f = secs_f*fsample;    
      
        %---standard fieldtrip definitios for trl---
        %{
        % define the trials
        trl(:,1) = stimulus_sample + pretrig;  % start of segment
        trl(:,2) = stimulus_sample + posttrig; % end of segment
        trl(:,3) = pretrig;                    % how many samples prestimulus

        % add the other information
        % these columns will be represented after ft_preprocessing in "data.trialinfo"
        % the last column in this example contains a "correct" boolean flag for each trial
        trl(:,4) = stimulus_value;
        trl(:,5) = response_value;
        trl(:,6) = reaction_time;
        %}
        %---
    
    %--------------------------------------------------------------------------
    % add baseline trials
    % cut at x seconds before the first stim trial occurs and split that into 1000 ms trials
    clear matt
    i56=1;
    for i55 = basei:basef
        matt(i56,1) = fsample*(i55-1) + 1 - bins_i;
        matt(i56,2) = fsample*(i55)+bins_i;
        matt(i56,3) = -bins_i;              %offset;
        matt(i56,4) = 0;                    %mark it as set0
        i56=i56+1;
    end

    %--------------------------------------------------------------------------
    % add stim trials based on stim_markers
    clear matt2
    matt2 = [stim_marker2(:,1)-bins_i      stim_marker2(:,1)+bins_f    0*ones(size(stim_marker2,1),1)-bins_i     stim_marker2(:,2) ];
    %--------------------------------------------------------------------------
    
    
    %cfg.trl = [matt; matt2];
    cfg.trl = matt2;
    
   
end
%--------------------------------------------------------------------------   





    
    






%% --------------------------------------------------------------------------
% do actual preprocessing
if flag == 4
    % for this one, we want to redefine trials according to cfg.trl
    cfgx=[];
    cfgx.trl = cfg.trl;
    data = cfg.data;
    pro = ft_redefinetrial(cfgx,data);
end
    
if flag==1 || flag==2 || flag ==3 || flag==5 || flag==10
    %here we want to preproces with filtering and CAR options
    pro = ft_preprocessing(cfg);
    
end










%% --------------------------------------------------------------------------
% 
% blank out stimulation and interpolate between points
pro2=pro;
if flag == 5
    
    blankingrange = [-0.003 0.005];
    
    for trial = 1:length(pro.time)
        
        for chan = 1:length(pro.label)
            bini = find(pro.time{trial} == blankingrange(1));
            binf = find(pro.time{trial} == blankingrange(2));

            v(:,1) = pro.trial{trial}(chan,bini-1);
            v(:,2) = pro.trial{trial}(chan,bini);
            v(:,3) = pro.trial{trial}(chan,binf);
            v(:,4) = pro.trial{trial}(chan,binf+1);
            
            x(:,1) = pro.time{trial}(bini-1);%*ones( size(pro.trial{trial},1),1);
            x(:,2) = pro.time{trial}(bini);%*ones( size(pro.trial{trial},1),1);
            x(:,3) = pro.time{trial}(binf);%*ones( size(pro.trial{trial},1),1);
            x(:,4) = pro.time{trial}(binf+1);%*ones( size(pro.trial{trial},1),1);
  
            xq = (blankingrange(1):1/fsample:blankingrange(2));%.*ones( size(pro.trial{trial},1),1);
            vq = interp1(x,v,xq,'spline');
            
            pro2.trial{trial}(chan,bini:binf) = vq;
            
        end
        
        
    end
    
    %dumb check
    %{
    figure(56)
    chan = 5;
    trial = 50;
    plot(pro.time{trial},pro.trial{trial}(chan,:)); hold on
    plot(pro2.time{trial},pro2.trial{trial}(chan,:)); hold on
    %}
    
    
    
    pro=pro2;
    
    %{
    cfg = [];
    %filter out line noise 
    cfg.bsfilter = 'yes';
    cfg.bsfreq = [59 61; 119 121; 179 181];
    
    %filter out low frequencies
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 12; 
    
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 25;
    
    pro3 = ft_preprocessing(cfg, pro2);
    %}
    
    %{
    figure(56)
    chan = 5;
    trial = 50;
    plot(pro3.time{trial},pro3.trial{trial}(chan,:)); hold on
    
    pro=pro3;
    %}
   
end












%%
%--------------------------------------------------------------------------
% if there are analog channel, need to eliminate last channels (they are analog channels)

if contains(pro.label{1},'analog') == 1 
    warning('There are analog channels, make sure we are dealing with them appropriately with addition of electrode table lookup')
    %pause
end



if size(pro.label,1) ~= size(pro.trial{1},1)
    warning('Number of labels ~= number of channels on pro.trials - attempting to manually change labels')

    % how many electrodes over?
    bigs   = size(pro.trial{1},1);
    smalls = size(pro.label,1);
    sizediff = bigs - smalls;

    %erase last sizediff electrodes
    for t = 1:length(pro.trial)
        pro.trial{t} = pro.trial{t}(1:bigs-sizediff,:);
    end


    %if size continues to be an issue, stop!
    if size(pro.label,1) ~= size(pro.trial{1},1)
        warning('Number of labels ~= number of channels on pro.trials - will cause problems later on (need to be equal)')
        %pause;
    else
        fprintf('Attempt successful \n')
    end


end



%% look at electrode table

zero_ch1to4=0;

%ideally, the first id electrode should always be ecog_101
try
    foldernamee1 = sprintf('/Volumes/Nexus/DBS/DBS%4.0f/Preprocessed Data/Sync/annot/DBS%4.0f_electrode.txt',pt,pt);
    T = readtable(foldernamee1);

    %go through each recording in pro
    for rawx = 1:size(pro.trial{1},1)

        
        %by default assign raw_i to raw_i
        for trialx = 1:size(pro.trial,2)
            pro.trial2{trialx}(rawx,:) = pro.trial{trialx}(rawx,:); 
        end
        

        %if cases where some ecogs were reversed
        try
            rawxname = sprintf('raw %1.0f',rawx);
            indx1 = find(strcmp(T.channel,rawxname));
            indx2 = find(strcmp(pro.label,rawxname));

            %switch channels
            if indx1 ~= indx2
                zero_ch1to4=1;
                for trialx = 1:size(pro.trial,2)
                    pro.trial2{trialx}(indx2,:) = pro.trial{trialx}(indx1,:); 
                end
            end
        catch
        end

    end
    
    
    
    %if we detected zero_ch1to4, then completely zero out channels 1:4
    if zero_ch1to4 == 1 && (pt>3020 && pt<3028)
        for indx2 = 1:4
            for trialx = 1:size(pro.trial,2)
                pro.trial2{trialx}(indx2,:) = 0*pro.trial{trialx}(indx2,:); 
            end
        end
    elseif zero_ch1to4 == 1
        fprintf('this patient shouldnt have ch_1:4 zeroed!!! look into it\n');
        pause
    end

    pro.trial=[];
    pro.trial = pro.trial2;
    pro.trial2=[];

    %rename "raw xx" labels with "ecog101" labeling system
    pro.label = T.electrode;
    
catch
    
end



%% --------------------------------------------------------------------------
% plot same figures that were plotted with raw to make sure channels are linning up

figure(45)
counter57=1;

indx = find(pro.trialinfo==cfg2.sett);

for channn = [1 26 52]
    subplot(3,1,counter57)
    tt = sprintf('channel %3.0f',channn);
    counter57=counter57+1;
    plot(pro.time{indx(1)}(:),  pro.trial{indx(1)}(channn,:)); hold on
    title(tt)
    
    
end




