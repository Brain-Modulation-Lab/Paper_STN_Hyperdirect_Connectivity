function [stim_marker2, stim_marker_time] = define_trials_ecog(cfg4,raw)

% function that will identify when stimulus occured by visual inspection 
% (i.e. when did the voltage increased suddenly?, implying that that is when stim happened)

% input = raw data, pt number
% output = stim_marker2 = bin when stim occured, stim_marker_time = time in
% ms when stim occured


% raw comes from:
% cfg4=[];
% cfg4.path=PATH_TRELLIS;
% cfg4.pattern='*.ns2';
% cfg4.filetype='trellis.ns2';
% info_ns2=bml_info_raw(cfg4);
% info_ns2 = bml_roi_table(info_ns2,'trellis');
% cfg4.chantype = 'lfp';
% cfg4.roi = info_ns2(fileidd,:);
% [raw, file_raw_map] = bml_load_continuous(cfg4);


% ajorge june 2018





%%
pt = cfg4.pt;

% identify channels that showcase a pronounced response to stimulation
% aka, which channel has the highest median value?
for ch = 1:63
    voltage_median(ch) = abs(median(raw.trial{1}(ch,:)));
end
[ ~ , indx] = max(voltage_median);
[sorted , indxs] = sort(voltage_median);
indx_median = indxs(30);


% patch a known "clean" electrode that shows clear stimulations
switch pt 
    case 3005
        indx = 60;
        
    case 3006
        indx = 60;%32+(161-129)+(288-258)
        
    case 3008
        indx = 5;
    case 3010
        indx = 85;
    case 3011
        indx = 20; %this channel was found to be a good one for identifying artifacts
    case 3012
        indx = 21;
    case 3018
        indx = 20;
    case 3022
        indx = 106; %this channel was found to be a good one for identifying artifacts
    case 3023
        indx = 106;
    case 3024
        indx = 17;
    case 3026
        indx = 106;
    case 3027
        indx = 81;
    case 3028
        indx = 85;
    case 3029
        %indx = find(strcmp(raw.label,'raw 137'));
        indx = 105;
        
    case 3030
        indx = 45;
    case 3031
        indx = 149;
    case 3032
        indx = 85;
        
    case 4068
        indx = 21;
        
    case 4079
        indx = 42;
        
    case 4088
        indx = 103;
        
end



% use this channel to identify places where voltage value changes
% dramatically with the idea that a high voltage delta in a short time is
% equivalent to a stimulation (double check this with plots)
if cfg4.sett == 999
    figure(22)
else
    figure(33)
end
subplot(2,1,2)
delta_lim = 1000;

%patch this pt
if pt == 3029
    delta_lim = 500;
end

plot(raw.time{1}, raw.trial{1}(indx,:),'k'); hold on
%plot(ones(length(raw.trial{1}(indx,:)),1)*delta_lim,'k'); hold on

subplot(2,1,1)
plot(raw.time{1}, raw.trial{1}(indx_median,:),'k'); hold on
%plot(ones(length(raw.trial{1}(indx_median,:)),1)*delta_lim,'k'); hold on

set(gcf,'Position',[ 111         676        1921         478])




%%

if cfg4.sett ~=25
    m=1; %lets define a trial as one single stim pulse
    sett=1; %lets define a sett as multiple trials of stim (each stim is a one 1hz pulse)
    delta_time_lim = 2500*raw.fsample/1000; %  1.5 seconds 
    
    %patch this pt
    if pt == 3005
        delta_time_lim = 1500*raw.fsample/1000;
    elseif pt == 3031
        delta_time_lim = 250*raw.fsample/1000;
    end

    clear stim_marker
    for i = 2:length(raw.trial{1}(indx,:))
        delta = abs(raw.trial{1}(indx,i)-raw.trial{1}(indx,i-1));
        if delta > delta_lim
            stim_marker(m,1) = i;
            stim_marker(m,2) = sett;        
            m=m+1;
        end
    end



    %eliminate markers that are too close to each other (closer than 50ms)
    j=1;
    stim_marker2=[];
    stim_marker2(j,1) = stim_marker(1,1);
    for i = 2:length(stim_marker)
        delta = abs(stim_marker(i,1) - stim_marker(i-1,1));
        if delta > 50*raw.fsample/1000
            j=j+1;
            stim_marker2(j,1) = stim_marker(i,1);
        end
    end

    %mark setts
    stim_marker2(1,2) = 1; 
    for i = 2:length(stim_marker2)
        delta_time = abs(stim_marker2(i,1) - stim_marker2(i-1,1));

        %if time between markers is bigger than 2.5 seconds, is a new set
        if delta_time > delta_time_lim
            sett=sett+1;
        end
        stim_marker2(i,2) = sett;  
    end


    % automatically check if  stim markers spaced every 1 second?
    % come back to this later


    %translate stim markers to time
    for i = 1:length(stim_marker2)
        stim_marker_time(i,1) = raw.time{1}(stim_marker2(i,1));
        %stim_marker_time(i,2) = stim_marker2(i,2);
        stim_marker_time(i,2) = cfg4.sett;
    end



    %% if sett has less than 10 trials, is not a set, eliminate and reorganzie


    if  cfg4.flag == 1
    %get rid of set 1,4,6,7,8,9,10,11

        %identify sets with less than 30 trials
        k10=1; k11=1; badsets=[];
        for k9 = 1:max(stim_marker2(:,2))
            dummyx = find(stim_marker2(:,2)==k9);
            if length(dummyx)<20
                badsets(k10) = k9; k10=k10+1;
            else
                godsets(k11) = k9; k11=k11+1;
            end
        end

        %manual patch
        if pt == 3024
            badsets = [  1     4     6     7     8      9    10];
            godsets = [  2     3     5     11    12    13];
        end

        if isempty(badsets) == 0
            for jj = badsets
                dummyindx = find(stim_marker2(:,2)==jj);
                stim_marker2(dummyindx,:) = [];
                stim_marker_time(dummyindx,:) = [];
            end

            c5=1;
            for jj = godsets
                dummyindx = find(stim_marker2(:,2)==jj);
                stim_marker2(dummyindx,2) = c5;
                stim_marker_time(dummyindx,2) = c5;
                c5=c5+1; 
            end
        end
    end
    %}




    %--------------------------
    % special case
    %{
    if pt == 3010
        temp_stim_marker2 = stim_marker2;
        temp_stim_marker_time = stim_marker_time;

        needed = 0;
        for jj = 1:max(stim_marker2(:,2))
            dummyindx = find(stim_marker2(:,2)==jj);
            if length(dummyindx)<=10
                temp_stim_marker2(dummyindx,:) = nan;
                needed = 1;
            end
        end  

        if needed == 1
            indx1 = isnan(temp_stim_marker2(:,2));
            indx2 = find(indx1==1);
            stim_marker2(indx2,:) = [];
            stim_marker_time(indx2,:) = [];
            stim_marker2(:,2) = stim_marker2(:,2) - min(stim_marker2(:,2)) + 1;
        end
    end
    %}

    %--------------------------
    % special case for dbs3024
    %{
    if pt == 3024 && cfg4.flag == 1
        warning('do we need this special case anymore?')
        pause
    %get rid of set 1,4,6,7,8,9,10,11   
        for jj = [1,4,6,7,8,9,10,11]
            dummyindx = find(stim_marker2(:,2)==jj);
            stim_marker2(dummyindx,:) = [];
            stim_marker_time(dummyindx,:) = [];
        end

        c5=1;
        for jj = [2 3 5 12 13 14]
            dummyindx = find(stim_marker2(:,2)==jj);
            stim_marker2(dummyindx,2) = c5;
            stim_marker_time(dummyindx,2) = c5;
            c5=c5+1; 
        end
    end
    %}

end






%% define baseline trials


if  cfg4.flag == 1 
    
    
    
    %how much baseline do we have?
    hmb = stim_marker2(1,1)/raw.fsample;
    kf=30;
    if hmb<32
        warning('we dont have much baseline')
        kf = 15;
        
        %patch some pts
        if pt == 3031
            kf = 10;
        end
        
        if pt == 3028 && cfg4.fileidd == 4
            kf=2; %use this for DBS3028 set 7 to 12
        end
        
    end

    %cut 30 trials
    clear stim_marker3 stim_marker_time3
    c = kf;
    for k = 1:kf
        stim_marker3(c,1) = stim_marker2(1,1) - raw.fsample - k*raw.fsample;
        stim_marker3(c,2) = 25;

        stim_marker_time3(c,1) = raw.time{1}(stim_marker3(c,1));
        stim_marker_time3(c,2) = cfg4.sett;

        c=c-1;

    end

    clear stim_marker4 stim_marker_time4
    stim_marker2 = [stim_marker3; stim_marker2];

    stim_marker_time = [stim_marker_time3; stim_marker_time];

end





%% define baseline trials (for 30kHz)

if cfg4.flag == 30 && cfg4.sett == 25

    
    for d = 1:30
        
        r = randn(1);
        
        stim_marker2(d,1) = 1 + 30000*(d-1) + r*30000;
        stim_marker2(d,2) = 0;
        
        stim_marker_time(d,1) = raw.time{1}(1) + (d-1) + r;
        stim_marker_time(d,2) = 0;
    end
    
    

end


%% plot markers


subplot(2,1,1)
plot(stim_marker_time(:,1),ones(length(stim_marker_time),1),'r*')
tt = sprintf('DBS%4.0f-sett:%2.0f median voltage - channel %3.0f',pt,cfg4.sett,indx_median);
title(tt)

subplot(2,1,2)
tt = sprintf('chosen channel to gather stim markers - channel %3.0f',indx);
title(tt)
plot(stim_marker_time(:,1),ones(length(stim_marker_time),1),'r*')


xlabel('time (secs)')
ylabel('Voltage (uV)')




%% plot selected channels for visual downstream checks

if cfg4.sett ~= 999 && cfg4.sett ~=25

    counter57=1;
    for channn = [1 26 52]

        figure(44)
        tt = sprintf('DBS%4.0f-sett:%2.0f all trials channel %3.0f',pt,cfg4.sett,channn);
        subplot(3,1,counter57)
        plot(raw.trial{1}(channn,:),'k'); hold on
        hold on
        plot(stim_marker2(:,1),ones(length(stim_marker2),1),'r*')
        hold on
        title(tt)
        set(gcf,'Position',[57         780        2142         414])

        figure(45)
        subplot(3,1,counter57)
        tt = sprintf('DBS%4.0f-sett:%2.0f 1st trial channel %3.0f',pt,cfg4.sett,channn);
        ti1 = [stim_marker2(1,1)-1500 stim_marker2(2,1)+1500];
        plot(raw.trial{1}(channn,ti1(1):ti1(2)),'k'); hold on
        hold on
        plot([1500 30000+1500]  ,[1 1],'r*'); hold on
        ylim([-400 400])
        counter57=counter57+1;
        title(tt)

        ylim([-400 400])
        title(tt)
        set(gcf,'Position',[743    41   909   728])

    end
    
    % save plots
    cd('/Users/ajorge/ajorge/data/ecog_stnpathways/stim_markers/')
    figure(22)
    filenamee = sprintf('%4.0f_set%2.0f_stim_marker_chosenchannel.png',pt,cfg4.sett);
    saveas(gcf,filenamee)
    figure(44)
    filenamee = sprintf('%4.0f_set%2.0f_stim_marker_all.png',pt,cfg4.sett);
    saveas(gcf,filenamee)
    figure(45)
    filenamee = sprintf('%4.0f_set%2.0f_stim_marker_firsttrial.png',pt,cfg4.sett);
    saveas(gcf,filenamee)
    
end





%% go thru each set, count number of trials per set. it should be around 30,
%otherwise throw error code
if cfg4.sett ~= 999
    for k7 = 1:max(stim_marker2(:,2))
        tot_trials = length(find(stim_marker2(:,2)==k7));
        
        
        %if patient/trial is one with 10Hz stim data:
        if pt == 3031
            if tot_trials<290 || tot_trials>320
                warning('trial numbers may be incorrect - check plots\n')
                %pause
            end
        else
            if tot_trials<25 || tot_trials>45
                warning('trial numbers may be incorrect - check plots\n')
                %pause
            end
        end
        
        
        
    end
end
    

















