% SCRIPT
%
%
% make individual channel plots
%
% 190926 this code is still on the works.




%% set up
clear all; clc

pt = 3028;
sett = 1;

%load 
cd('/Users/ajorge/ajorge/src/stn_pathways/')
load('ecog_brodmann_lookupt','ecog_brodmann_lookupt')

%% calculate baseline t-statistic once
[tstatb_975, tstatb_025] = pro2avg_dummy_baseline(pt);    


%% load pro (all trials for that sett and patient combo)

clear pro
filenamee = sprintf('pro_DBS%4.0f_sett%02.0f',pt,sett);
cd /Users/ajorge/ajorge/data/ecog_stnpathways/pro/
load(filenamee,'pro')

%find all the trials that involved this sett
indx = find(pro.trialinfo==sett);

%go thru all channels
choix = length(pro.label);
if choix < 100
    choixx = 1:63;
elseif choix > 100
    choixx = [1:63 65:127];
end
    
    
%%

choi    = 10;
trace   = 1;

clc

indx2 = find(pro.time{indx(1)}>0.002  & pro.time{indx(1)}<=0.015);     %this time in seconds
clear y_smoothed2

%go thru all trials for that channel
for i = 1:length(indx)


    clear x y p mu f_y  y_detrended y_smoothdetrend_mean y_detrended y_smoothdetrend_mean y_detrended_mean
    %detrend the raw line
    x = pro.time{indx(i)}(indx2);
    y = pro.trial{indx(i)}(choi,indx2);
    [p,~,mu] = polyfit(x,y,8);
    f_y = polyval(p,x,[],mu);
    y_detrended(i,:) = y - f_y;

    %smooth
    y_smoothed2(i,:) = smooth(y_detrended(i,:),30);

    %{
    figure(54)
    subplot(3,1,1)
    plot(x,y); hold on; plot(x,f_y);
    subplot(3,1,2)
    %plot(x,y_detrended); 
    hold on; plot(x,y_smoothed2);
    %}


    %if delta voltage 2ms - 4ms is not >1uV, then take out trial
    %{
    indx3 = find(x>0.002  & x<=0.004);
    miny = min(y_smoothed2(i,indx3));
    maxy = max(y_smoothed2(i,indx3));
    if abs(miny) + abs(maxy) > 1
    else
        y_smoothed2(i,:) = smooth(y_detrended(i,:),30)*nan;
    end
    %}
    
    
    figure(324)
    plot(x, y_smoothed2(i,:))
    hold on

end


%%

%------------------------------------------------------------------
% perform t-statistic at every bin
clear tstat_sett
for bin = 1:size(y_smoothed2,2)
    [h,p,ci,stats] = ttest(y_smoothed2(:,bin));
    tstat_sett(bin) = stats.tstat;
end
%figure(67)
%plot(tstat_sett); hold on; plot(tstatb_975(choi)*ones(length(tstat_sett),1)); plot(tstatb_025(choi)*ones(length(tstat_sett),1));



%------------------------------------------------------------------
%get average trace of the detrended and smoothed line
%y_smoothed_mean2 = smooth(mean(y_smoothed2,1),2);
y_smoothed_mean2 = mean(y_smoothed2,1);


 subplot(3,1,3)
 plot(x,y_smoothed_mean2)
 hold on
%}



%------------------------------------------------------------------
%plot mean trace
%{
if choi <64
    figure(302)
    subplot(9,7,counter9)
elseif choi>63
    figure(303)
    subplot(9,7,counter9-63)
end
counter9=counter9+1;
plot(pro.time{indx(i)}(indx2), y_smoothed_mean2,'k'); hold on
titlee = sprintf('ch %3.0f',choi);
title(titlee)
ylim([-5 5])
xlim([0 0.010])
if choi == 1 || choi == 65
    xlabel('time (s)')
    ylabel('Voltage (uV)')
    titlee = sprintf('sett %3.0f',sett);
    title(titlee);
end
%}



%------------------------------------------------------------------
%save into avg_dummy

avg_dummy.data{trace} = y_smoothed_mean2;
avg_dummy.time{trace} = pro.time{indx(i)}(indx2);

avg_dummy.tstat_sett{trace} = tstat_sett;
avg_dummy.tstatb_975{trace} = tstatb_975(choi);
avg_dummy.tstatb_025{trace} = tstatb_025(choi);


%other bookkeeping
avg_dummy.channel{trace}              = pro.label{choi};
avg_dummy.electrode_n(trace)          = choi;
avg_dummy.sett(trace)                 = sett;
avg_dummy.pt(trace)                   = pt;

%avg_dummy.brodmann(trace)             = ecog_cortexlocs_manual(pt,choi); %deprecrated
avg_dummy.brodmann(trace)             = ecog_brodmann_lookupt(pt,choi); %to get this table, run: SCRIPT_makeBrodmannTable.m



[tract_was_chosen, stim_amp, stim_freq] = chosen_tract_during_sx(pt,sett);
avg_dummy.sett_was_chosen(trace)      = tract_was_chosen;
avg_dummy.stim_amp(trace)             = stim_amp;
avg_dummy.stim_freq(trace)            = stim_freq;








    
    
  
    