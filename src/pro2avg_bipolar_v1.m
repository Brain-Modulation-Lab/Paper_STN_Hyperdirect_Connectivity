function avg = pro2avg_bipolar_v1(pt,sett,avg,trace,tstatb_975, tstatb_025)

% transforms pro to avg
% (averages 30 trials into one average trace)
% it also calcualtes the t-statistic for the baseline of each channel AND
% the t-statistic for the 30 trials per channel











%% load the experimental sett

clear pro
filenamee = sprintf('pro_DBS%4.0f_sett%02.0f',pt,sett);
cd /Users/ajorge/ajorge/data/ecog_stnpathways/pro/
load(filenamee,'pro')

close all


%load 
cd('/Users/ajorge/ajorge/src/stn_pathways/')
load('ecog_brodmann_lookupt','ecog_brodmann_lookupt')


%find all the trials that involved this sett
indx = find(pro.trialinfo==sett);

%go thru all channels
choix = length(pro.label);
if choix < 100
    choixx = 1:63;
elseif choix > 100
    choixx = [1:63 65:127];
end







%% baseline correct all trials and all channels
for choi = choixx
 
    close all
    clc
    fprintf('Calculating pt=%4.0f, sett=%2.0f, on choi = %3.0f\n',pt,sett,choi)
    
    indx2 = find(pro.time{indx(1)}>0.001  & pro.time{indx(1)}<=0.010);     %this time in seconds

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
        y_smoothed2(i,:,choi) = smooth(y_detrended(i,:),30);
        
    end
    
    
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

end













%% average on a 3x3 grid  (first grid)

%if y_smoothed2 is "noisy", then nan it
for k = 1:size(y_smoothed2,3)

    if sum(max(y_smoothed2(:,50:270,k))>20)>1
        y_smoothed2(:,:,k) = nan;
    end

end




for i = 1:7
    
    choi(1:3) = [1  2  3]  + 3*(i-1);
    choi(4:6) = [22 23 24] + 3*(i-1);
    choi(7:9) = [43 44 45] + 3*(i-1);
    
    y_3x3(:,:,i) = nanmean(y_smoothed2(:,:,choi),3);
    
    %{
    figure(31)
    for j = 1:30
        subplot(1,7,i)
        plot(x, y_3x3(j,:,i))
        hold on
        
        plot([0.003 0.003],[-6 6],'--k')
        hold on
        plot([0.005 0.005],[-6 6],'--k')
    
    end
    %}
    
    
    
end



%% bipolar

for i = 1:6
    
    y_bipol(:,:,i) = y_3x3(:,:,i) - y_3x3(:,:,i+1);
    
    figure(32)
    subplot(1,6,i)
    plot(x, y_bipol(:,:,i),'Color',[0.7 0.7 0.7])
    hold on
    
    plot(x, nanmean(y_bipol(:,:,i),1),'k','LineWidth',3)
    ylim([-3 3])
    plot([0.003 0.003],[-3 3],'--k')
    plot([0.005 0.005],[-3 3],'--k')
   
    
    avg.data{trace} = nanmean(y_bipol(:,:,i),1);
    avg.time{trace} = pro.time{indx(i)}(indx2);
    
    avg.y_3x3{trace} = y_3x3;
    avg.y_bipol{trace} = y_bipol;
    
    %other bookkeeping
    avg.ch3x3(trace)                = i;
    avg.sett(trace)                 = sett;
    avg.pt(trace)                   = pt;
    
    clear choi
    choi(4:6) = [22 23 24] + 3*(i-1);
    choi = choi(6);
    avg.channel{trace}              = pro.label{choi};
    avg.electrode_n(trace)          = choi;
    avg.brodmann(trace)             = ecog_brodmann_lookupt(pt,choi);
    
    if isnan(avg.brodmann(trace))==0
            titlee = brodman_n2s(avg.brodmann(trace));
            title(titlee);
     end
    
    
    [tract_was_chosen, stim_amp, stim_freq] = chosen_tract_during_sx(pt,sett);
    avg.sett_was_chosen(trace)      = tract_was_chosen;
    avg.stim_amp(trace)             = stim_amp;
    avg.stim_freq(trace)            = stim_freq;

    trace=trace+1;

end



filenamee = sprintf('chosen%1.0f_DBS%4.0f_3x3_bipolared_array%1.0f_sett%02.0f.eps',avg.sett_was_chosen(trace-1),pt,1,sett);
cd('/Users/ajorge/ajorge/data/ecog_stnpathways/bipolared_plots/')
set(gcf,'Position',[ 879         823        1332         233])

saveas(gcf,filenamee,'epsc')

close all













%% average on a 3x3 grid  (first grid)

if size(y_smoothed2,3)>63

    %if y_smoothed2 is "noisy", then nan it
    for k = 1:size(y_smoothed2,3)

        if sum(max(y_smoothed2(:,50:270,k))>20)>1
            y_smoothed2(:,:,k) = nan;
        end

    end



    for i = 1:7

        choi(1:3) = [1  2  3]  + 3*(i-1) + 64;
        choi(4:6) = [22 23 24] + 3*(i-1) + 64;
        choi(7:9) = [43 44 45] + 3*(i-1) + 64;     

        y_3x3(:,:,i) = nanmean(y_smoothed2(:,:,choi),3);

        %{
        figure(31)
        for j = 1:30
            subplot(1,7,i)
            plot(x, y_3x3(j,:,i))
            hold on

            plot([0.003 0.003],[-6 6],'--k')
            hold on
            plot([0.005 0.005],[-6 6],'--k')

        end
        %}

    end



    %% bipolar

    for i = 1:6

        y_bipol(:,:,i) = y_3x3(:,:,i) - y_3x3(:,:,i+1);

        figure(32)
        subplot(1,6,i)
        plot(x, y_bipol(:,:,i),'Color',[0.7 0.7 0.7])
        hold on

        plot(x, nanmean(y_bipol(:,:,i),1),'k','LineWidth',3)
        ylim([-3 3])
        plot([0.003 0.003],[-3 3],'--k')
        plot([0.005 0.005],[-3 3],'--k')


        avg.data{trace} = nanmean(y_bipol(:,:,i),1);
        avg.time{trace} = pro.time{indx(i)}(indx2);

        avg.y_3x3{trace} = y_3x3;
        avg.y_bipol{trace} = y_bipol;

        %other bookkeeping
        avg.ch3x3(trace)                = i;
        avg.sett(trace)                 = sett;
        avg.pt(trace)                   = pt;

        clear choi
        choi(4:6) = [22 23 24] + 3*(i-1) + 64;
        choi = choi(6);
        avg.channel{trace}              = pro.label{choi};
        avg.electrode_n(trace)          = choi;
        avg.brodmann(trace)             = ecog_brodmann_lookupt(pt,choi);

        if isnan(avg.brodmann(trace))==0
            titlee = brodman_n2s(avg.brodmann(trace));
            title(titlee);
        end


        [tract_was_chosen, stim_amp, stim_freq] = chosen_tract_during_sx(pt,sett);
        avg.sett_was_chosen(trace)      = tract_was_chosen;
        avg.stim_amp(trace)             = stim_amp;
        avg.stim_freq(trace)            = stim_freq;

        trace=trace+1;

    end



    filenamee = sprintf('chosen%1.0f_DBS%4.0f_3x3_bipolared_array%1.0f_sett%02.0f.eps',avg.sett_was_chosen(trace-1),pt,2,sett);
    cd('/Users/ajorge/ajorge/data/ecog_stnpathways/bipolared_plots/')
    set(gcf,'Position',[ 879         823        1332         233])

    saveas(gcf,filenamee,'epsc')


end














