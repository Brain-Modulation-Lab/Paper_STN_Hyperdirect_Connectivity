% SCRIPT explore pro
%
% explore how the structure "pro" looks like
%
% ajorge november 2019
%
%
% 


clear all
close all
clc


%% 
clear all
close all
clc

%pt = 3018;
pt = 4085;

%design filter
d3 = designfilt('bandpassfir','FilterOrder',20, 'CutoffFrequency1',300,'CutoffFrequency2',400, 'SampleRate',30000);

for sett = 1% [1 2 4 5 6 7 8 9]
    
    clear x y f_y y_detrended y_detrended4 y_raw y_smoothed y_smoothed2 y_smoothed4 
    
    %load pro data
    filenamee = sprintf('pro_DBS%4.0f_sett%02.0f',pt,sett);
    cd /Users/ajorge/ajorge/data/ecog_stnpathways/pro/
    load(filenamee,'pro')

    % pick your option
    pplotalltrials = 0 ; %pplotalltrials=0 --> plot all channels, pplotalltrials=1 --> plot all trials for one channel,
    
    %select channel
    if pplotalltrials == 0
        choix = [1:22]; arraynum=1; linee=1;
        %choix = [22:42]; arraynum=1; linee=2;
        %choix = [43:63]; arraynum=1; linee=3;
        %choix = [65:127]; arraynum=2; linee=4;
        %choix = [65+21:127+21]; arraynum=2; linee=5;
        %choix = [65+21*2:127+21*2]; arraynum=2; linee=6;
        plotuptotrials = length(pro.trial);
    elseif pplotalltrials == 1
        choix = 10;
        %plotuptotrials = 30;
        plotuptotrials = length(pro.trial);
    end
    
    %which tract (e.g. trialinfo==0 is the baseline tract, 1 is the first tract, etc.)
    indx = find(pro.trialinfo==sett);



    choicount=1;
    %go through each channel
    for choi = choix

        %go thru every trial
        i2=1; 
        for i = 1:plotuptotrials

            %just plot 30 trials
            if i2>30; i2=30; end

            %decide what period of time you want to average with (in seconds)
            
            %settings 1
            %indx2_f = 0.025;
            %indx2 = find(pro.time{indx(i)}>0.0005  & pro.time{indx(i)}<=indx2_f); 
            
            %settings 2
            indx2_f = 0.100;
            indx2 = find(pro.time{indx(1)}>0.001  & pro.time{indx(1)}<=indx2_f);     %this time in seconds

            %index to plot raw data
            indx3 = find(pro.time{indx(i)}>=-0.010 & pro.time{indx(i)}<=indx2_f); 
            %indx3 = find(pro.time{indx(i)}>=-0.010 & pro.time{indx(i)}<=-0.0005); 
            
            %index to plot baseline up to stim
            indx4 = find(pro.time{indx(i)}>=-0.010 & pro.time{indx(i)}<=-0.001); 
            %index to plot stim artifact
            indx5 = find(pro.time{indx(i)}>=-0.001 & pro.time{indx(i)}<=0.0005); 

            

            if sett == 25
                indx2=indx3;
            end


            if pplotalltrials == 1
                %FIRST PLOT plot raw signal (only singal within indx2)
                figure(1000+choi)
                subplot(5,6,i2)

                %option 1 (raw)
                plot(pro.time{indx(i)}(indx3), pro.trial{indx(i)}(choi,indx3),'Color',[128,177,211]./256); hold on

                % option 2 (with filters)
                y55 = pro.trial{indx(i)}(choi,indx2); hold on
                y56 = filtfilt(d3,y55);
                %plot(pro.time{indx(i)}(indx2), y55,'k-.'); hold on
                plot(pro.time{indx(i)}(indx2), y56,'r'); hold on

                titlee = sprintf('trial %2.0f/30',i);
                title(titlee)
                if i == 30
                    xlabel('Time (s)')
                    ylabel('Voltage (uV)')
                end
                xlim([-0.005 indx2_f])
                y5 =  pro.trial{indx(i)}(choi,indx2);
                vmed = median(y5);
                viq1 = quantile(y5,0.01);
                viq3 = quantile(y5,0.99);
                ylim([viq1 viq3])
                %ylim([-30 30])
                set(gcf,'Position',[ 1         623        1161         578])
            end


            %{
            % try filters
            if pplotalltrials == 1
                %highpass 
                x = pro.time{indx(i)}(indx2);
                y = (pro.trial{indx(i)}(choi,indx2));
                %d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',0.15,'DesignMethod','butter');
                %d1 = designfilt('lowpassfir','FilterOrder',4, 'PassbandFrequency',120, 'SampleRate',30000);
                %d1 = designfilt('lowpassfir','PassbandFrequency',0.25, 'StopbandFrequency',0.35,'PassbandRipple',0.5, 'StopbandAttenuation',65,'DesignMethod','kaiserwin');
                %d1 =  designfilt('bandstopiir','FilterOrder',20, 'HalfPowerFrequency1',320,'HalfPowerFrequency2',360, 'SampleRate',30000);
                %d1 = designfilt('bandpassiir','FilterOrder',10, 'HalfPowerFrequency1',35,'HalfPowerFrequency2',55, 'SampleRate',30000);

                y_filtered = filtfilt(d1,y);
                figure(1000+choi)
                subplot(5,6,i)
                plot(x,y_filtered,'r')
            end
            %}



            %smoothed without detrending
            %{
            figure(2000+choi)
            subplot(5,6,i)
            plot(pro.time{indx(i)}(indx2), smooth(pro.trial{indx(i)}(choi,indx2),30) ,     'LineWidth',2); hold on
            %}

            %detrend the raw line
            x = pro.time{indx(i)}(indx2);

            %option 1
            y = (pro.trial{indx(i)}(choi,indx2));

            %option 2 (worked up to march 2019) 
            [p,~,mu] = polyfit(x,y,8);  % using 8 (dec2018-jan2019)
            f_y = polyval(p,x,[],mu); 
            y_detrended(i,:,choi) = y  - f_y;
            y_raw(i,:,choi) = y;
            
                %option 2 (for baseline to stim)
                x4 = pro.time{indx(i)}(indx4);
                y4 = (pro.trial{indx(i)}(choi,indx4));
                [p4,~,mu4] = polyfit(x4,y4,8);  % using 8 (dec2018-jan2019)
                f_y4 = polyval(p4,x4,[],mu4); 
                y_detrended4(i,:,choi) = y4  - f_y4;
                y_raw4(i,:,choi) = y4;
                
                x5 = pro.time{indx(i)}(indx5);
                y5 = (pro.trial{indx(i)}(choi,indx5));
                [p5,~,mu5] = polyfit(x5,y5,8);  % using 8 (dec2018-jan2019)
                f_y5 = polyval(p5,x5,[],mu5); 
                y_detrended5(i,:,choi) = y5  - f_y5;
                y_raw5(i,:,choi) = y5;
                

            %{
            %option 3 (get more baseline)
            %fit a polynomial with indx2
            x = pro.time{indx(i)}(indx2);
            y = (pro.trial{indx(i)}(choi,indx2)); 
            [p,~,mu] = polyfit(x,y,8);  % using 8 (dec2018-jan2019)
            f_y = polyval(p,x,[],mu); 
            %fill values
            y_detrended(i,indx3,choi) = pro.trial{indx(i)}(choi,indx3);
            %refill values
            y_detrended(i,indx2,choi) = y - f_y;
            indx2=indx3;
            %}
            
            % option 3 (with filters)
            %{
            y55 = pro.trial{indx(i)}(choi,indx2); hold on
            y56 = filtfilt(d3,y55);
            y_detrended(i,:,choi) = y56  - f_y;
            y_raw(i,:,choi) = y56;
            %}


            
            
            if pplotalltrials == 1
                %plot raw, detrend line, and detrended signal
                figure(3000+choi)
                subplot(5,6,i2)
                plot(pro.time{indx(i)}(indx2), y, 'Color',[128,177,211]./256); hold on
                plot(pro.time{indx(i)}(indx2), f_y, '--r', 'LineWidth',2);


                %d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',0.15,'DesignMethod','butter');
                %d1 = designfilt('bandstopiir','FilterOrder',20, 'HalfPowerFrequency1',50,'HalfPowerFrequency2',1500,  'SampleRate',30000);
                %y_detrended(i,:,choi) = filtfilt(d1,y_detrended(i,:,choi));
                plot(pro.time{indx(i)}(indx2), y_detrended(i,:,choi),'Color',[102,194,165]./256);


                titlee = sprintf('trial %2.0f/30',i);
                title(titlee)
                if i == 30
                    xlabel('Time (s)')
                    ylabel('Voltage (uV)')
                    legend('raw','detrend','detrended')
                end
                xlim([-0.005 indx2_f])
                %ylim([-30 30])
                set(gcf,'Position',[1158         618        1141         573])
            end



            %smoothed the detrended   
            %y_smoothed(i,:) = sgolayfilt(y_detrended(i,:),7,201);
            y_smoothed2(i,:,choi) = smooth(y_detrended(i,:,choi),20);
            
            y_smoothed4(i,:,choi) = smooth(y_detrended4(i,:,choi),20);
            y_smoothed5(i,:,choi) = smooth(y_detrended5(i,:,choi),20);
            
            
            if pplotalltrials == 1
                figure(4000+choi)
                subplot(5,6,i2)
                plot(pro.time{indx(i)}(indx2), y_detrended(i,:,choi),'Color',[102,194,165]./256); hold on
                %plot(pro.time{indx(i)}(indx2), y_smoothed(i,:), 'LineWidth',2);
                plot(pro.time{indx(i)}(indx2), y_smoothed2(i,:,choi), '-k','LineWidth',2);
                ylim([-30 30])
                xlim([0 indx2_f])
                titlee = sprintf('trial %3.0f',i);
                title(titlee);
                if i == 30
                    xlabel('Time (s)')
                    ylabel('Voltage (uV)')
                    %legend('detrended','smoothed detrended')
                    titlee = sprintf('channel  %3.0f',choi);
                    title(titlee)
                    set(gcf,'Position',[1           1        1156         548])
                end
                
                %plot all trials onto one figure
                figure(6000+choi)
                plot(pro.time{indx(i)}(indx2), y_smoothed2(i,:,choi), '-k','LineWidth',0.1);
                hold on
                plot([0.003 0.003],   [-30 30],'--k')
                plot([0.0055 0.0055], [-30 30],'--k')
                plot([0.0105 0.0105], [-30 30],'--k')
                set(gcf,'Position',[  1158           1        1147         543])
            end           
  
            i2=i2+1;
        end




        %plot the average of the smoothed lines (average all trials)
        y_smoothed_mean2 = mean(y_smoothed2(:,:,choi),1);
        if pplotalltrials == 1
            figure(6000+choi)
            plot(pro.time{indx(i)}(indx2), y_smoothed_mean2,'k','LineWidth',2); hold on 
        end
        %y_raw_mean = mean(y_raw(:,:,choi),1);
        %y_raw_mean_smoothed = smooth(y_raw_mean,30);

        % detrend, then average, then smooth
        y_detrended_mean5 = mean(y_detrended,1);

        if pplotalltrials == 0

            if choicount < 22

                figure(5000+linee)
                subplot(7,3,choicount)

                %plot y_smoothed

                for i = 1:size(y_smoothed2,1)
                %   plot(pro.time{indx(i)}(indx2)*1000, y_smoothed2(i,:,choi),'Color', [150 150 150]./256 ,'LineWidth',0.1); hold on
                end
                %}

                %plot(pro.time{indx(i)}(indx2)*1000, y_raw_mean,'k','LineWidth',2); hold on
                %plot(pro.time{indx(i)}(indx2)*1000, y_raw_mean_smoothed,'k','LineWidth',2); hold on

                %plot(pro.time{indx(i)}(indx2)*1000, y_smoothed_mean2,'r','LineWidth',2); hold on      %chosen plot circa dec2018       
                plot(pro.time{indx(i)}(indx2)*1000, y_smoothed_mean2(),'LineWidth',2); hold on 

                %plot(pro.time{indx(i)}(indx2)*1000, y_detrended_mean5, 'Color', [200 200 200]./256); hold on
                %legend('Detrended Mean','Smoothed')

                xlabel('Time (ms)')
                ylabel('Voltage (uV)')
                xlim([0 100])%indx2_f*1000])

                ylimi = min(y_smoothed_mean2);
                ylimf = max(y_smoothed_mean2);
                if ylimi == ylimf
                    ylimf = ylimi+1;
                end
                ylim([ylimi ylimf])

                ylim([-10 10])

                set(gcf,'Position',[  901   575   325   148])
                titlee = sprintf('channel  %3.0f',choi);
                title(titlee)

            end
            choicount = choicount+1;
        end





    end


end







%% Figure 2 ("trace samples")

choi = choix;

figure(2)
for i = 1:30
    plot(pro.time{indx(i)}(indx2) ,     y_smoothed2(i,:,choi), 'Color', [0.5 0.5 0.5 ]  )
    hold on
    plot(pro.time{indx(i)}(indx4) ,     y_smoothed4(i,:,choi), 'Color', [0.5 0.5 0.5 ]  )
    hold on
    plot(pro.time{indx(i)}(indx5) ,     y_smoothed5(i,:,choi), 'Color', [0.5 0.5 0.5 ]  )
    
end

figure(3)
clear y2 y4 y5 e2 e4 e5

%stdshade(y_smoothed2(:,:,choi), 0.2)
y2 = nanmean(y_smoothed2(:,:,choi),1);
e = nanstd(y_smoothed2(:,:,choi),[],1);
confplot(pro.time{indx(i)}(indx2),  y2, e)
hold on

e = nanstd(y_smoothed4(:,:,choi),[],1);
confplot(pro.time{indx(i)}(indx4)+0.0006,  nanmean(y_smoothed4(:,:,choi),1), e)
hold on

e = nanstd(y_smoothed5(:,:,choi),[],1);
confplot(pro.time{indx(i)}(indx5),  nanmean(y_smoothed5(:,:,choi),1), e)


figure(4)
clear ytotal xtotal ym e 
ytotal = [  y_smoothed4(:,:,choi)                y_smoothed5(:,31:end,choi)                  y_smoothed2(:,5:end,choi)];
xtotal = [  pro.time{indx(i)}(indx4)+0.001       pro.time{indx(i)}(indx5(31:end))          pro.time{indx(i)}(indx2(5:end))  ];
ym = smooth(nanmean(ytotal,1),2);
e = smooth(nanstd(ytotal),2);
confplot(xtotal, ym, e)


xlim([-0.001 0.012])
ylim([-6 6])
xlabel('Time (sec)')
ylabel('Evoked Potential (/uV)')

set(gcf,'Position',[1782         810         323         141])

filenamee = sprintf('figure2_DBS%4.0f_choi%3.0f',pt,choix);
saveas(gcf,filenamee,'epsc')




















%% EXPLORATORY - detrend, CMR (using all), smooth
%
% nothing came out of this as of Nov 2019


y_median = mean(y_smoothed2(:,:,[8:9]),3);
figure(9001)
for trialx = 1:30
    subplot(6,5,trialx)
    plot(y_median(trialx,:));
end


y_cared = y_smoothed2 - y_median;

y_trialavg = mean(y_cared,1);

for choi = 1:21
    figure(5000)
    subplot(7,3,choi)
    y_trialavg2 = smooth(y_trialavg(1,:,choi),10);
    
    plot(pro.time{indx(i)}(indx2)*1000, y_trialavg2,'g.', 'LineWidth',2); hold on
    
    %plot(pro.time{indx(i)}(indx2)*1000, y_median,'g.', 'LineWidth',2); hold on
    
    titlee = sprintf('channel  %3.0f',choi);
    title(titlee)
    xlim([0 20])
    ylim([-5 5])
    
    if choi == 8 %|| choi == 9
        figure(5001)
        plot(pro.time{indx(i)}(indx2)*1000, y_trialavg2,'g.', 'LineWidth',2); hold on
        
        
        %plot(pro.time{indx(i)}(indx2)*1000, mean(y_median,1),'y', 'LineWidth',2); hold on
        
         xlabel('Time (ms)')
        ylabel('Voltage (uV)')
        xlim([0 18])
    end
    
end
%}




