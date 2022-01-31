function [tstatb_975, tstatb_025] = pro2avg_baseline(pt)



%% load the baseline sett, calculate the t-statistic

%load file
filenamee = sprintf('pro_DBS%4.0f_sett%02.0f',pt,25);
cd /Users/ajorge/ajorge/data/ecog_stnpathways/pro/
load(filenamee,'pro')

%find all the trials that involved this sett (for bookkepping)
indx = find(pro.trialinfo==25);

%figure out how many ecog arrays (1 or 2 arrays)
counter9=1;
choix = length(pro.label);
if choix < 100
    choixx = 1:63;
elseif choix > 100
    choixx = [1:63 65:127];
end

%go thru all channels
for choi = choixx

    clc
    fprintf('Calculating baseline stats, pt=%4.0f, sett=%2.0f, on choi = %3.0f\n',pt,25,choi)

    indx2 = find(pro.time{indx(1)}>0.001  & pro.time{indx(1)}<=0.020);     %this time in seconds

    %go thru all trials for that channel
    clear y_smoothed2 y_detrended y_smoothdetrend_mean y_detrended y_smoothdetrend_mean y_detrended_mean
    for i = 1:length(indx)

        %detrend the raw line
        x = pro.time{indx(i)}(indx2);
        y = pro.trial{indx(i)}(choi,indx2);
        [p,~,mu] = polyfit(x,y,8);
        f_y = polyval(p,x,[],mu);
        y_detrended(i,:) = y - f_y;

        %smooth
        y_smoothed2(i,:,choi) = smooth(y_detrended(i,:),30);
        
    end
    
    %calculate t-statistic
    
    clear tstat
    for k = 1:1000
        %pick random trial
        xmin=2;
        xmax=size(y_smoothed2,1)-1;
        n=1;
        troi=round(xmin+rand(1,n)*(xmax-xmin));

        %pick 30 random voltage values within the trial
        xmin=2;
        xmax=length(y_smoothed2(troi,:,choi))-1;
        n=300;
        valoi=round(xmin+rand(1,n)*(xmax-xmin));
        vals = y_smoothed2(i,valoi,choi);
        
        [h,p,ci,stats] = ttest(vals);
        tstat(k) = stats.tstat;
    end
    
    tstatb_975(choi) = quantile(tstat,0.975);    
    tstatb_025(choi) = quantile(tstat,0.025); 

end

















