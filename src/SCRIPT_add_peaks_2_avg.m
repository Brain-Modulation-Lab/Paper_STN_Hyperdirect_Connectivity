% SCRIPT avg to avg w/ peak info AND 
%
%
% OPTION 1:
% only adds peaks to avg.EPv and avg.EPt if they are significant
%
% OPTION 2:
% adds all peaks found to avg.EPv and avg.EPt
%



%%
clear all
close all
clc

cd('/Users/ajorge/ajorge/data/ecog_stnpathways/avg/')

%load('avg_190817')
%load('avg_191125')

load('avg_191209_4000s')







%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTION 1

%% automatically find significant peaks

for trace = 1:length(avg.data)

    fprintf('trace %5.0f\n',trace)
    
    [pks,locs,w,p] = findpeaks(avg.data{trace},avg.time{trace},'MinPeakDistance',0.0025,'Annotate','extents');
    
    [pks2,locs2,w2,p2] = findpeaks(-avg.data{trace},avg.time{trace},'MinPeakDistance',0.0025,'Annotate','extents');

    %compare tstat_sett at peak vs tstat_baseline
    for m = 1:length(locs)

        %look for the bin location related to the time_secs of peak
        clear locsb
        locsb(m) = find(avg.time{trace}==locs(m));

        %compare that tstat vs tstat_baseline, only record significant peak if tstat_trace>tstat_975%
        if avg.tstat_sett{trace}(locsb(m)) > avg.tstatb_975{trace}
            try
                avg.EPv{trace}(m) = pks(m);
                avg.EPt{trace}(m) = locs(m);
            catch
                avg.EPv{trace}(m) = nan;
                avg.EPt{trace}(m) = nan;
            end
        else
            avg.EPv{trace}(m) = nan;
            avg.EPt{trace}(m) = nan;
        end        
    end

end










%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTION 2

%quick vizualize an arbitrary trace
figure(44)
trace = 400;
plot(  avg.time{trace} ,   avg.data{trace} )

%%


for trace = 1:length(avg.data)
    
    avg.EPv{trace} = [];
    avg.EPt{trace} = [];

    fprintf('trace %5.0f\n',trace)
    
    %find peaks
    [pks,locs,w,p] = findpeaks(smooth(avg.data{trace}),avg.time{trace},'MinPeakDistance',0.001,'Annotate','extents');
    
    %find throughs
    [pks2,locs2,w2,p2] = findpeaks(-smooth(avg.data{trace}),avg.time{trace},'MinPeakDistance',0.001,'Annotate','extents');

    
    %eliminate pks and pks2 with small prominences
    %inx  = find(p<0.25);  pks(inx)   = [];   locs(inx)  = [];
    %inx2 = find(p2<0.25); pks2(inx2) = [];   locs2(inx2) = [];
    
    %save all amplitudes
    try
        avg.EPv{trace}(1) = p(1);
        avg.EPt{trace}(1) = locs(1);
        
    catch
        avg.EPv{trace}(1) = nan;
        avg.EPt{trace}(1) = nan;
    end
    
    try
        avg.EPv{trace}(2) = pks(2) - -pks2(1);
        avg.EPt{trace}(2) = locs(2);
        avg.EPtrough{trace}(2) = locs2(1);
    catch
        avg.EPv{trace}(2) = nan;
        avg.EPt{trace}(2) = nan;
        avg.EPtrough{trace}(2) = nan;
    end
    
    try
        avg.EPv{trace}(3) = pks(3) - -pks2(2);
        avg.EPt{trace}(3) = locs(3);
        avg.EPtrough{trace}(3) = locs2(2);
    catch
        avg.EPv{trace}(3) = nan;
        avg.EPt{trace}(3) = nan;
        avg.EPtrough{trace}(3) = nan;
    end
    
    try
        avg.EPv{trace}(4) = pks(4) - -pks2(3);
        avg.EPt{trace}(4) = locs(4);
        avg.EPtrough{trace}(4) = locs2(3);
    catch
        avg.EPv{trace}(4) = nan;
        avg.EPt{trace}(4) = nan;
        avg.EPtrough{trace}(4) = nan;
    end
    

end























%% save avg, now with peak info

avg2 = avg;
clear avg;
avg = avg2;
cd('/Users/ajorge/ajorge/data/ecog_stnpathways/avg/')
%save('avg_191209_4000s','avg')



%% visualize EPv


for t = 1:length(avg.data)
    mat(t,1:length(avg.EPv{t})) = avg.EPv{t};   
end

figure(31)
mapp = [247,247,247;253,219,199;244,165,130;214,96,77;178,24,43;103,0,31]./256;
mapp = [0 0 0; 1 1 1];

colormap(mapp)
caxis([0 0.001])
imagesc(mat)



%% plot histograms

clear mat_v mat_t
tt = size(avg.data,2);
mat_v = nan*ones(tt,10);
mat_t = nan*ones(tt,10);

% 3 = primary sensory cortex (including 1,2 and 3 for simplicity purposes)
% 4 = primary motor cortex
% 6 = premotor
% 21 = 22 = 41 = 42 = STG
% 44 = pars opercularis
% 43 = pars triangularis


%choose wich EP
epx = 1;
setwaschosen = 1;

c=1;
for areaa = [6 4 3 22    45 44 43   40]
    
    try
        indx = find(avg.brodmann==areaa)% & avg.sett_was_chosen==1);
        for v = 1:length(indx)
            mat_v(v,c) = avg.EPv{indx(v)}(epx); 
            mat_t(v,c) = avg.EPt{indx(v)}(epx); 
        end
        c=c+1;
        catch
    end
end
    


%take out outliers
indx = find(mat_v>10 | mat_v<0.1);
mat_v(indx) = nan;
mat_t(indx) = nan;


%take out if outside than a time interval
indx = find(mat_t<0.002 | mat_t>0.010);
mat_v(indx) = nan;
mat_t(indx) = nan;


%{
figure(48)
subplot(2,1,1)
boxplot(mat_v)
title('voltages')
hold on
x1 = nanmean(mat_v);
x2 = nanstd(mat_v);


%{
plot([1 2 3 4 5 6 7 8],x1,'or'); hold on
plot([1 2 3 4 5 6 7 8],x1+x2,'ob'); hold on
plot([1 2 3 4 5 6 7 8],x1-x2,'ob'); hold on
%}

subplot(2,1,2)
boxplot(mat_t)
title('seconds')
hold on
x1 = nanmean(mat_t);
x2 = nanstd(mat_t);
hold on
%{
plot([1 2 3 4 5 6 7 8],x1,'or'); hold on
plot([1 2 3 4 5 6 7 8],x1+x2,'ob'); hold on
plot([1 2 3 4 5 6 7 8],x1-x2,'ob'); hold on
%}


%{
clc
[h,p1,ci,stats] = ttest2(mat_v(:,2),mat_v(:,1));
[h,p2,ci,stats] = ttest2(mat_v(:,2),mat_v(:,3));
[h,p3,ci,stats] = ttest2(mat_v(:,2),mat_v(:,4));

[h,p4,ci,stats] = ttest2(mat_t(:,2),mat_t(:,1));
[h,p5,ci,stats] = ttest2(mat_t(:,2),mat_t(:,3));
[h,p6,ci,stats] = ttest2(mat_t(:,2),mat_t(:,4));

fprintf('amplitudes pvals: %4.3f %4.3f %4.3f    latencies pvals:%4.3f %4.3f %4.3f\n',p1,p2,p3,p4,p5,p6)
%}


%paired t-test
%{
mat_v_21 = mat_v(:,2)-mat_v(:,1);
mat_v_23 = mat_v(:,2)-mat_v(:,3);
mat_v_24 = mat_v(:,2)-mat_v(:,4);

mat_t_21 = mat_t(:,2)-mat_t(:,1);
mat_t_23 = mat_t(:,2)-mat_t(:,3);
mat_t_24 = mat_t(:,2)-mat_t(:,4);

[h,p1,ci,stats] = ttest(mat_v_21);
[h,p2,ci,stats] = ttest(mat_v_23);
[h,p3,ci,stats] = ttest(mat_v_24);

[h,p4,ci,stats] = ttest(mat_t_21);
[h,p5,ci,stats] = ttest(mat_t_23);
[h,p6,ci,stats] = ttest(mat_t_24);

figure(35)
subplot(2,1,1)
boxplot([mat_v_21, mat_v_23, mat_v_24])
subplot(2,1,2)
boxplot([mat_t_21, mat_t_23, mat_t_24])

fprintf('amplitudes pvals: %4.3f %4.3f %4.3f    latencies pvals:%4.3f %4.3f %4.3f\n',p1,p2,p3,p4,p5,p6)

%}
%}

[p,tbl,stats] = anova1(mat_v);
multcompare(stats);


[p,tbl,stats] = anova1(mat_t);
multcompare(stats);














































































