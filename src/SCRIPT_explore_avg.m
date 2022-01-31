%% Explore avgWITHpeaks




%% load avg
clear all
close all
clc

cd('/Users/ajorge/ajorge/data/ecog_stnpathways/avg/')

%load('avg_190625') %0.1 to 20 ms polynomial
%load('avg_190818') %0.1 to 15 ms polynomial
%load('avg_190817') %0.1 to 10 ms polynomial     0.001 to 0.0100    CHOSEN ONE until 191125

load('avg_191125')  %0.1 to 10 ms polynomial, selecting all peaks (not only the significant ones)

%load('avg_191209_4000s'); %4000s

%load('avg_191210_orthodromics')

load('avg_auxiliary')



%% create avg_aux 
avg_aux.cortical_names{6}  = 'Premotor';    %
avg_aux.cortical_names{4}  = 'M1';
avg_aux.cortical_names{3}  = 'S1';
avg_aux.cortical_names{43} = 'Subcentral';
avg_aux.cortical_names{22} = 'STG';
avg_aux.cortical_names{44} = 'Pars O.';
avg_aux.cortical_names{45} = 'Pars T.'; 
avg_aux.cortical_names{40} = 'Supramarginal';
%avg_aux.cortical_colors([6 4 3 22 43 44 45],:) = [118 89 150;60 96 150;70 145 121;83 77 154;154 78 150;150 44 83;211 226 184]./256;
avg_aux.cortical_colors([6 4 3 43 22 44 45 40],:) = [118 89 150;60 96 150;70 145 121;83 77 154;154 78 150;150 44 83;164 187 137; 96 95 56]./256;

%save('avg_auxiliary','avg_aux')








%% how many significant peaks per patient

clc

clear inx setts_available
ptoi = unique(avg.pt);
ptcounter = zeros(length(ptoi),1);
for i = 1:length(ptoi)
   
   %find all traces related to this patient
   inx = find(avg.pt==ptoi(i));
    
   setts_available(i) = length(unique(avg.sett(inx)));
   ecog_available(i)  = length(unique(avg.electrode_n(inx)));
    
    
    %go thru each trace
    for j = 1:length(inx)
        %go thru each significant peak (at 2ms, 5ms, 10ms, etc)
        for m = 1:length(avg.EPv{inx(j)})
            if avg.EPt{inx(j)}(m)>0.0001 && avg.EPt{inx(j)}(m)<0.014     &&      avg.EPv{inx(j)}(m)>0.01 && avg.EPv{inx(j)}(m)<20
                ptcounter(i,1) =  ptcounter(i,1)+1;
            end
        end
    end
    
end


figure(301)
subplot(3,2,1)
bar(ptcounter)
title('all peaks')
subplot(3,2,2)
boxplot(ptcounter)
%ylim([0 1500])
subplot(3,2,3)
bar(ptcounter./setts_available')
title('all peaks/setts')
subplot(3,2,4)
boxplot(ptcounter./setts_available')
%ylim([0 250])
subplot(3,2,5)
bar(ptcounter./setts_available'./ecog_available')
title('all peaks/setts/ecog_n')
subplot(3,2,6)
boxplot(ptcounter./setts_available'./ecog_available')
%ylim([0 2.5])

% four peaks = the EP0 and EP1-3








%% quick plot - traces with significant peaks in reddots, sorted by area

plot_rnd_individual_avg(avg,avg_aux)


%% plot all traces for a cortical area
brodmannoi = 6;
%plotperarea(brodmannoi)
clc
pts = unique(avg.pt);

c=1;
figure(320)
for ptoi = pts
    subplot(4,5,c)
    inx = find(avg.brodmann==brodmannoi & avg.pt==ptoi);
   
    for i = inx
        if max(avg.data{i})<10;
            plot(avg.time{i}, avg.data{i})
            hold on
        end
    end
    
    ylim([-20 20])
    titlee = sprintf('DBS%4.0f',ptoi);
    title(titlee)
    c=c+1;
end






%% Kluster analysis (DEPRECRATED)
% Plot peak histogram (peak time vs count) CLUSTER ANALYSIS
%peakhistograms_kclusters(avg,avg_aux)

% simple histograms per cortical area + anova1 and multicomparisons
%histogram_percorticalarea_kcluster(avg,avg_aux)





%% average of all traces per subject
% figure 3. short latency cortical evoked potentials
plot_avg_alltraces(avg,avg_aux)


%% simple histograms per cortical area + anova1 and multicomparisons
histogram_percorticalarea(avg,avg_aux)









%% Figure - stim amplitude vs cEP1 (1mA 2mA 3mA) comparison

figure(31)
clear yy
yy = nan(2000,3);
c=1;
EPoi = 3;

%plot each area
for brodmannoi = [45 44 6 4 43 3 22 40]
    
    for stimoi = [1 2 3]
        inx = find(avg.stim_amp==stimoi & avg.brodmann==brodmannoi);
        for i=1:length(inx)
            
            if avg.EPv{inx(i)}(2)>0 && avg.EPv{inx(i)}(2)<10
                yy(i,stimoi)            = avg.EPv{inx(i)}(EPoi);
            end
        end
    end
    
    
    subplot(3,3,c); 
    if c == 1
        ylabel('cEP_1  (\muV)')
        xlabel('Stimulation Amplitude (mA)')
    end
    title(avg_aux.cortical_names{brodmannoi})
    hold on
    
    x = [1 2 3];
    y = nanmean(yy,1);
    e = nanstd(yy,1);
    
    errorbar(x,y,e,'k')
    xlim([0.5 3.5])
    ylim([0 4.5])
    
    c=c+1;
    
end


%plot all areas
clear yy g1 g2 g3
yy = nan(9000,3);

d=1;
for stimoi = [1 2 3]
    inx = find(avg.stim_amp==stimoi);
    for i=1:length(inx)
        if avg.EPv{inx(i)}(EPoi)>0 && avg.EPv{inx(i)}(EPoi)<10
            
            yy(i,stimoi)            = avg.EPv{inx(i)}(EPoi);
            
            g1(d) = avg.EPv{inx(i)}(EPoi);
            g2{d} = sprintf('%s',num2str(avg.pt(inx(i))));
            g3(d) = sprintf('%s',num2str(avg.stim_amp(inx(i))));
            d=d+1;
            
        end
    end
end
subplot(3,3,c); 
title('All Gyri')
hold on

x = [1 2 3];
y = nanmean(yy,1);
e = nanstd(yy,1);
errorbar(x,y,e,'k')
xlim([0.5 3.5])
ylim([0.1 4.5])

set(gcf,'Position',[1724         458         392         390])

figure(31)
cd /Users/ajorge/Desktop/hyperdirectfigs
filenamee = sprintf('stimamp_vs_cEP1.eps');
set(gcf,'Position',[  1660         416         392         390])
saveas(gcf,filenamee,'epsc')
n = sum(isnan(yy(:,3))==0)




% statiscitcs on all cEP (per stimulation amplitude)

figure(40)
%[p,tbl,stats] = anova1(yy,{'1mA','2mA','3mA'});
[p, tbl, stats] = kruskalwallis(yy,{'1mA','2mA','3mA'});
set(gcf,'Position',[180   577   568   410])
%multiple comparisons


figure(41)
[c_v,m_v,~,~] = multcompare(stats,'CType','hsd','alpha',0.05);
xlabel('Amplitude (\muV)')
set(gcf,'Position',[  1837         710        284         177])
namee = sprintf('kluster%1.0f_voltage',2);
saveas(gcf,namee,'epsc')

figure(43)
subplot(1,3,EPoi-1)
xlabel('Stimulation Amplitude (mA)')
labell = sprintf('cEP_%1.0f (uV)',EPoi-1);
ylabel(labell)
ylim([0 6])
hold on
boxplot(yy,'Labels',{'1 mA','2 mA','3 mA'},'OutlierSize',0)

cd /Users/ajorge/Desktop/hyperdirectfigs
filenamee = sprintf('stimamp_vs_cEP1_simple.eps');
saveas(gcf,filenamee,'epsc')

% statiscitcs on all gyri, compare pt interaction

clc

[~,~,stats] = anovan(g1,{g2})% g3})%,'model',2,'varnames',{'pt','stim'});


results = multcompare(stats,'Dimension',[1 2])






%% orthodromic plots

clear y
y = nan(10000,2970);
clc

figure(100)

c2=1;
%for pt = [        3005        3008        3010        3011        3012        3016        3017        3018        3019        3023        3024        3026         3027        3028        3029        3030        3032]
    pt=3005
    for choss = [0 1]
        
        clear y
        y = nan(10000,2970);

        inx = find(avg.sett_was_chosen==choss & avg.stim_amp==3 & avg.brodmann==4);

        c=1;
        for i = inx

            if max( avg.data{i}(100:2500))<200
                %plot(avg.time{i}, avg.data{i},'Color',[0.8 0.8 0.8])
                y(c,1:length(avg.data{i})) = smooth(avg.data{i});
                %hold on
                c=c+1;
            end

        end

        %subplot(5,5,c2)
        
        if choss == 1; plot(avg.time{i},smooth(nanmean(y,1)),'r');   hold on; ylim([-10 10]); end
        if choss == 0; plot(avg.time{i},(nanmedian(y,1)),'b'); hold on; ylim([-20 20]); end
 
    end

    
    c2=c2+1;
%end

































































