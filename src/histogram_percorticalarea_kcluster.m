function histogram_percorticalarea_kcluster(avg,avg_aux)







clc
close all
clear dummymat_t dummymat_v

stimoi = 3; %either 1mA, 2mA or 3mA
inx = find(avg.stim_amp==stimoi)% & avg.sett_was_chosen==1);

c=1
for klusteri = 1%[1 2 3]

    dummymat_t = nan*ones(10000,50);
    dummymat_v = nan*ones(10000,50);

    %go thru each trace
    for trace = 1:length(inx)
        trace

        %on a specific trace, go thru each significant peak (at 2ms, 5ms, 10ms, etc)  
        for m = 1:length(avg.EPt{inx(trace)})

            %if that significant peak is within the klusteclr timerange...
            if avg.EPt{inx(trace)}(m)>=kluster(klusteri,1,stimoi)/1000 && avg.EPt{inx(trace)}(m)<=kluster(klusteri,2,stimoi)/1000
                
                %if no outliers within klusteroi (less tha 10 uV)
                i8 = find(avg.time{inx(trace)} > kluster(klusteri,1,stimoi)/1000 & avg.time{inx(trace)} < kluster(klusteri,2,stimoi)/1000);
                if max(abs(avg.data{inx(trace)}(i8))) < 10 & max(abs(avg.data{inx(trace)}(i8))) > 0.1

                    %figure out its brodmann area
                    if avg.brodmann(inx(trace))>0 && avg.brodmann(inx(trace))<100
                        dummymat_t(c,avg.brodmann(inx(trace))) = avg.EPt{inx(trace)}(m);
                        dummymat_v(c,avg.brodmann(inx(trace))) = avg.EPv{inx(trace)}(m);

                        %plot the subset of traces we are analyzig here
                        %plottrace2subplot(avg.brodmann(inx(trace)), avg.time{inx(trace)}, avg.data{inx(trace)}, avg.EPt{inx(trace)}(m),  avg.EPv{inx(trace)}(m));
                        
                        c=c+1;
                    end
                    
                else
                end
                
            end
        end
    end

    figure(43)
    subplot(2,3,klusteri+3)
    boxplot(dummymat_t(:,[6 4 3 22 43 44 45 40]),'Labels',{'Premotor','M1','S1','STG','Subcentral','Pars T.','Pars O.','Supramarginal'})
    
    mu = nanmean(dummymat_t(:,[6 4 3 22 43 44 45 40]))*1000
    si = nanstd(dummymat_t(:,[6 4 3 22 43 44 45 40]))*1000
    
    
    
    
    subplot(2,3,klusteri)
    dummy_v{stimoi} = dummymat_v;
    boxplot(dummymat_v(:,[6 4 3 22 43 44 45 40]),'Labels',{'Premotor','M1','S1','STG','subcentral','Pars T.','Pars O.','Supramarginal'})
    ylim([0 10])
    set(gcf,'Position',[  789   180   993   789])
    
     mu = nanmean(dummymat_v(:,[6 4 3 22 43 44 45 40]))*1000
    si = nanstd(dummymat_v(:,[6 4 3 22 43 44 45 40]))*1000

    
    % 1-way anova
    figure(40)
    [p,tbl,stats] = anova1(dummymat_v(:,[6 4 3 22 43 44 45 40]),{'Premotor','M1','S1','STG','Subcentral','Pars T.','Pars O.','Supramarginal'});
    %[p, tbl, stats] = kruskalwallis(dummymat_v(:,[6 4 3 22 43 44 45 40]),{'Premotor','M1','S1','STG','Subcentral','Pars T.','Pars O.','Supramarginal'});
    set(gcf,'Position',[180   577   568   410])
    %multiple comparisons
    figure(20)
    [c_v{klusteri},m_v,~,~] = multcompare(stats,'CType','hsd');
    xlabel('Amplitude (\muV)')
    set(gcf,'Position',[  1837         710        284         177])
    namee = sprintf('kluster%1.0f_voltage',klusteri);
    saveas(gcf,namee,'epsc')


    % 1-way anova
    figure(34)
    [p,tbl,stats] = anova1(dummymat_t(:,[6 4 3 22 43 44 45 40]),{'Premotor','M1','S1','STG','Subcentral','Pars T.','Pars O.','Supramarginal'});
    %[p,tbl,stats] = kruskalwallis(dummymat_t(:,[6 4 3 22 43 44 45 40]),{'Premotor','M1','S1','STG','Subcentral','Pars T.','Pars O.','Supramarginal'});
    set(gcf,'Position',[188    82   560   420])
    %multiple comparisons
    figure(21)
    [c_t{klusteri},m_t,~,~] = multcompare(stats,'CType','hsd');
    xlabel('Latency (ms)')
    set(gcf,'Position',[  1820         226       284         177])
    namee = sprintf('kluster%1.0f_time',klusteri);
    saveas(gcf,namee,'epsc')
%}
    
end

