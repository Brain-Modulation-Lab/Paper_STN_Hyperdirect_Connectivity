function histogram_percorticalarea(avg,avg_aux)



%% simple per cortical area comparison

clc
%close all


for condition = [1 3 5]

    matt = nan*ones(2000,50);
    matv = nan*ones(2000,50);

    m=2;

    aoi = [6 4 3 22 43 44 45 40];

    for coi = aoi
        c=1;

        %tracei = 1:length(avg.stim_amp);
        %tracei = find(avg.stim_amp==3);
        if condition == 1
            tracei = find(avg.brodmann==coi & avg.stim_amp==3);
        elseif condition == 3
            tracei = find(avg.brodmann==coi & avg.sett_was_chosen==1);
        elseif condition == 5
            tracei = find(avg.brodmann==coi & avg.sett_was_chosen==0	);
        end

        %go thru each trace...
        for i = 1:length(tracei)
            matt(c,coi) = avg.EPt{tracei(i)}(m);
            matv(c,coi) = avg.EPv{tracei(i)}(m);
            c=c+1;        
        end

    end

    indx = find(matv<-0.01 | matv>10);%&      matt>0.002 & matt<0.014);
    matt(indx) = nan;
    matv(indx) = nan;


    %quick figure
    figure(30)
    subplot(3,2,condition)
    %boxplot(matt(:,[6 4 3 22 43 44 45 40])*1000,'Labels',{'Premotor','M1','S1','STG','Subcentral','Pars T.','Pars O.','Supramarginal'},'PlotStyle','compact','Colors','rgbm')
    boxplot(matt(:,[6 4 3 22 43 44 45 40])*1000,'Labels',{'Premotor','M1','S1','STG','SC','Pars T.','Pars O.','SM'},'PlotStyle','compact','Colors',avg_aux.cortical_colors([6 4 3 22 43 44 45 40],:),'OutlierSize',1)
    %set(gca,'XTickLabel',{' '})
    mut = nanmedian(nanmean(matt))*1000;
    hold on;
    plot([0 9],[mut mut],'k--')
    ylabel('Latency (ms)')
    set(gca,'Xlabel',[]) 



    subplot(3,2,condition+1)
    boxplot(matv(:,[6 4 3 22 43 44 45 40]),'Labels',{'Premotor','M1','S1','STG','SC','Pars T.','Pars O.','SM'},'PlotStyle','compact','Colors',avg_aux.cortical_colors([6 4 3 22 43 44 45 40],:),'OutlierSize',1)
    %set(gca,'XTickLabel',{' '})
    %set(gca, 'YScale', 'log')
    muv = nanmedian(nanmean(matv));
    hold on;
    plot([0 9],[muv muv],'k--')
    ylabel('Amplitude (\muV)')
    set(gca,'Xlabel',[]) 


    
end



print -depsc -tiff -r300 -painters filename.eps







%% statiscics (voltage amplitude)

mu = nanmean(matv(:,[6 4 3 22 43 44 45 40]),1);
st = nanstd(matv(:,[6 4 3 22 43 44 45 40]),1);
for z = 1:8
    zz = [6 4 3 22 43 44 45 40];
    fprintf('brodmann%2.0f   %2.1f%2.1f\n',zz(z), mu(z),st(z))
end

mu = quantile(matv(:,[6 4 3 22 43 44 45 40]),0.50);
st = quantile(matv(:,[6 4 3 22 43 44 45 40]),0.25);
for z = 1:8
    zz = [6 4 3 22 43 44 45 40];
    fprintf('brodmann%2.0f   %2.1f %2.1f\n',zz(z), mu(z),st(z))
end



figure(40)
[p,tbl,stats] = anova1(matv(:,[6 4 3 22 43 44 45 40]),{'Premotor','M1','S1','STG','Subcentral','Pars T.','Pars O.','Supramarginal'});
%[p, tbl, stats] = kruskalwallis(matv(:,[6 4 3 22 43 44 45 40]),{'Premotor','M1','S1','STG','Subcentral','Pars T.','Pars O.','Supramarginal'});
set(gcf,'Position',[180   577   568   410])
%multiple comparisons
figure(41)

[c_v,m_v,~,~] = multcompare(stats,'CType','hsd','alpha',0.05);
xlabel('Amplitude (\muV)')
set(gcf,'Position',[  1837         710        284         177])
namee = sprintf('kluster%1.0f_voltage',2);
saveas(gcf,namee,'epsc')







%%  statiscics (latency )

% I want to drive home the fact that the "variance" of values in latency is
% not that big. So select values only in the 2.5 to 3.5 ms range and do
% stats

matt2 = nan(2000,50);
for cols = 1:50
    clear inx2
    inx2 = find(matt(:,cols)>2.5*10^-3 & matt(:,cols)<3.5*10^-3);
    matt2(1:length(inx2),cols) = matt(inx2,cols);
end

mu = nanmean(matt2(:,[6 4 3 22 43 44 45 40]),1);
st = nanstd(matt2(:,[6 4 3 22 43 44 45 40]),1);
for z = 1:8
    zz = [6 4 3 22 43 44 45 40];
    fprintf('brodmann%2.0f   %2.1f%2.1f\n',zz(z), mu(z),st(z))
end

mu = quantile(matt2(:,[6 4 3 22 43 44 45 40]),0.50);
st = quantile(matt2(:,[6 4 3 22 43 44 45 40]),0.25);
for z = 1:8
    zz = [6 4 3 22 43 44 45 40];
    fprintf('brodmann%2.0f   %2.1f %2.1f\n',zz(z), mu(z),st(z))
end



figure(40)
%[p,tbl,stats] = anova1(matt2(:,[6 4 3 22 43 44 45 40]),{'Premotor','M1','S1','STG','Subcentral','Pars T.','Pars O.','Supramarginal'});
[p, tbl, stats] = kruskalwallis(matt2(:,[6 4 3 22 43 44 45 40]),{'Premotor','M1','S1','STG','Subcentral','Pars T.','Pars O.','Supramarginal'});
set(gcf,'Position',[180   577   568   410])
%multiple comparisons
figure(41)

[c_v,m_v,~,~] = multcompare(stats,'CType','hsd','alpha',0.05);
xlabel('Amplitude (\muV)')
set(gcf,'Position',[  1837         710        284         177])
namee = sprintf('kluster%1.0f_voltage',2);
saveas(gcf,namee,'epsc')

















%% complex ANOVA N

clc
c2=1;
clear yy subject chosen cortex stim
m=2;

allsubj = unique(avg.pt);

for subj = allsubj%[3005 3008 3010 3011 3017 3018 3023 3030]
    
    %tracei = 1:length(avg.stim_amp);
    %tracei = find(avg.stim_amp==3);
    tracei = find(avg.pt==subj & avg.stim_amp==3);

    %go thru each trace...
    for i = 1:length(tracei)
        
        if avg.EPv{tracei(i)}(m) > 1 & avg.EPv{tracei(i)}(m) < 5 & isnan(avg.EPv{tracei(i)}(m)) == 0 
            yy(c2)        = avg.EPv{tracei(i)}(m);
            subject(c2)   = avg.pt(tracei(i));
            chosen(c2)    = avg.sett_was_chosen(tracei(i));
            cortex(c2)    = avg.brodmann(tracei(i));
            stim(c2)      = avg.stim_amp(tracei(i));
            c2=c2+1;
        end
       
    end


end




%% anovan
clc
close all

figure(31)
[p,table,stats] = anovan(yy,{subject,chosen,cortex,stim},'model',1,'varnames',{'subject','chosen','cortex','stim'})

%[p,table,stats] = anovan(yy,{subject},'model',1,'varnames',{'subject'})

figure(32)
results = multcompare(stats,'Dimension',[2 3])







