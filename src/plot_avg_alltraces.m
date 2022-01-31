function plot_avg_alltraces(avg,avg_aux)
%
%
% ajorge
% 
% last updated Feb 8 2021
%



%%

clc
close all

clear traceavg
c=1;
inxc = 0;





% pick condition
condition = 1; %1=all traces, 2=chosen tract, 3=non chosen tract ALL, 4=non chosen tract randn





for pt = [3012 3016 3019]%[3005        3008        3010        3011        3012        3016        3017        3018        3019        3023        3024        3026  3027        3028        3029        3030        3032]

%for pt = [       4079        4081        4085]
    
sprintf('using pt=%4.0f\n',pt)

    
    if condition == 1
        inx = find(avg.pt==pt );
        colorr=1;
        filenamee = sprintf('grand_avg_allpts_all.eps');
    elseif condition == 2
        inx = find(avg.pt==pt & avg.sett_was_chosen==1);
        colorr=2;
        filenamee = sprintf('grand_avg_allpts_chosentract.eps');
    elseif condition == 3
        inx = find(avg.pt==pt & avg.sett_was_chosen==0);
        colorr=3;
        filenamee = sprintf('grand_avg_allpts_nonchosentract.eps');
    elseif condition == 4
        
        %randomly select 30% of traces from avg.sett_was_chosen==0
        %{
        inx = find(avg.pt==pt & avg.sett_was_chosen==0);
        howmany = round(length(inx)*0.30);
        a = 1;
        b = length(inx);
        r = round((b-a).*rand(howmany,1) + a);
        clear inx2
        inx2 = inx(r);
        inx = inx2;
        %save
        load('avg_auxiliary.mat')
        avg_aux.grand_avg_settings.inx{pt} = inx;
        save('avg_auxiliary.mat','avg_aux');
        %}
        inx = avg_aux.grand_avg_settings.inx{pt};
        colorr=3;
        filenamee = sprintf('grand_avg_allpts_nonchosentract_randn.eps');
    end
    
    
    

    
    inxc = inxc+length(inx);
    
    clear traceavg
    lenn = length(avg.data{inx(1)});
    
    figure(4)
    for i = 1:length(inx)
        if max(avg.data{inx(i)}(50:200))<5 & mean(avg.data{inx(i)})>0.01
            traceavg(i,1:lenn) = avg.data{inx(i)};
            
            plot(avg.time{1}, smooth(traceavg(i,:),20),'k')
            ylim([-5 5])
            hold on
        else
            traceavg(i,1:lenn) = nan(1,lenn);
        end
    end

    figure(200)
    subplot(4,5,c)
    plot(avg.time{1}, smooth(nanmean(traceavg,1),20),'g')
    titlee = sprintf('DBS%4.0f',pt);
    title(titlee)
    hold on
    
    traceavg2(c,:) = nanmean(traceavg,1);

    ylim([-3 3])
    
    plot([0.003 0.003],[-3 3],'--k')
    plot([0.005 0.005],[-3 3],'--k')
    
    c=c+1;
    
    
end






%%
%{
close all
figure(300)

clc

%stich the baseline, stim artifact and latent period
npts = 17 % size(traceavg2,1);

switch condition
    case 1
        x1 = avg_aux.burntends.xtotal(1,1:299);
        y1 = avg_aux.burntends.ytotal{3}(1:30,1:299)*0.05;
        
        x2 = avg.time{1}(1:253);
        y2 = nan*ones(30,253);
        %y2(1:npts,1:253) = traceavg2([5 8 11],1:253);
        y2(1:npts,1:253) = traceavg2(:,1:253);
        % couldn't figure out what numbers were used for original plot.
        % plotting this below
        
        
        x3 = avg_aux.burntends.xtotal(1,552:end);
        y3 = avg_aux.burntends.ytotal{3}(1:30,552:end)*0.05;
        
        x = [x1 x2 x3];
        %y = smooth(nanmean([y1 y2 y3]),15);
        y = smooth(nanmean([y1 y2 y3]),15);
        
        %e = smooth(nanstd([y1 y2 y3]),15)*2;
        e = smooth(nanstd([y1 y2 y3]),15)*2;
        
        
    case 2
        x1 = avg_aux.burntends.xtotal(1,1:299);
        y1 = avg_aux.burntends.ytotal{19}(1:30,1:299)*0.05;
        x2 = avg.time{1}(1:232);
        y2 = nan*ones(30,232);
        y2(1:17,1:232) = traceavg2(:,1:232);
        x3 = avg_aux.burntends.xtotal(1,532:end);
        y3 = -avg_aux.burntends.ytotal{19}(1:30,532:end)*0.05;
        x = [x1 x2 x3];
        y = smooth(nanmean([y1 y2 y3]),20);
        e = smooth(nanstd([y1 y2 y3]),20);
    case 3
        x1 = avg_aux.burntends.xtotal(1,1:299);
        y1 = avg_aux.burntends.ytotal{15}(1:30,1:299)*0.05;
        x2 = avg.time{1}(1:253);
        y2 = nan*ones(30,253);
        y2(1:17,1:253) = traceavg2(:,1:253);
        x3 = avg_aux.burntends.xtotal(1,552:end);
        y3 = avg_aux.burntends.ytotal{15}(1:30,552:end)*0.05;
        x = [x1 x2 x3];
        y = smooth(nanmean([y1 y2 y3]),15);
        e = smooth(nanstd([y1 y2 y3]),15);
    case 4
        x1 = avg_aux.burntends.xtotal(1,1:299);
        y1 = avg_aux.burntends.ytotal{27}(1:30,1:299)*0.05;
        x2 = avg.time{1}(1:253);
        y2 = nan*ones(30,253);
        y2(1:17,1:253) = traceavg2(:,1:253);
        x3 = avg_aux.burntends.xtotal(1,552:end);
        y3 = avg_aux.burntends.ytotal{27}(1:30,552:end)*0.05;
        x = [x1 x2 x3];
        y = smooth(nanmean([y1 y2 y3]),15);
        e = smooth(nanstd([y1 y2 y3]),15);
end

        


plot(x,y)
hold on
plot(x,e)

figure(400+condition)
confplot(x,  y', e',colorr)
%ylim([-1.1 1.1])
ylim([-1 1])
xlim([-0.001 0.012])
inxcc = sprintf('%4.0f traces',inxc);
text(0.008,-0.75,inxcc)
set(gcf,'Position',[  298   738   369   152])
xlabel('Time (sec)')
ylabel('Evoked Potential (\muV)')
hold on

saveas(gcf,filenamee,'epsc')

%}






%%


clear indxs


clc

if condition == 1
    
    
    
    
    for m = [2 3 4]  %evoked potentiaal of iterestet
        c=1;
        clear matt matv mu matthrough
        
        tracei = 1:length(avg.stim_amp);
        %tracei = find(avg.stim_amp==3);
        %tracei = find(avg.brodmann==22);
        
        %go thru each trace...
        for i = 1:length(tracei)
            matt(c) = avg.EPt{tracei(i)}(m);
            matv(c) = avg.EPv{tracei(i)}(m);
            matthrough(c) = avg.EPtrough{tracei(i)}(m);
            c=c+1;
        end



        %--------------------------------------------------------------------------
        % quick histogram plots
        figure(21)
        subplot(2,2,1)
        histogram(matt,50)
        subplot(2,2,3)
        histogram(matv,50)
        %--------------------------------------------------------------------------

        %unselect voltage values that are outliers
        indx = find(matv>-0.1 & matv<10);%&      matt>0.002 & matt<0.014);
        matt2 = matt(indx);
        matv2 = matv(indx);
        
        %save indx of interest for each peak of interest
        indxs{m} = indx;

        %--------------------------------------------------------------------------
        % quick histogram plots
        figure(21)
        subplot(2,2,2)
        histogram(matt2,50)
        subplot(2,2,4)
        histogram(matv2,50)
        %--------------------------------------------------------------------------



        ep_colors = [228,26,28; 55,126,184; 77,175,74]./256;

        figure(54)
        if m==2
            subplot(3,3,m-1)
            histogram(matv2,50,'FaceColor',ep_colors(1,:))
            ylabel('Count')
            xlabel('P1-T1 Amplitude (\muV)')
            set(gca, 'YScale', 'log')
            subplot(3,3,m-1+3)
            histogram(matt2,50,'FaceColor',ep_colors(1,:))
            ylabel('Count')
            xlabel('P1 Latency (s)')
            %set(gca, 'YScale', 'log')
            subplot(3,3,m-1+6)
            %histogram(matthrough,50,'FaceColor',ep_colors(1,:))
            ylabel('Count')
            xlabel('T1 Latency (s)')
            hold on
        elseif m==3

            subplot(3,3,m-1)
            histogram(matv2,50,'FaceColor',ep_colors(2,:))
            ylabel('Count')
            xlabel('P2-T2 Amplitude (\muV)')
            set(gca, 'YScale', 'log')
            subplot(3,3,m-1+3)
            histogram(matt2,50,'FaceColor',ep_colors(2,:))
            ylabel('Count')
            xlabel('P2 Latency (s)')
            %set(gca, 'YScale', 'log')
            subplot(3,3,m-1+6)
            %histogram(matthrough,50,'FaceColor',ep_colors(2,:))
            ylabel('Count')
            xlabel('T2 Latency (s)')
        elseif m==4

            subplot(3,3,m-1)
            histogram(matv2,50,'FaceColor',ep_colors(3,:))
            ylabel('Count')
            xlabel('P3-T3 Amplitude (\muV)')
            set(gca, 'YScale', 'log')
            subplot(3,3,m-1+3)
            histogram(matt2,50,'FaceColor',ep_colors(3,:))
            ylabel('Count')
            xlabel('P3 Latency (s)')
            %set(gca, 'YScale', 'log')
            subplot(3,3,m-1+6)
            %histogram(matthrough,50,'FaceColor',ep_colors(3,:))
            ylabel('Count')
            xlabel('T3 Latency (s)')

        end
        
        
        %statistics
        
        mu = nanmean(matt2);
        si = nanstd(matt2);
        iq1 = quantile(matt2,0.25);
        med = quantile(matt2,0.50);
        iq3 = quantile(matt2,0.75);
        
        fprintf('mean (std) %2.1f %2.1f seconds  \n',mu*1000,2*si*1000)
        fprintf('median=%2.1f (IQR=%2.1f-%2.1f) seconds \n',med*1000,iq1*1000,iq3*1000)
        
        mu = nanmean(matv2);
        si = nanstd(matv2);
        iq1 = quantile(matv2,0.25);
        med = quantile(matv2,0.50);
        iq3 = quantile(matv2,0.75);
        
        fprintf('mean (std) %2.1f %2.1f uvolts \n',mu,2*si)
        fprintf('median=%2.1f (IQR=%2.1f-%2.1f) seconds \n\n',med,iq1,iq3)

        
        
       

    end
    
    
          

end




saveas(gcf,'grand_histogram.eps','epsc')

print -depsc -tiff -r300 -painters filename.eps







%% preliminary plots on mean and median

figure(41)

%inxcombo0 = union(indxs{2},indxs{3});
%inxcombo = union(inxcombo0,indxs{4});
inxcombo = indxs{2};

clear  dum5 dum6 dum7
for i = 1:length(inxcombo)
    %plot(avg.time{i}, avg.data{inxcombo(i)}); hold on
    dum5(i,:) = avg.data{inxcombo(i)};
    dum6(i) = avg.EPv{i}(2);
    dum7(i) = avg.EPtrough{i}(2);
end


figure(41)
subplot(3,1,1)
histogram(dum6,1000)
set(gca, 'YScale', 'log')
xlim([0 10])
subplot(3,1,2)
histogram(dum7*1000,50)
set(gca, 'YScale', 'log')
%xlim([0 10])


figure(42)
plot(avg.time{1}, nanmean(dum5,1))
hold on
plot(avg.time{1}, nanmean(dum5,1) + 1*nanstd(dum5,[],1))
plot(avg.time{1}, nanmean(dum5,1) - 1*nanstd(dum5,[],1))
ylim([-3 3])

figure(43)
plot(avg.time{1}, nanmedian(dum5,1)); hold on
plot(avg.time{1}, quantile(dum5,0.25))
plot(avg.time{1}, quantile(dum5,0.75))
ylim([-2 2])

%standard error
figure(44)
plot(avg.time{1}, nanmean(dum5,1))
hold on
plot(avg.time{1}, nanmean(dum5,1) + 10*nanstd(dum5,[],1)/sqrt(size(dum5,1)))
plot(avg.time{1}, nanmean(dum5,1) - 10*nanstd(dum5,[],1)/sqrt(size(dum5,1)))
ylim([-3 3])




%% plot formal picture of average evoked potential

close all
figure(300)

clc

%stich the baseline, stim artifact and latent period
npts = 17 % size(traceavg2,1);

x1 = avg_aux.burntends.xtotal(1,1:299);
y1 = avg_aux.burntends.ytotal{3}(1:30,1:299)*0.05;



x2 = avg.time{1}(1:253);
%y22 = nanmean(dum5,1);
%y2 = y22(1:253);
y2 = dum5(:,3:255);

x3 = avg_aux.burntends.xtotal(1,552:end);
y3 = avg_aux.burntends.ytotal{3}(1:30,552:end)*0.05;

x = [x1 x2 x3];

y = smooth([nanmean(y1) nanmean(y2) nanmean(y3)],15);
e = smooth([nanstd(y1)*5 nanstd(y2) nanstd(y3)*3],50);


plot(x,y)
hold on
plot(x,e)

figure(400+condition)
confplot(x,  y', e',colorr)
%ylim([-1.1 1.1])
ylim([-1.5 1.5])
xlim([-0.001 0.012])
inxcc = sprintf('%4.0f traces',9450);
text(0.008,-0.75,inxcc)
set(gcf,'Position',[  941   588   471   152])
xlabel('Time (sec)')
ylabel('Evoked Potential (\muV)')
hold on

saveas(gcf,filenamee,'epsc')







