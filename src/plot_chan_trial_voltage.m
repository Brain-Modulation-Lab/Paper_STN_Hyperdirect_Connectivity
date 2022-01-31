function plot_chan_trial_voltage(pro,sett,stim_marker2)




% NEED TO COMPLETETLY EDIT THIS FUNCTION


%% prestim
for ch = 1:21
    array(ch,:) = raw.trial{1}(ch,1:stim_marker2(1)-1);   
end

figure(301)
imagesc(zscore(array))

mapp = [5,48,97;33,102,172;67,147,195;146,197,222;209,229,240;247,247,247;253,219,199;244,165,130;214,96,77;178,24,43;103,0,31;]./256;
colormap(mapp)
caxis([-3 3])
%clabel('voltag (uV)')
xlabel('time (ms)')
ylabel('channel')
h = colorbar;
ylabel(h, 'voltage (z-score)')



%% first row of electrodes
clear array2 array3
for set = sett%1:stim_marker2(end,2)
    indx = find(stim_marker2(:,2)==set);
    begg = stim_marker2(indx(1),1);
    endd = stim_marker2(indx(end),1);
    for ch = 1:21
        array2(ch,:) = raw.trial{1}(ch, begg:endd); 
    end
end
figure(302)
subplot(3,1,1)
imagesc(zscore(array2,[],2))

mapp = [5,48,97;33,102,172;67,147,195;146,197,222;209,229,240;247,247,247;253,219,199;244,165,130;214,96,77;178,24,43;103,0,31;]./256;
colormap(mapp)
caxis([-3 3])
%clabel('voltage (uV)')
xlabel('time (ms)')
ylabel('channel')
h = colorbar;
ylabel(h, 'voltage (z-score)')



%% middle row
for set = sett%1:stim_marker2(end,2)
    indx = find(stim_marker2(:,2)==set);
    begg = stim_marker2(indx(1),1);
    endd = stim_marker2(indx(end),1);
    for ch = 22:42
        array3(ch,:) = raw.trial{1}(ch, begg:endd); 
    end
end
figure(302)
subplot(3,1,2)
imagesc(zscore(array3(22:42,:),[],2))

mapp = [5,48,97;33,102,172;67,147,195;146,197,222;209,229,240;247,247,247;253,219,199;244,165,130;214,96,77;178,24,43;103,0,31;]./256;
colormap(mapp)
caxis([-3 3])
%clabel('voltage (uV)')
xlabel('time (ms)')
ylabel('channel')
h = colorbar;
ylabel(h, 'voltage (z-score)')

yticklabels({'26','31','36','41'})



%% last row of electrodes
for set = sett%1:stim_marker2(end,2)
    indx = find(stim_marker2(:,2)==set);
    begg = stim_marker2(indx(1),1);
    endd = stim_marker2(indx(end),1);
    for ch = 43:63
        array4(ch,:) = raw.trial{1}(ch, begg:endd); 
    end
end
figure(302)
subplot(3,1,3)
imagesc(zscore(array4(43:63,:),[],2))

mapp = [5,48,97;33,102,172;67,147,195;146,197,222;209,229,240;247,247,247;253,219,199;244,165,130;214,96,77;178,24,43;103,0,31;]./256;
colormap(mapp)
caxis([-3 3])
%clabel('voltage (uV)')
xlabel('time (ms)')
ylabel('channel')
h = colorbar;
ylabel(h, 'voltage (z-score)')
yticklabels({'47','52','57','62'})


%%










































