function peakhistograms_kclusters(avg,avg_aux)




%%
clc
c=1;
clear matt matv mu
 
tracei = find(avg.stim_amp==3);

%
%go thru each trace...
for i = 1:length(tracei)
        
    %collect all significant peak voltage and time info
    for m = 1:length(avg.EPt{tracei(i)})
        matt(c) = avg.EPt{tracei(i)}(m);
        matv(c) = avg.EPv{tracei(i)}(m);
        c=c+1;
    end
    
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
indx = find(matv>0.001 & matv<20    &      matt>0.002 & matt<0.014);
matt2 = matt(indx);
matv2 = matv(indx);

%--------------------------------------------------------------------------
% quick histogram plots
figure(21)
subplot(2,2,2)
histogram(matt2,50)
subplot(2,2,4)
histogram(matv2,50)
%--------------------------------------------------------------------------

figure(50)
subplot(3,3,[4,5,6])
histogram(matv2,40,'FaceColor','k')
ylabel('Count')
xlabel('cEP amplitude (uV)')

avg_aux.iq1_v_all = quantile(matv2,0.25);
avg_aux.iq2_v_all = quantile(matv2,0.50);
avg_aux.iq3_v_all = quantile(matv2,0.75);
cd('/Users/ajorge/ajorge/data/ecog_stnpathways/avg/')
%save('avg_auxiliary','avg_aux')


subplot(3,3,[1,2,3])
%histogram(matt2,40,'FaceColor','k')
ylabel('Count')
xlabel('cEP latency (ms)')



%k-means the time (ms)
clear idx C sumd D
[idx,C,sumd,D] = kmeans(matt2',3);
% only run kmeans once (it changes cuz random seed)
%cd('/Users/ajorge/ajorge/data/ecog_stnpathways/kmeans')
%load('kmeans_3mA_190610.mat')

clear f
f = [round(40/3) round(40/3) round(40/3)];
%f{1} = [4.20:0.25:8.0]*0.001;
%f{2} = [8.0:0.25:14]*0.001;
%f{3} = [2.0:0.25:4.25]*0.001;

ep_colors = [228,26,28; 55,126,184; 77,175,74]./256;

for i= 1:3
    figure(50)
    subplot(2,1,1)
    inx5 = find(idx==1);
    histogram(matt2(1,inx5),f{i},'FaceColor',ep_colors(i,:))
    hold on
    mu(i) = nanmean(matt2(idx==i));
    st(i) = nanstd(matt2(idx==i));
end
ylabel('Count')
xlabel('cEP latency (ms)')

thefirst = find(mu<0.0042);
themiddl = find(mu>0.0042 & mu<0.008);
thelast  = find(mu>0.008);

%plot voltage distribution for each k-mean
%{
subplot(3,3,7)
histogram(matv2(idx==thefirst),round(40/3),'FaceColor',ep_colors(thefirst,:))
ylabel('Count')
xlabel('EP amplitude (uV)')
subplot(3,3,8)
histogram(matv2(idx==themiddl),round(40/3),'FaceColor',ep_colors(themiddl,:))
ylabel('Count')
xlabel('EP amplitude (uV)')
subplot(3,3,9)
histogram(matv2(idx==thelast),round(40/3),'FaceColor',ep_colors(thelast,:))
ylabel('Count')
xlabel('EP amplitude (uV)')
%}

% plot voltage distribution for each k-mean 
subplot(3,3,[7 8 9]) 
rowvals = [1 2 3]';
maxlen = max(  [ length(matv2(idx==thefirst))     length(matv2(idx==themiddl))    length(matv2(idx==thelast)) ]) ;
mat56 = nan*ones(3,maxlen);
mat56(1,1:length(matv2(idx==thefirst))) = matv2(idx==thefirst);
mat56(2,1:length(matv2(idx==themiddl))) = matv2(idx==themiddl);
mat56(3,1:length(matv2(idx==thelast))) = matv2(idx==thelast);
boxplot(mat56','ori','horizontal','positions',rowvals,'PlotStyle','compact')
ylabel('K-Cluster')
xlabel('cEP amplitude (uV)')



%[p,tbl,stats]  = kruskalwallis(mat56');
%c = multcompare(stats);


fprintf('the k-means mu = %5.3f %5.3f %5.3f   and   std = %5.3f %5.3f %5.3f\n',mu*1000,st*2*1000)

minmaxmat = [
    min(matt2(idx==1)) max(matt2(idx==1)) 
    min(matt2(idx==2)) max(matt2(idx==2)) 
    min(matt2(idx==3)) max(matt2(idx==3))
    ];

fprintf('the k-means range = %5.3f-%5.3f   %5.3f-%5.3f    %5.3f-%5.3f        \n',minmaxmat'*1000)

cd /Users/ajorge/Desktop/
print -depsc -tiff -r300 -painters filename.eps














%% defined time k-clusters

clear kluster

% clustering for 1mA
kluster(1,:,1) = [2.033 4.100];
kluster(2,:,1) = [4.133 8.200 ];
kluster(3,:,1) = [8.233 13.900];

% clustering for 2mA
kluster(1,:,2) = [2.067 4.067 ];
kluster(2,:,2) = [4.567 8.600    ];
kluster(3,:,2) = [9.500 13.367  ];

% clustering for 3mA
kluster(1,:,3) = [2.033 4.200];
kluster(2,:,3) = [4.233 8.033];
kluster(3,:,3) = [8.067 13.967];
