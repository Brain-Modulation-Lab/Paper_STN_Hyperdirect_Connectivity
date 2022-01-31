%SCRIPT number of electrodes, by Brodmann area


%% load
cd('/Users/ajorge/Desktop')
load('avg','avg')

unique(avg.pt)

%%

brod = 0;

%find unique patients
uniquept = unique(avg.pt);

for k = 1:length(uniquept)
    %find the index for i unique pt
    inx1 = find(avg.pt == uniquept(k));

    %how many sets
    uniquesets = unique(avg.sett(inx1));

    %pick the first set to localete unique channels associated with this ptient
    %(only do once)
    inx2 = find(avg.pt == uniquept(k) & avg.sett == uniquesets(1));

    %identify brodmann area
    brod = [brod avg.brodmann(inx2)];
end

inx5 = find(brod>50)
brod(inx5)=[];

figure(41)
histogram(brod)


%%
figure(42)
clear n
clc
m=1;
for k = [45 44 6 4 43 3 22 40]
    n(m) = length(find(brod==k));
    m=m+1;
end
c = categorical({'Pars T.','Pars O.','Premotor','M1','Subcentral','S1','STG','SM'});
c = reordercats(c,{'Pars T.','Pars O.','Premotor','M1','Subcentral','S1','STG','SM'});
bar(c,n)

ylabel('Electrode Count')













