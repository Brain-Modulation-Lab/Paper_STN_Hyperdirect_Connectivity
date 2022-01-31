%SCRIPT_normalize_cEP_by_array
%
%
% ajorge march 28 2021
% the way Chen does it, not used on my paper






%% load avg
clear all
close all
clc

cd('/Users/ajorge/ajorge/data/ecog_stnpathways/avg/')

%load('avg_190817')

load('avg_191125')
load('avg_auxiliary')



%% write strip number (according to electrode number info)

arraycounter=1;
for i = 1:9450

    gf = avg.electrode_n(i);
    
    avg.array(i) = arraycounter;
    
    if gf == 63 || gf == 127
        arraycounter=arraycounter+1;
    end
    
    
end

%check work
clc
i=8000;
f=9450;
[avg.pt(i:f)' avg.electrode_n(i:f)' avg.array(i:f)']


%% plot ERP voltages by array

clear dum
for i = 1:150
    
    inx = find(avg.array==i);
    
    for j=1:63
        dum(j,i) = avg.EPv{inx(j)}(2);
    end
    
end

figure(55)
boxplot(dum)

%% normalize ERP by 95 quantile

for i = 1:150;
    q = quantile(dum(:,i),0.95);
    dum2(:,i) = dum(:,i)./q;
end

figure(56)
boxplot(dum2)
ylim([0 1])


































