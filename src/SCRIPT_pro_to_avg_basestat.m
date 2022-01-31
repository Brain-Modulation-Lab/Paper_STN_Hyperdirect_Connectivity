% SCRIPT_pro_to_avg

% takes the entire 30 trials and averages them into one 
%
% input:
% pro.trial{30:90}(choi,:) contains 30 trials
%
% output:
% avg.data{trace} contains averaged data for the 30 trials below
%
% 

clear all
close all
clc


%% Get avg
% Save averaged, smoothed traces in avg along with electrode info



%ptoi = load_patient_list_STNhyperdirect(2);

%patients to load 
%ptoi = [3005        3008        3010        3011        3012        3016        3017        3018        3019        3023   3024        3026        3027        3028        3029        3030 3032];
%[4079 4081 4085]

ptoi = [4079 4081 4085]

%QUICK avg clear
%{
cd('/Users/ajorge/Desktop')
avg = [];
save('avg','avg')
trace=1
%}

tic
for pt=ptoi
    
    switch pt
        case 3005
            setts = [1 2 3 4 5 6];
        case 3006
            setts = [1 2 3];
        case 3008
            setts = [1 2 3];
        case 3010
            setts = [1 2 3 4 5 6];
        case 3011
            setts = [1 2 3 4];
        case 3012
            setts = [1 2 3];
        case 3016
            setts = [4 5 6];
        case 3017
            setts = [1 2];
        case 3018
            setts = [2 3 4];
        case 3019
            setts = [1 2 3];
        case 3023
            setts = [1 2 3 4 5 6 7 8 9];
        case 3024
            setts = [1 2 3 4 5 6];
        case 3025
            setts = [1 2 3 4 5 6];
        case 3026
            setts = [1 2 3 4 5 6];
        case 3027
            setts = [1 2 3 4 5 6];
        case 3028
            setts = [2 4 6 8 10 12];
        case 3029
            setts = [2 4 5 8 10 12];
        case 3030
            setts = [2 4 6 8 10 12];
        case 3032
            setts = [1 2 3 4 5 6];
            
        case 4079 
            setts = [1 2 3];
        case 4081 
            setts = [1 2 3 4 5 6 7 8 9];
        case 4085 
            setts = [1 2 3 4 5 6];
            
            
            
        otherwise
            warning('need to specify setts')
            pause
    end
    
    
    
%calculate baseline t-statistic once
[tstatb_975, tstatb_025] = pro2avg_baseline(pt);        

%add more traces to avg (this will also calculate the baseline t-statistic)
for sett = setts
    
    clear avg
    cd('/Users/ajorge/Desktop')
    load('avg')
    if size(avg)>0
        trace = size(avg.data,2)+1;
    else
        trace = 1;
    end
    
    avg = pro2avg(pt,sett,avg,trace,tstatb_975,tstatb_025);
    %avg = pro2avg_bipolar_v1(pt,sett,avg,trace,tstatb_975,tstatb_025);
    %avg = pro2avg_orthodromic(pt,sett,avg,trace,tstatb_975,tstatb_025);
    
    
    
    % save avg
    cd('/Users/ajorge/Desktop')
    save('avg','avg')

end

fprintf('time so far %4.0f\n',toc/60)

end


% how many traces per brodmann area
figure(401)
inx = find(avg.brodmann<200);
histogram(avg.brodmann(inx),200)

unique(avg.pt)

toc



%% manually save
% grab avg.mat from desktop and manually rename it

%% Load avg
cd('/Users/ajorge/Desktop')
load('avg','avg')

unique(avg.pt)























