function CortElecLoc2spreadsheet(pt)



%% load patient-specific electrode locations
if cfg.pt == 3019
    clear pathtofile
    pathtofile = sprintf('/Volumes/Nexus/DBS/DBS%4.0f/Anatomy/FreeSurfer/preop/Electrode_Locations/archive/Dengyu',cfg.pt);
    cd(pathtofile)
    load('CortElecLocL_SMC_eq.mat')
else
    clear pathtofile
    pathtofile = sprintf('/Volumes/Nexus/DBS/DBS%4.0f/Anatomy/FreeSurfer/preop/Electrode_Locations',cfg.pt);
    cd(pathtofile)
    load('CortElecLocL_eq.mat')
end

%only proceed with first 1:63 electrodes for now
for chann = 1:126
    CortElecLoc2{chann} = CortElecLoc{chann};
end
clear CortElecLoc
CortElecLoc = CortElecLoc2;


%some patients have "inverted" electrode arrays (eg. LFP 1 would be on LFP 21, and vice versa)
if cfg.pt == 3017
    CortElecLoc(1:21)  = flip(CortElecLoc(1:21));
    CortElecLoc(22:42) = flip(CortElecLoc(22:42)); 
    CortElecLoc(43:63) = flip(CortElecLoc(43:63)); 
end

%% assign values

for i = 1:63
    
    mat(i,1) = CortElecLoc{i}(1);
    mat(i,2) = CortElecLoc{i}(2);
    mat(i,3) = CortElecLoc{i}(3);
    
end

mat(64,1) = NaN;
mat(64,2) = NaN;
mat(64,3) = NaN;


for i = 65:126

    mat(i,1) = CortElecLoc{i}(1);
    mat(i,2) = CortElecLoc{i}(2);
    mat(i,3) = CortElecLoc{i}(3);

end

for i = 127:163

    mat(i,1) = NaN;
    mat(i,2) = NaN;
    mat(i,3) = NaN;

end



% MANUALLY copy "mat" to the spreadsheet result



