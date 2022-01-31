% SCRIPT add STN mni and Cortex mni to avg





%% load
clear all
close all
clc

cd('/Users/ajorge/ajorge/data/ecog_stnpathways/avg/')
load('avg_190818')






%% add STN mni avg

% import table
% import the following document leaddbsdepth2mnihyperdirectv3.xlsx
cd('/Users/ajorge/ajorge/data/ecog_stnpathways')
%manually import table
tablee = leaddbsdepth2mnihyperdirectv4(:,1:10);

for trace = 1:length(avg.data)

    fprintf('trace %5.0f\n',trace)
    
    %find out patient and assoc. row
    onpt    = avg.pt(trace);
    roi_pt  = tablee.subject==onpt;
    
    %find out set and assoc. row
    onsett  = avg.sett(trace);
    roi_sett= tablee.sett==onsett;
    
    %find out final single row
    roi_final = find(roi_pt==1 & roi_sett==1);
    if length(roi_final)~=1
        warning('error')
        pause
    end
    
    avg.STN_mni{trace} = table2array(tablee(roi_final,[8 9 10]));
    
end






%% add Cortex mni to avg

ptoix = load_patient_list_STNhyperdirect(2)
 
for ptoi = ptoix

    tofolder = sprintf('/Users/ajorge/ajorge/data/ecog_stnpathways/freesurfer_preop_ElecLoc/DBS%4.0f_preop_ElecLoc/MNI_on_cortex',ptoi);
    cd(tofolder);
    load('CortElecLocL.mat');

    %all the traces associated with ptoi
    inx = find(avg.pt == ptoi);

    %for each trace, associate with an ECoG mni location
    for i = 1:length(inx)

        ecogn = avg.electrode_n(inx(i));
        if ecogn >63
            ecogn=ecogn-1;
        end
        
        try
            avg.ECOG_mni(inx(i),:) = CortElecLoc{ecogn};
        catch
            avg.ECOG_mni(inx(i),:) = [nan nan nan];
        end

    end


end


%% plot all ecogs locations onto cortex (last updated 2021 12 12) (Figure 1 making)


figure(201)

%load MNI brain
cd('/Users/ajorge/ajorge/data/ecog_stnpathways')
load('cortex_MNI.mat')

%plot MNI brain
atlas=2; % specify Desikan-Killiany atlas
V = zeros(length(BS1.Vertices),3); % initialize color matrix   6 4 3 43 22 44 45 40
for region=1:length(BS1.Atlas(atlas).Scouts)
    if region == 31 || region == 57 || region == 55 || region == 15 || region == 67 || region == 25 || region == 29 || region == 51 
        V(BS1.Atlas(atlas).Scouts(region).Vertices,:) = repmat(BS1.Atlas(atlas).Scouts(region).Color,length(BS1.Atlas(atlas).Scouts(region).Vertices),1);
    else
        V(BS1.Atlas(atlas).Scouts(region).Vertices,:) = repmat([0.7 0.7 0.7],length(BS1.Atlas(atlas).Scouts(region).Vertices),1);
    end
end


figure(301)
patch('vertices',BS1.Vertices,'faces',BS1.Faces,'FaceVertexCData',V,'edgecolor','none','FaceColor','interp');
axis equal
%camlight('headlight','infinite');
camlight('left','infinite');
fh(1)=gcf;
axis off;
hold on

for i = 1:length(avg.sett)
    plot3(avg.ECOG_mni(i,1), avg.ECOG_mni(i,2), avg.ECOG_mni(i,3),'ko')
    plot3(avg.ECOG_mni(i,1), avg.ECOG_mni(i,2), avg.ECOG_mni(i,3),'w.', 'MarkerSize',10)
end



view([-90 0])



%% save
avg2 = avg;
clear avg;
avg = avg2;
cd('/Users/ajorge/ajorge/data/ecog_stnpathways/avg/')
save('avg_190818','avg')









