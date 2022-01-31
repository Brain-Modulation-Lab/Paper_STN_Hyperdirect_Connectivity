function [array1, array2] = ecog_cortexlocs_auto(pt,plotornot)

% input:  pt number (only the number from the DBSxxxx format and electrode
% output: 21x3 matrix (x2 arrays) containing channel's Broadmann area location




%% load cortex mni
load('/Users/ajorge/ajorge/analysis/ecog_stnpathways/brodmann_and_ecogs/cortex_MNI.mat');

%% specify atlas
% 2 = original atlas used by Dengyu
% 3 = Desikan-Killiany atlas (prefered by Dengyu)
% 4 = Brodmann
atlasnum=2; % specify Desikan-Killiany atlas

%% only select gyrus
Ind_Gy_vert = [];
ind_region = [];
for i = 1:length(BS1.Atlas(atlasnum).Scouts)    
    if any(BS1.Atlas(atlasnum).Scouts(i).Label == 'G')
        ind_region = [ind_region,i];
        Ind_Gy_vert = [Ind_Gy_vert, BS1.Atlas(atlasnum).Scouts(i).Vertices];
    end    
end

Ind_Gy_vert = sort(Ind_Gy_vert);
Gyrus = BS1.Vertices(Ind_Gy_vert,:);

%% match patient specific MNI coordinates to atlas index

clear min_val idx CortElecLoc_MNI filenamee
clc

% load patient's CortElecLoc (either MNI or eq)
%filenamee = sprintf('/Volumes/Nexus/DBS/DBS%4.0f/Anatomy/FreeSurfer/preop/Electrode_Locations/*_eq_Session3.mat',pt);
%filenamee = sprintf('/Volumes/Nexus/DBS/DBS%4.0f/Anatomy/FreeSurfer/preop/Electrode_Locations/*_MNI.mat',pt);
%Dir_Elec = dir([filenamee]);
%load([Dir_Elec.folder,filesep,Dir_Elec.name]);

filenamee = sprintf('/Users/ajorge/ajorge/data/ecog_stnpathways/freesurfer_preop_ElecLoc/DBS%4.0f_preop_ElecLoc/MNI_on_cortex/CortElecLoc_MNI.mat',pt);
load(filenamee)

if pt == 3019 || pt == 3022
    CortElecLoc_MNI = CortElecLoc;
end


%find the gyrus with the minimum distance to ECoG electrode
[min_val,idx] = min(pdist2(cell2mat(CortElecLoc_MNI'),Gyrus),[],2);
%[min_val,idx] = min(pdist2((CortElecLoc),Gyrus),[],2);


%how far away can the ecog electrode be from cortex before we mark it
tolerance = 7;
if max(min_val)>tolerance
    warning('distance between ecog and gyrus is >6mm')
    pause
end

elec2mnicort_indx = [Ind_Gy_vert(idx)',min_val];


%% associate mni coordinates (in the form on a single idx number) with atlas region
% look for your index in all cortical regions
% prefrontal
clear roi
roi = BS1.Atlas(atlasnum).Scouts(:);
roi_idx = [];
roi_lab = [];
for i = 1:length(roi)
    roi_idx = [roi_idx, roi(i).Vertices];
    roi_lab = cat(1,roi_lab, i'*ones(length( roi(i).Vertices),1));
end
roi_lab = roi_lab';


for k=1:length(elec2mnicort_indx)
    areaindx(k) = roi_lab(find(roi_idx==elec2mnicort_indx(k,1)));
    areaname{k} = BS1.Atlas(atlasnum).Scouts(areaindx(k)).Label;
    areacolor(k,:) = BS1.Atlas(atlasnum).Scouts(areaindx(k)).Color;
end


%% plot atlas
V = zeros(length(BS1.Vertices),3); % initialize color matrix
for region=1:length(BS1.Atlas(atlasnum).Scouts)
    V(BS1.Atlas(atlasnum).Scouts(region).Vertices,:) = repmat(BS1.Atlas(atlasnum).Scouts(region).Color,length(BS1.Atlas(atlasnum).Scouts(region).Vertices),1);
end


if plotornot == 1 
    figure(450); 
    patch('vertices',BS1.Vertices,'faces',BS1.Faces,'FaceVertexCData',V,'edgecolor','none','FaceColor','interp','FaceAlpha',0.7);
    axis equal
    camlight('headlight','infinite');
    fh(1)=gcf;
    axis off;
end


%% plot ecog electrodes on atlas

if plotornot == 1 
    hold on
    figure(450)
    clc

    for elecnum = 1:length( CortElecLoc_MNI)
        clear dummy1
        dummy1 = CortElecLoc_MNI{elecnum};
        plot3(dummy1(1), dummy1(2), dummy1(3), '*', 'MarkerSize',10,'Color', areacolor(elecnum,:))
        hold on
    end

    view([-99 15])
    set(gcf,'Position',[  914          83        1356        1049])
end

%% give Broadmann areas

% 3 = primary sensory cortex (including 1,2 and 3 for simplicity purposes)
% 4 = primary motor cortex
% 6 = premotor
% 43 = gustatory


% 21 = 22 = 41 = 42 = STG
% 44 = pars opercularis
% 43 = pars triangularis

% 91 = difficult to tell between two areas, could run downstream analysis with value from electrode before or electrode after
% 99 = difficult to tell the area (based on the ecog localization, might as well code these as NaNs)


%% first array
clear array1_a

array1_a = nan*ones(63,1);




if length( CortElecLoc_MNI)>=63
    ul = 63;
else
    ul = length( CortElecLoc_MNI);
end

for m=1:ul
    
    if min_val(m)<tolerance
        if strcmp('G_temp_sup-Lateral L',areaname(m)) == 1
            array1_a(m) = 22;
        
        elseif strcmp('G_and_S_subcentral L',areaname(m)) == 1
            array1_a(m) = 43; %subcentral gyrus = gustatory, semantic analysis
        
        elseif strcmp('G_postcentral L',areaname(m)) == 1
            array1_a(m) = 3;
        elseif strcmp('G_precentral L',areaname(m)) == 1
            array1_a(m) = 4;
        elseif strcmp('G_front_middle L',areaname(m)) == 1
            array1_a(m) = 6;

        elseif strcmp('G_front_inf-Opercular L',areaname(m)) == 1
            array1_a(m) = 44;
        elseif strcmp('G_front_inf-Triangul L',areaname(m)) == 1
            array1_a(m) = 45;
            
        elseif strcmp('G_pariet_inf-Supramar L',areaname(m)) == 1
            array1_a(m) = 40; %supramarginal = secondary sensory cortex
            
        else
            array1_a(m) = 999;
        end
    else
        array1_a(m) = 998;
    end
    
end



array1 = array1_a;
%{
try
    array1_a = reshape(array1_a,[21,3]);
    array1 = array1_a';
catch
    array1 = array1_a;
end
%}

%figure
%imagesc(array1)



%% second array
clear array1_a
array1_a = nan*ones(63,1);



try
for m=1:63
    
    n=m+63;
    
    if min_val(m)<tolerance
        if strcmp('G_temp_sup-Lateral L',areaname(n)) == 1
            array1_a(m) = 22;
        elseif strcmp('G_and_S_subcentral L',areaname(n)) == 1
            array1_a(m) = 43;
        elseif strcmp('G_postcentral L',areaname(n)) == 1
            array1_a(m) = 3;
        elseif strcmp('G_precentral L',areaname(n)) == 1
            array1_a(m) = 4;
        elseif strcmp('G_front_middle L',areaname(n)) == 1
            array1_a(m) = 6;

        elseif strcmp('G_front_inf-Opercular L',areaname(n)) == 1
            array1_a(m) = 44;
        elseif strcmp('G_front_inf-Triangul L',areaname(n)) == 1
            array1_a(m) = 45;
        else
            array1_a(m) = 999;
        end
    else
        array1_a(m) = 999;
    end
end


%{
array1_a = reshape(array1_a,[21,3]);
array2 = array1_a';
%}
array2 = array1_a;


catch
    array2=nan*ones(1,63)';
end





%figure
%imagesc(array2)








