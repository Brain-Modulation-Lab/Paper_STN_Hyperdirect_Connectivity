function plot_ecog2cortex(cfg,values)
%
% dependencies:
% getColor.m
% projects2verts.m
%
% input:
% cfg.pt = xxxx e.g. 3012
% cfg.mni = 0 (plot native cortex) OR = 1 (plot MNI cortex)
% cfg.hot = 1 (hot color scheme) OR = 0 (cold color scheme)
% values (1x63 array) or (1x126 array)
%
% output:
% plot values onto brain cortex
% check that ECoG channels are correctly numbered
%
% ajorge from witold's code sep 2018




%% load patient-specific cortex AND patient-specific e- loc


if cfg.mni == 0
    Subject = sprintf('DBS%4.0f',cfg.pt);

    %load patient-specific cortex
    clear pathtofile
    pathtofile = sprintf('/Volumes/Nexus/DBS/DBS%4.0f/Anatomy/FreeSurfer/preop',cfg.pt);
    cd(pathtofile)
    load('cortex_indiv.mat')

    try`
        try
            %load patient-specific electrode locations
            if cfg.pt == 3019
                clear pathtofile
                pathtofile = sprintf('/Volumes/Nexus/DBS/DBS%4.0f/Anatomy/FreeSurfer/preop/Electrode_Locations/archive/Dengyu',cfg.pt);
                cd(pathtofile)
                load('CortElecLocL_SMC_eq.mat')
            elseif cfg.pt == 3022
                clear pathtofile
                pathtofile = sprintf('/Volumes/Nexus/DBS/DBS%4.0f/Anatomy/FreeSurfer/preop/Electrode_Locations/MNI',cfg.pt);
                cd(pathtofile)
                load('CortElecLoc_eq_Session3.mat')    
            elseif cfg.pt == 3023
                clear pathtofile
                pathtofile = sprintf('/Volumes/Nexus/DBS/DBS%4.0f/Anatomy/FreeSurfer/preop/Electrode_Locations/',cfg.pt);
                cd(pathtofile)
                load('CortElecLoc_MNI.mat') 
                CortElecLoc = CortElecLoc_MNI;
            elseif cfg.pt == 3025
                clear pathtofile
                pathtofile = sprintf('/Volumes/Nexus/DBS/DBS%4.0f/Anatomy/FreeSurfer/preop/Electrode_Locations',cfg.pt);
                cd(pathtofile)
                load('CortElecLoc_eq_Session3.mat') 

            else
                clear pathtofile
                pathtofile = sprintf('/Volumes/Nexus/DBS/DBS%4.0f/Anatomy/FreeSurfer/preop/Electrode_Locations',cfg.pt);
                cd(pathtofile)
                load('CortElecLocL_eq.mat')
            end

        catch
            clear pathtofile
            pathtofile = sprintf('/Volumes/Nexus/DBS/DBS%4.0f/Anatomy/FreeSurfer/preop/Electrode_Locations',cfg.pt);
            cd(pathtofile)
            load('CortElecLocL.mat')
        end
    catch
        clear pathtofile
        pathtofile = sprintf('/Volumes/Nexus/DBS/DBS%4.0f/Anatomy/FreeSurfer/Electrode_Locations',cfg.pt);
        cd(pathtofile)
        load('CortElecLocL.mat')
    end



    clear CortElecLoc2
    %only proceed with first 1:63 electrodes for now
    if cfg.pt == 3010 || cfg.pt == 3018 || cfg.pt == 3020 || cfg.pt>=3024
        for chann = 1:126
            CortElecLoc2{chann} = CortElecLoc{chann};
        end
    elseif cfg.pt == 3025
        for chann = 1:114
            CortElecLoc2{chann} = CortElecLoc(chann,:);
        end
    else
        for chann = 1:63
            CortElecLoc2{chann} = CortElecLoc{chann};
        end
    end
    clear CortElecLoc
    CortElecLoc = CortElecLoc2;


    %some patients have "inverted" electrode arrays (eg. LFP 1 would be on LFP 21, and vice versa)
    if cfg.pt == 3017
        CortElecLoc(1:21)  = flip(CortElecLoc(1:21));
        CortElecLoc(22:42) = flip(CortElecLoc(22:42)); 
        CortElecLoc(43:63) = flip(CortElecLoc(43:63)); 
        warning('check the other ecog array')
        %pause
    end

    
    
    
    
    
    
elseif cfg.mni == 1
    %load MNI cortex
    clear pathtofile
    cd('/Users/ajorge/ajorge/analysis/ecog_stnpathways/brodmann_and_ecogs');
    load('cortex_MNI.mat')
    cortex.vert = BS1.Vertices;
    cortex.tri = BS1.Faces;
    
    %load patient's MNI coordinates
    Subject = sprintf('DBS%4.0f',cfg.pt);
    clear pathtofile
    pathtofile = sprintf('/Volumes/Nexus/DBS/DBS%4.0f/Anatomy/FreeSurfer/preop/Electrode_Locations/MNI_on_cortex',cfg.pt);
    cd(pathtofile)
    load('CortElecLocL.mat')
    
end





% load values that you want to plot, whatever they might be
Val = values;

%Val needs to be a 63x1 array
if length(Val) == 126
else
    warning('inconsistent size')
    
end

% normalize as desired (Val numbers need to be between 0 and 1)
ValueRange(1)   = nanmin(nanmin(Val))-0.00001;  %having a hard time plotting 0, seems that it needs to be (0,1), not [0,1]
ValueRange(2)   = nanmax(nanmax(Val));

%ValueRange      = [-nanmax(abs(ValueRange)), nanmax(abs(ValueRange))];
Val             = (Val - ValueRange(1))/(ValueRange(2)-ValueRange(1)); 



%% projects an ECoG location to the closest point on mapped cortex
[ a, ~ ] = project2verts( CortElecLoc, cortex.vert );
a = cell2mat(a);

%% pick colormap

% ECoG electrode graphical characteristics
d   = 2.8;
de  = 2.6;
tau = 1.2;
cm  = colormap;

V = zeros(length(cortex.vert),1);
V_color = 1*ones(length(cortex.vert),3);

%{
for v = 1:length(a)
    aeconn = find(pdist2(cortex.vert(a(v),:), cortex.vert)<=de);
    aconn = find(pdist2(cortex.vert(a(v),:), cortex.vert)>de & pdist2(cortex.vert(a(v),:), cortex.vert)<=d);
    dconn = pdist2(cortex.vert(a(v),:), cortex.vert(aconn,:));
%     V(aeconn) = V(aeconn)+Val(v);
%     V(aconn) = V(aconn)+Val(v)*exp(-(dconn-de)/tau)';
    V(aeconn) = arrayfun(@(x) max(Val(v), x), V(aeconn));
    V(aconn) = arrayfun(@max, Val(v)*exp(-(dconn-de)/tau)', V(aconn));
end

for v = find(V~=0)'
    V_color(v,:) = getColor(V(v), [0 1], cm);
end
%}


%% plot ECoG values as a SINGLE color where the ECoG electrode is supposed to be

figure(301)
%subplot(2,1,2)

if cfg.last==1
    %{
    % plot the cortex
    Hp = patch('vertices',cortex.vert,'faces',cortex.tri,'FaceVertexCData', V_color,'edgecolor','none','FaceColor','interp',...
       'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5);
    axis equal
    camlight('headlight','infinite');
    axis off;
    hold on

    %colormap jet;
    if cfg.hot==1
        mapp = ([255,247,236;254,232,200;253,212,158;253,187,132;252,141,89;239,101,72;215,48,31;179,0,0;127,0,0])./256;
    else
        mapp = ([255,247,251;236,231,242;208,209,230;166,189,219;116,169,207;54,144,192;5,112,176;4,90,141;2,56,88])./256;
    end
    colormap(mapp);
    %}
    
    
    %% load MNI brain
    load('/Users/ajorge/ajorge/analysis/ecog_stnpathways/brodmann_and_ecogs/cortex_MNI.mat');

    %% specify atlas
    % 2 = original atlas used by Dengyu
    % 3 = Desikan-Killiany atlas (prefered by Dengyu)
    % 4 = Brodmann
    atlas=2; % specify Desikan-Killiany atlas

    %% associate mni coordinates (in the form on a single idx number) with atlas region
    % look for your index in all cortical regions
    % prefrontal
    clear roi
    roi = BS1.Atlas(atlas).Scouts(:);
    roi_idx = [];
    roi_lab = [];
    for i = 1:length(roi)
        roi_idx = [roi_idx, roi(i).Vertices];
        roi_lab = cat(1,roi_lab, i'*ones(length( roi(i).Vertices),1));
    end
    roi_lab = roi_lab';

    areaindx = roi_lab(find(roi_idx==100803));
    areaname = BS1.Atlas(atlas).Scouts(areaindx).Label


    %% plot atlas
    V = zeros(length(BS1.Vertices),3); % initialize color matrix
    for region=1:length(BS1.Atlas(atlas).Scouts)
        V(BS1.Atlas(atlas).Scouts(region).Vertices,:) = repmat(BS1.Atlas(atlas).Scouts(region).Color,length(BS1.Atlas(atlas).Scouts(region).Vertices),1);
    end

    figure(301)
    patch('vertices',BS1.Vertices,'faces',BS1.Faces,'FaceVertexCData',V,'edgecolor','none','FaceColor','interp');
    axis equal
    camlight('headlight','infinite');
    fh(1)=gcf;
    axis off;
    
    
    
end

% plot the electrodes
elec = reshape(cell2mat(CortElecLoc),3,length(CortElecLoc))';
cm = colormap;
%frange = [min(Val), max(Val)];
for e=1:length(CortElecLoc)

    hold on; plot3(elec(e,1), elec(e,2), elec(e,3), '.', 'markersize', 20, 'color', 'black');

end

% some ECoG electrodes are "incorrectly buried" in the cortex, so make cortex semi transparent
alpha 0.75 

% set the camera angle parameters
if cfg.pt == 3017
    DispCamPos.cp   = [  -594.2715   48.8825  115.2896];
    DispCamPos.cva  =  7.1057;
    DispCamPos.ct   = [ -18.1868 72.7562 40.9107];
    DispCamPos.uv   = [    0     0     1];
elseif cfg.pt == 3018
    DispCamPos.cp   = [-851.9560 -27.8936 110.2414];
    DispCamPos.cva  =  7.5246;
    DispCamPos.ct   = [-20.3551 6.5690 2.8726];
    DispCamPos.uv   = [0.1316 -0.1026 0.9860];
else
    DispCamPos.cp   = [-876.9612   58.2398  166.9879];
    DispCamPos.cva  =  7.5246;
    DispCamPos.ct   = [-20.3551 6.5690 2.8726];
    DispCamPos.uv   = [0.1316 -0.1026 0.9860];

end

set(gca,'CameraPosition',DispCamPos.cp,...
    'CameraTarget',DispCamPos.ct,...
    'CameraViewAngle',DispCamPos.cva,...
    'CameraUpVector',DispCamPos.uv);

%set(gcf,'Position',cfg.pos1)
set(gcf,'Position',[     1617          19         618         478])
%}


%% Save figures
if cfg.last==1
    cd /Users/ajorge/Desktop/
    saveas(gcf, [Subject, '_magnitude.fig'], 'fig');
    %saveas(gcf, [Subject, '_magnitude.png']);
    % no luck with saving png or pdf. manually save for now
end






