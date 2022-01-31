function U01_plot_cortex_Magnitude(cfg,values)
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

    try
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





%% load values that you want to plot, whatever they might be
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




%% build up CortElecLoc

%{
CortElecLoc1 = CortElecLoc;
values1      = Val;

CortElecLoc2 = CortElecLoc;
values2      = Val;

CortElecLoc3 = CortElecLoc;
values3      = Val;

CortElecLoc4 = CortElecLoc;
values4      = Val;

CortElecLoc = [CortElecLoc1 CortElecLoc2 CortElecLoc3 ]
Val         = [values1 values2 values3 ]


indx = find(Val>0.8 | Val<0.2)
Val(indx) = nan;
%}





%% projects an ECoG location to the closest point on mapped cortex
[ a, ~ ] = project2verts( CortElecLoc, cortex.vert );
a = cell2mat(a);


%% pick colormap
close all
%colormap jet;
if cfg.hot==1
    mapp = ([255,247,236;254,232,200;253,212,158;253,187,132;252,141,89;239,101,72;215,48,31;179,0,0;127,0,0])./256;
else
    mapp = ([255,247,251;236,231,242;208,209,230;166,189,219;116,169,207;54,144,192;5,112,176;4,90,141;2,56,88])./256;
end


mapp = [
202,0,32
244,165,130
146,197,222
5,113,176]./256;

colormap(mapp);



% ECoG electrode graphical characteristics
d   = 2.8;              % diameter of decay bigger circle
de  = 2.6;              % diameter of electrode (small circle)
tau = 1.2;              % decay fudge factor
cm  = colormap;         %

%{
% original parameters
d   = 2.4;
de  = 1.2;
tau = 1.2;
cm  = colormap;
%}

% initialize color vector with [1 1 1]'s -- gray (background)
V = zeros(length(cortex.vert),1);
V_color = 1*ones(length(cortex.vert),3);

%% for each electrode, color surrounding vertices according to exponential decay

for v = 1:length(a)
    % find nearby vertices within a given radius
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

% plot ECoG values as a BROAD color on the surface of the brain


figure(201)
Hp = patch('vertices',cortex.vert,'faces',cortex.tri,'FaceVertexCData', V_color,'edgecolor','none','FaceColor','interp',...
    'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5);
axis equal
camlight('headlight','infinite');
axis off;
hold on

%colormap jet;
%{
if cfg.hot==1
    mapp = ([255,247,236;254,232,200;253,212,158;253,187,132;252,141,89;239,101,72;215,48,31;179,0,0;127,0,0])./256;
else
    mapp = ([255,247,251;236,231,242;208,209,230;166,189,219;116,169,207;54,144,192;5,112,176;4,90,141;2,56,88])./256;
end
colormap(mapp);
%}




colormap(mapp);


%set the camera angle parameters
% to get current camera options (get current options, then fix values):
%{
DispCamPos.cp = campos
DispCamPos.cva = camva
DispCamPos.ct = camtarget
DispCamPos.uv = camup
%}


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
    DispCamPos.cp   = [-851.9560 -27.8936 110.2414];
    DispCamPos.cva  =  7.5246;
    DispCamPos.ct   = [-20.3551 6.5690 2.8726];
    DispCamPos.uv   = [0.1316 -0.1026 0.9860];
    
end

set(gca,'CameraPosition',DispCamPos.cp,...
    'CameraTarget',DispCamPos.ct,...
    'CameraViewAngle',DispCamPos.cva,...
    'CameraUpVector',DispCamPos.uv);

%set(gcf,'Position',cfg.pos1)
   set(gcf,'Position',[ 1617         514         625         416])
hold on
%}

%% plot ECoG values as a SINGLE color where the ECoG electrode is supposed to be

figure(301)
%subplot(2,1,2)

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

% plot the electrodes
elec = reshape(cell2mat(CortElecLoc),3,length(CortElecLoc))';
cm = colormap;
%frange = [min(Val), max(Val)];
for e=1:length(CortElecLoc)
    %plot circles with colors
    
    %if isnan(Val(e))==1
        hold on; plot3(elec(e,1), elec(e,2), elec(e,3), 'o', 'markersize', 6, 'color', 'black');
    %else
    %    hold on; plot3(elec(e,1), elec(e,2), elec(e,3), 'o', 'markersize', 6, 'color', getColor( Val(e), [0 1], cm ));
    %end
    
    %plot ECoG number
    if e == 1 || e == 26 || e == 52
        hold on; text(elec(e,1), elec(e,2), elec(e,3), num2str(e),'FontSize',5);
    end
end

% some ECoG electrodes are "incorrectly buried" in the cortex, so make cortex semi transparent
alpha 0.75 

% mark a few electrodes to make sure they are in the correct orientation
%hold on; plot3( elec(1,1),  elec(1,2),  elec(1,3), '.', 'markersize', 4, 'color', 'k');
%hold on; plot3(elec(26,1), elec(26,2), elec(26,3), '.', 'markersize', 4, 'color', 'k');
%hold on; plot3(elec(52,1), elec(52,2), elec(52,3), '.', 'markersize', 4, 'color', 'k');



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
cd /Users/ajorge/Desktop/

saveas(gcf, [Subject, '_magnitude.fig'], 'fig');
%saveas(gcf, [Subject, '_magnitude.png']);
% no luck with saving png or pdf. manually save for now







