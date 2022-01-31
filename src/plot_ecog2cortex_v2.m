function plot_ecog2cortex_v2(cfg,values)






%%
clc

%load patient ECOG mni coordinates
tofolder = sprintf('/Users/ajorge/ajorge/data/ecog_stnpathways/freesurfer_preop_ElecLoc/DBS%4.0f_preop_ElecLoc/MNI_on_cortex',cfg.ptoi);

cd(tofolder);

if cfg.ptoi==4081
    load('CortElecLoc_MNI.mat'); 
    CortElecLoc = CortElecLoc_MNI;
elseif cfg.ptoi==4085
    load('CortElecLoc.mat'); 
else
    load('CortElecLocL.mat');
end

figure(301)
%plot MNI coordinates
for i = 1:63
    plot3(CortElecLoc{i}(1), CortElecLoc{i}(2), CortElecLoc{i}(3),'k.','MarkerSize',25)
    
    if cfg.enumerate==1
        if i == 1 || i == 5 || i == 10 || i == 15 || i == 20 
            plot3(CortElecLoc{i}(1), CortElecLoc{i}(2), CortElecLoc{i}(3),'r.','MarkerSize',25)
            namee = sprintf('%2.0f',i);
            text(CortElecLoc{i}(1), CortElecLoc{i}(2), CortElecLoc{i}(3),namee,'Color','red','FontSize',14)
        end
    end
    
    hold on
end
try
for i = 64:126
    plot3(CortElecLoc{i}(1), CortElecLoc{i}(2), CortElecLoc{i}(3),'k.','MarkerSize',25)
    
    if cfg.enumerate==1
        if i == 1+63 || i == 5+63 || i == 10+63 || i == 15+63 || i == 20+63 
            plot3(CortElecLoc{i}(1), CortElecLoc{i}(2), CortElecLoc{i}(3),'r.','MarkerSize',25)
            namee = sprintf('%2.0f',i);
            text(CortElecLoc{i}(1), CortElecLoc{i}(2), CortElecLoc{i}(3),namee,'Color','red','FontSize',14)
        end
    end
    
    hold on
end  
catch
end







if cfg.plotbrain==1
    %load MNI brain
    cd('/Users/ajorge/ajorge/data/ecog_stnpathways')
    load('cortex_MNI.mat')

    %plot MNI brain
    atlas=2; % specify  atlas (I am using atlas 2, Witek uses 3 (as of Nov 2019)
    V = zeros(length(BS1.Vertices),3); % initialize color matrix
    
    %color tested cortical regions (e.g. M1, S1...)
    for region=1:length(BS1.Atlas(atlas).Scouts)
        
       % if region == 29 || region == 25 || region == 31 || region == 57 || region == 15 || region == 55 || region == 67 || region == 51 
            V(BS1.Atlas(atlas).Scouts(region).Vertices,:) = repmat(BS1.Atlas(atlas).Scouts(region).Color,length(BS1.Atlas(atlas).Scouts(region).Vertices),1);
        %else
           % V(BS1.Atlas(atlas).Scouts(region).Vertices,:) =  repmat([211 211 211]./256,length(BS1.Atlas(atlas).Scouts(region).Vertices),1);
       % end
    
    end
    %gray out untested cortical regions
    
    
    figure(301)
    
    
    %patch('vertices',BS1.Vertices,'faces',BS1.Faces,'FaceVertexCData',V,'edgecolor','none','FaceColor','interp');
    patch('vertices',BS1.Vertices,'faces',BS1.Faces,'FaceVertexCData',V,'edgecolor','none','FaceColor','interp','facelighting', 'gouraud')%, 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5);
    
    axis equal
    camlight('headlight','infinite');
    fh(1)=gcf;
    axis off;
    view([-90 0])
    set(gcf,'Position',[1196         487         959         688])
    
    

end












