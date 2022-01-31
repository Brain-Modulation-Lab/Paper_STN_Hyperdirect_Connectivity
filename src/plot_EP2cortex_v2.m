function plot_EP2cortex(subject_list, val_mat_list, xe_indxed)
%
%
% ajorge march 2021
%
%




%% load MNI brain
load('/Users/ajorge/ajorge/analysis/ecog_stnpathways/brodmann_and_ecogs/cortex_MNI.mat');

%% specify atlas
% 2 = original atlas used by Dengyu
% 3 = Desikan-Killiany atlas (prefered by Dengyu)
% 4 = Brodmann
atlas=2; % specify Desikan-Killiany atlas



%% ---------------------------------------------------------------------------
% get tt (vertices) and V
%{
%close all
clc



%change subject list to string subject list
for i = 1:length(subject_list)
    subject_list2{i} = sprintf('DBS%4.0f',subject_list(i));
end

tt = val2cortex(subject_list2, val_mat_list, 'median');


% process vertices
%ValueRange(1)   = nanmin(tt);  %having a hard time plotting 0, seems that it needs to be (0,1), not [0,1]
%ValueRange(2)   = nanmax(tt);
%tt_scaled = (tt - ValueRange(1))/(ValueRange(2)-ValueRange(1));

% assign color mapp
mapp = ([255,247,236;254,232,200;253,212,158;253,187,132;252,141,89;239,101,72;215,48,31;179,0,0;127,0,0])./256;

mapp = flipud([178,24,43
214,96,77
244,165,130
253,219,199
247,247,247
209,229,240
146,197,222
67,147,195
33,102,172])./256;

% preassign black to all vertices
V = ones(length(tt),3);

% divide vertices values into mapp values 
% h1 = histogram(tt,length(mapp));         %9 colors

howmany = 1:length(mapp);
for i = 1:length(howmany)
     indexx{i} = find(tt < prctile(tt,(100/9) * i ) & tt >= prctile(tt, (100/9) * (i-1)));
     
     if i == 9
         indexx{i} = find(tt <= prctile(tt,(100/9) * i ) & tt >= prctile(tt, (100/9) * (i-1)));
     end
     
end



% assign each binned vertix above a color from mapp
for vert_id = 1:length(tt)
    if isnan(tt(vert_id))
    else
        color_idx = find(cell2mat(cellfun(@(x) ismember(vert_id,x),indexx,'UniformOutput',0)));
        
        V(vert_id,:) = mapp(color_idx,:);
    end
end

%}





%% or don't get V and plot raw circle for voltages

V = zeros(length(BS1.Vertices),3); % initialize color matrix
    for region=1:length(BS1.Atlas(atlas).Scouts)
        %V(BS1.Atlas(atlas).Scouts(region).Vertices,:) = repmat(BS1.Atlas(atlas).Scouts(region).Color,length(BS1.Atlas(atlas).Scouts(region).Vertices),1);
        V(BS1.Atlas(atlas).Scouts(region).Vertices,:) = repmat([200 200 200]/.256,length(BS1.Atlas(atlas).Scouts(region).Vertices),1);

    end

    
%---------------------------------------------------------------------------
% plot brain
figure(444)
cd('/Users/ajorge/ajorge/analysis/ecog_stnpathways/brodmann_and_ecogs');
load('cortex_MNI.mat')
cortex.vert = BS1.Vertices;
cortex.tri = BS1.Faces;

% plot the cortex
Hp = patch('vertices',cortex.vert,'faces',cortex.tri,'FaceVertexCData', V,'edgecolor','none','FaceColor','interp','facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5);
axis equal
camlight('headlight','infinite');
axis off;
hold on

% some ECoG electrodes are "incorrectly buried" in the cortex, so make cortex semi transparent
alpha 0.75 
%alpha 1

% mark a few electrodes to make sure they are in the correct orientation
%hold on; plot3( elec(1,1),  elec(1,2),  elec(1,3), '.', 'markersize', 4, 'color', 'k');
%hold on; plot3(elec(26,1), elec(26,2), elec(26,3), '.', 'markersize', 4, 'color', 'k');
%hold on; plot3(elec(52,1), elec(52,2), elec(52,3), '.', 'markersize', 4, 'color', 'k');

% set the camera angle parameters
DispCamPos.cp   = [  -594.2715   48.8825  115.2896];
DispCamPos.cva  =  7.1057;
DispCamPos.ct   = [ -18.1868 72.7562 40.9107];
DispCamPos.uv   = [    0     0     1];

view([-100 7])

%}




%% plot voltage values on top




clear colorss
colorss = hot(256);
colorss = flipud(colorss);
%colorss = [  colorss(100:1:256,:) ];
%colorss = [  colorss(1:1:66,:);  colorss(171:1:256,:)   ];

%colorss = diverging_map([linspace(0,1,256)],[0.230, 0.299, 0.754],[0.706, 0.016, 0.150]);
%colorss = diverging_map([linspace(0,1,256)],[209,229,240]./256,[215,48,31]./256);

%hh = logspace(0,1,256)/10

%colorss = [  colorss(150:1:256,:) ];
%colorss = [  colorss(1:10:64,:);  colorss(65:5:126,:); colorss(127:2:190,:); colorss(190:1:256,:)   ];

figure(444)
clear inx9

%organize from low voltage to high voltage
%[val_mat_list_2, inxs] = sort(val_mat_list);

for i = 1:length(xe_indxed)
    
   
    
    inx8 = round((length(colorss)-1)/1*val_mat_list(i)+1);
        
    %wolfram alpha to fit a log funciton on {{0.25,30},{.50,180},{.75,180},{1.0,200}}
    %if val_mat_list(i)<0; x=
    %inx8 = round(61*log(0.14*(abs(val_mat_list(i)*1000))));
    
    
    if inx8<1; inx8=1; end
    if inx8>length(colorss); inx8=length(colorss); end
    
    
    inx9(i) = inx8;
    
    %plot(xe_indxed(i,1),xe_indxed(i,2),'r.')
    %plot3(xe_indxed(i,1),xe_indxed(i,2),xe_indxed(i,3),'o','MarkerSize',5,'MarkerFaceColor',colors(inx7,:),'MarkerEdgeColor',colors(inx7,:));
    plot3(xe_indxed(i,1),xe_indxed(i,2),xe_indxed(i,3),'o','MarkerSize',5,'MarkerFaceColor',colorss(inx8,:),'MarkerEdgeColor',colorss(inx8,:));
    hold on
   
  

end


figure(445)
imagesc(1:length(colorss))
colormap(colorss)
xticks([1 256])
set(gca,'xticklabel',[0 1])
set(gca,'yticklabel',{[]})
xlabel([0 1])
xlabel('normalized EP voltage')

figure(446)
histogram(inx9)
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})



%%




x = -1 + 2*rand(1,10000);
y = sqrt(2)*erfinv(x);

figure(1)
histogram(x)

figure(2)
h = histogram(val_mat_list)


x2 = 1/2 * erfc(-val_mat_list/sqrt(2))
figure(3)
histogram(x2)








%% oficial histogram for this plot

clc

%oficial histogram
val_mat_list(find(val_mat_list<0))=[];
figure(445)
h1 = histogram(val_mat_list,50);
figure(447)
b = bar(h1.BinEdges(1:(end-1)), h1.Values, 0.9);
colorss = ([5,113,176; 146,197,222; 244,165,130; 202,0,32])./256;




b.FaceColor = 'flat';

iq(1) = quantile(val_mat_list,0.25);
iq(2) = quantile(val_mat_list,0.50);
iq(3) = quantile(val_mat_list,0.75);
iq(4) = quantile(val_mat_list,1.00);


for i = 1:length(h1.BinEdges(1:(end-1)))
    if h1.BinEdges(i) < iq(1)
        b.CData(i,:) = [colorss(1,:)];
    elseif h1.BinEdges(i) <= iq(2)
        b.CData(i,:) = [colorss(2,:)];
    elseif h1.BinEdges(i) <= iq(3)
       b.CData(i,:) = [colorss(3,:)];
    elseif h1.BinEdges(i) < iq(4)
        b.CData(i,:) = [colorss(4,:)];

    end
end

set(gcf,'Position',[  171   482   216   154])
%ylim([0 10])


cd /Users/ajorge/Desktop/hyperdirectfigs
filenamee = sprintf('barplot_%2.0f',1);
saveas(gcf,filenamee,'epsc')

xlabell = sprintf('Cortical EP_%1.0f Amplitude (uV)',1);
xlabel(xlabell )
ylabel('Count')

















