%SCRIPT peaks to cortex v2






%% load data

clear all
close all
clc

cd('/Users/ajorge/ajorge/data/ecog_stnpathways/avg/')

load('avg_190609')

load('avg_auxiliary')





%---------------------------------------------------------------------------
%% get voltage values of interest
clear yd subject_list val_mat_list


ptoi = load_patient_list_STNhyperdirect(2);
yd_count=[];

ptoix =  3016%[        3012        3016        3017        3018        3019        3023        3024        3027        3028        3029        3030        3032  ]  ;
%ptoix2 = [                                              3018                   3023              3026  3027        3028        3029        3030        3032 ];



c2=1;
for ptoi = ptoix
    
    clear yd
    
    %figure out if patient had one or two ECoG arrays
    inx = find(avg.pt==ptoi);
    maxx = max(avg.electrode_n(inx));
    if maxx<65
        yd = nan*ones(1,64);
        maxx = 63;
    elseif ptoi == 3012
        yd = nan*ones(1,63+4);
        maxx = 67;
    else
        yd = nan*ones(1,128);
        maxx = 126;
    end
    
    
    %what ecogs are we going to plot
    index1 = find(avg.pt==ptoi & avg.stim_amp==3  & avg.sett_was_chosen==1 & (avg.brodmann==4 | avg.brodmann==6 | avg.brodmann==3 | avg.brodmann==22 | avg.brodmann==43  | avg.brodmann==45 ));
    
    if length(index1)>0
        EPoi = 1;
        
        %go thru each trace
        for i = 1:length(index1)

            traceoi = index1(i);

            %if there is a significant peak for EPoi (the EP#, e.g. EP2=the second peak associated with ~4-8ms)
            if isempty(avg.EPv{traceoi})==0
                c3 = 1;
                c=avg.electrode_n(traceoi);
                
                if isempty(avg.EPt{traceoi}(EPoi)) == 0 && isnan(avg.EPt{traceoi}(EPoi)) == 0
                    
                    
                    
                    c=avg.electrode_n(traceoi);
                    yd(1,c)   = avg.EPv{traceoi}(EPoi);
                    c3=c3+1

                else
                    
                    
                    
                    yd(1,c)   = nan;
                    c3=c3+1
                    
                end
            end

        end


        inx = find(yd<0.1 | yd>10);
        yd(inx) = nan;
        yd(inx) = nan;
        
        %inx = find(yd<0.0045 | yd>0.0065);
        %yd(inx) = nan;
        %yd(inx) = nan;
        

        %figure
        %histogram(yd)
        
        if maxx==126
        	yd([64 128]) = [];
        elseif ptoi == 3012
            %yd([68 69]) = [];
        elseif maxx==63
            yd([64:end]) = [];
        end

        %build input
        subject_list{c2} = sprintf('DBS%4.0f',ptoi);
        val_mat_list{c2} = yd(1,1:maxx);
        c2=c2+1;
        
        yd_count = [yd_count yd(1,1:maxx)];
        
    else
    end


end

figure(303)
histogram(yd_count,40)
xlabel('Cortical Evoked Potential (\muV)')
ylabel('ECoG Electrode Count')

    
    
%

%

%% ---------------------------------------------------------------------------
% get tt (vertices) 
%close all
clc


    
tt = val2cortex(subject_list, val_mat_list, 'median');


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


%---------------------------------------------------------------------------
% plot
figure
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








%%


index1 = find(avg.pt==3018 & avg.sett_was_chosen==1 & avg.stim_amp==3 & avg.brodmann==6 )

avg.electrode_n(index1)-63
















