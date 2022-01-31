%% Script plot ECOG EP magnitude/latency on the STN
%
% plots EPs on the cortex
%
% options to normalize EPs on the cortex
%
% ahmed jorge march 28 2021
%



%% load avg
clear all
close all
clc

cd('/Users/ajorge/ajorge/data/ecog_stnpathways/avg/')

%load('avg_190817')

load('avg_191125')
load('avg_auxiliary')






%% 
clc
clear index1 centerofvolt_save centerofvolt_save2
close all

c6=1;
plotfactor2=1;

bigredbubble = nan(50,3);


%
cd /Users/ajorge/Desktop/hyperdirectfigs



val_mat_list = []
subject_list = []
xe_indxed = []
x52=[]
x53=[]



%for ptoi = [    3012        3016        3017        3018        3019        3023 3024        3026        3027        3028        3029        3030        3032]
%

for brodmannoi = [44 45 6 4 3 43 22 40]



    %% build array
    EPoi = 2;
    plotfactor = 100;

    index1 = find(avg.brodmann==brodmannoi & avg.stim_amp==3)% & avg.sett_was_chosen==0 );

    c=1; clear yd xd xe
    %go thru each trace of interest
    for i = 1:length(index1)

        traceoi = index1(i);

        %if there is a significant peak for EPoi (the EP#, e.g. EP2=the second peak associated with ~4-8ms)
        if isempty(avg.EPv{traceoi})==0
            if isempty(avg.EPt{traceoi}(EPoi)) == 0 && isnan(avg.EPt{traceoi}(EPoi)) == 0
                yd(c,1)   = avg.EPv{traceoi}(EPoi); 
                td(c,1)   = avg.EPt{traceoi}(EPoi);
                ydata(c,:) = avg.data{traceoi};
            else
                yd(c,1)   = nan;
                td(c,1)   = nan;
                ydata(c,:) = nan*ones(1,270);

            end
            xd(c,:)   = avg.STN_mni{traceoi};
            xe(c,:)   = avg.ECOG_mni(traceoi,:);
            xen(c,:)  = avg.electrode_n(traceoi);
            pt(c,:)   = avg.pt(traceoi);
            
            c=c+1;
        end

    end
    
    %get rid of nan yd
    inx = find(isnan(yd));
    yd(inx) = [];
    td(inx) = [];
    ydata(inx,:) = [];
    xd(inx,:) = [];
    xe(inx,:) = [];
    xen(inx,:) = [];
    pt(inx,:) = [];

    %get rid of voltage outliers
    inx = find( yd>10)% | yd<1);
    yd(inx) = [];
    td(inx) = [];
    ydata(inx,:) = [];
    xd(inx,:) = [];
    xe(inx,:) = [];
    xen(inx,:) = [];
    pt(inx,:) = [];
    
    %get rid of latencies <2 and >4.5
    %{
    inx = find(td<0.002 | td>0.004);
    yd(inx) = [];
    td(inx) = [];
    ydata(inx,:) = [];
    xd(inx,:) = [];
    xe(inx,:) = [];
    xen(inx,:) = [];
    pt(inx,:) = [];
    %}
    
    %{
    % simple histogram
    figure(31)
    histogram(yd,50)
    figure(30)
    boxplot(yd)
    
    %plot
    figure(32)
    plot(avg.time{1}, ydata')
    hold on
    plot(td,yd,'rs')
    %}
    
    
    
     %% normalize yloc (implemeted on March 28 2021)
    ydd = yd./max(yd);





    %% voltage on STN (organized)
    %--------------------------------------------------------------------------
    % find out unique STN stim locations
    clear ulocs y_loc y_loc_n y_loc_n_all y_loc_n_sig
    ulocs = unique(xd(:,1:3),'rows');
    if length(unique(ulocs(:,1))) ~= length(ulocs)
        sprintf('deal with this')
        pause
    end
    %plot unique locs
    %{
    figure(55)
    hold on
    plot3(ulocs(:,1), ulocs(:,2), ulocs(:,3), 'r+','MarkerSize',8)
    %}

    % --------------------------------------------------------------------------
    % find all voltage values associated with each unique location, and average into one 
    for i = 1:length(ulocs)
        inx = find(xd(:,1) == ulocs(i,1));
        
        y_loc(i)       = nanmean(yd(inx));
        y_loc_n_all(i) = length(inx);   %e.g. all ECoGs on motor cortex 
        y_loc_n_sig(i) = length(inx) - sum(isnan(yd(inx))); %e.g. only ECoGs that are sig in motor cortex 
        
        %y_loc(i) = y_loc(i)*y_loc_n_sig(i);
    end
    
    %get rid of nans
    inx             = find(isnan(y_loc));
    y_loc(inx)      = [];
    y_loc_n_all(inx)= [];
    y_loc_n_sig(inx)= [];
    ulocs(inx,:)    = [];
    
    
    
    %% normalize yloc
    %y_loc = y_loc./max(y_loc);
   
    figure(31)
    subplot(2,1,1)
    histogram(y_loc_n_all,30)
    subplot(2,1,2)
    histogram(y_loc_n_sig,30)

    %% --------------------------------------------------------------------------
    % define colorss
    
    method = 4; %4=interquantile, 3=intertertile
    
    %interquantiles based on each subset of the data (dynamic)   
    if method == 4
        colorss = flipud([215,25,28;253,174,97;171,217,233;44,123,182])./256;
        med = nanmedian(y_loc);
        iq1 = quantile(y_loc,0.25);
        iq3 = quantile(y_loc,0.75);
        %}
        
        %interquantiles based on all traces (static)
        %{
        iq1 = avg_aux.iq1_v_all;
        med = avg_aux.iq2_v_all;
        iq3 = avg_aux.iq3_v_all;
        %}
        
    elseif method == 3
        colorss = flipud([215,25,28;253,174,97;44,123,182])./256;
        iq1 = quantile(y_loc,0.33);
        med = quantile(y_loc,0.67);
        iq3 = nan;

    end
    
    
    
    
    
    %% plot histogram of voltage distirbutions (really a bar plot cuz .eps and illustratror are annoying)
    plothistovoltdistribution(y_loc, y_loc_n_sig, y_loc_n_all, colorss, brodmannoi, plotfactor2, c6, EPoi, iq1, med, iq3, method)



    
    %%
    [outt, vector_tail, vector_head] = plotstnandcolorss(plotfactor,brodmannoi,y_loc,y_loc_n_sig,y_loc_n_all,colorss,EPoi,iq1,med,iq3,ulocs,avg_aux,1,method);   
    bigredbubble(brodmannoi,:) = outt.bigredbubble(brodmannoi,:);
    
    [outt, vector_tail, vector_head] = plotstnandcolorss(plotfactor,brodmannoi,y_loc,y_loc_n_sig,y_loc_n_all,colorss,EPoi,iq1,med,iq3,ulocs,avg_aux,2,method);

    [outt, vector_tail, vector_head] = plotstnandcolorss(plotfactor,brodmannoi,y_loc,y_loc_n_sig,y_loc_n_all,colorss,EPoi,iq1,med,iq3,ulocs,avg_aux,3,method);

    
    
    
    
    %% get center of cortical STN for this region
    corticalcenter(brodmannoi,:) = mean(xe,1);
        

    %% --------------------------------------------------------------------------
    % voltage vs Radius
    %voltage_vs_radiusfrommotorSTN(ulocs,plotfactor,avg_aux,c6,brodmannoi)


    
    %% Project invidivual STN stim locations onto STN long axis
    
    clear distancesi c xx yy
    c=1;
    for i = 1:length(y_loc)
        if isnan(y_loc(i))~=1
            figure(81)
            xx(c) = STN2lipskiaxis( vector_tail, vector_head, avg_aux, ulocs(i,:), brodmannoi, 0);
            yy(c) = y_loc(i);
            c=c+1;
            
            %match for brodmann=6
            cd('/Users/ajorge/ajorge/src/stn_pathways/' )

            %if brodmann==6 && 
            
            
        end
    end
    
    %xx=-xx
    %remove outliers
    %inx = find(yy>5);
        %xx(inx) = [];
    %yy(inx) = [];
    
    
    % plot individual points
    figure(660+brodmannoi)
    %plot(xx,yy,'.','MarkerSize',40,'Color',avg_aux.cortical_colors(brodmannoi,:))
    plot(xx,yy,'.k','MarkerSize',20)
    hold on

    %fit model
    mdl = fitglm((xx),yy);
    hold on
    plot(xx,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2).*xx,'k','LineWidth',3)
    [rho,pval_spearman] = corr(xx',yy','Type','Spearman');
    xlabel('Distance from STN_m_o_t_o_r (mm)')
    ylabel('Cortical Evoked Potential (\muV)')
    textt = sprintf('rho=%3.2f p=%4.3f',rho,pval_spearman);
    text(1,8,textt);
    cd /Users/ajorge/Desktop/hyperdirectfigs
    filenamee = sprintf('relatioship_%2.0f.eps',brodmannoi);
    saveas(gcf,filenamee,'epsc')

    
    
    
    
    
    %% make Chen plot (cortex)
    
    %get all 1dimensional STN locations
    
    %vector_head  = [-12.6901  -15.0986   -7.1118];  %STN_assc center
    %vector_tail  = [-10.4126  -11.7141   -7.6281];  %STN_moto center
    
    clear x5 inx
    for i = 1:length(xd)
        x5(i) = STN2lipskiaxis( vector_tail, vector_head, avg_aux, xd(i,:), brodmannoi, 0);
    end
    
    %get all "anterior STN locations"
    %{
    if brodmannoi==4 || brodmannoi==3 || brodmannoi==22 || brodmannoi==43 || brodmannoi==44
        inx = find(x5<=-2);
    elseif brodmannoi==6 || brodmannoi==40 || brodmannoi==45
        inx = find(x5>1);
    end
    
    %}
    
    %get all "posterior STN locations"
    %{
    if brodmannoi==4 || brodmannoi==3 || brodmannoi==22 || brodmannoi==43 || brodmannoi==44
        inx = find(x5>0);
    elseif brodmannoi==6 || brodmannoi==40 || brodmannoi==45
        inx = find(x5<1);
    end
    %}
    
    
    %get all locations
    
    inx = 1:length(ydd);
    %}
    
   
    x52 = [x52 x5];
    x53 = [x53 x5(inx)];
    val_mat_list = [val_mat_list; ydd(inx)];
    subject_list = [subject_list; pt(inx)];
    xe_indxed = [xe_indxed; xe(inx,:)];
    
    
    
    %plot on brain
    %plot_EP2cortex(subject_list, val_mat_list, xe_indxed)
    %}


end





% make Chen plot (cortex)

%plot on brain
plot_EP2cortex_v2(subject_list, val_mat_list, xe_indxed)
set(gcf,'pos',[ 1285         483         982         718])
%}


































%% center of red-bubble-voltage on STN



%plot the red bubbles 
for brodmannoi = [44 45 6 4 43 3 22 40]
    figure(1001)
    plot3(bigredbubble(brodmannoi,1), bigredbubble(brodmannoi,2), bigredbubble(brodmannoi,3),'ok','MarkerSize',50,'MarkerFaceColor', avg_aux.cortical_colors(brodmannoi,:))
    hold on
end

%plot STN
cd('/Users/ajorge/software2/lead/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL Minimal (Ewert 2017)')
load('atlas_index.mat') 
sidee = 2;
patch('vertices',atlases.fv{1,sidee}.vertices,'faces', atlases.fv{1,sidee}.faces,'facecolor',atlases.colormap(1,:),'FaceAlpha',0.05,'EdgeColor','none'); hold on
patch('vertices',atlases.fv{2,sidee}.vertices,'faces', atlases.fv{2,sidee}.faces,'facecolor',atlases.colormap(2,:),'FaceAlpha',0.05,'EdgeColor','none')
patch('vertices',atlases.fv{3,sidee}.vertices,'faces', atlases.fv{3,sidee}.faces,'facecolor',atlases.colormap(3,:),'FaceAlpha',0.05,'EdgeColor','none')
view([0 90])   %axial view

set(gcf,'Position',[21   781   560   420])

%for formal plot
formalplot=1;
set(gcf,'Position',[  690           1        1193        1200])
xlabel('x mni (mm)','FontSize',40)
ylabel('y mni (mm)','FontSize',40)
a = get(gca,'XTickLabel');
set(gca,'fontsize',40)

%legend
figure(1003)
c=1;
for brodmannoi = [6 40 45 44 3 22 43 4]
    plot(0,c,'ok','MarkerSize',8,'MarkerFaceColor', avg_aux.cortical_colors(brodmannoi,:)); hold on
    textt = sprintf('%s',avg_aux.cortical_names{brodmannoi});
    text(0.2,c,textt)
    set(gca,'fontsize',2)
    c=c-1;
end
axis off
set(gcf,'Position',[    1426         778         213          82])
cd /Users/ajorge/Desktop/hyperdirectfigs
filenamee = sprintf('legend.eps');
saveas(gcf,filenamee,'epsc')






figure(1002)
stn_a_center  = [-12.6901  -15.0986   -7.1118];
stn_m_center  = [-10.4126  -11.7141   -7.6281];

plot3(stn_a_center(1), stn_a_center(2), stn_a_center(3), '*r','MarkerSize',30); hold on
plot3(stn_m_center(1), stn_m_center(2), stn_m_center(3), '*g','MarkerSize',30)

clear distances xx yy
c=1;
for brodmannoi = [44 45 6 4 43 3 22 40]
    
    figure(1003)
    distances(brodmannoi) = STN2lipskiaxis(stn_m_center, stn_a_center, avg_aux, bigredbubble(brodmannoi,:), brodmannoi,1);
    set(gcf,'Position',[560   781   560   420])
    xx(c) = (corticalcenter(brodmannoi,2));
    xx2(c)= abs(corticalcenter(brodmannoi,2));
    yy(c) = distances(brodmannoi);
    
    %for formal plot, comment the following line out, run block, then
    %uncomment line
    dummy(brodmannoi) = STN2lipskiaxis(stn_m_center, stn_a_center, avg_aux, bigredbubble(brodmannoi,:), brodmannoi,4);
    
    
    figure(533)
    plot((xx(c)),    yy(c),   'ok','MarkerSize',20,'MarkerFaceColor', avg_aux.cortical_colors(brodmannoi,:)) 
    set(gcf,'Position',[  1112         781         560         420])
    hold on
    
    figure(534)
    plot((xx2(c)),    yy(c),   'ok','MarkerSize',20,'MarkerFaceColor', avg_aux.cortical_colors(brodmannoi,:)) 
    set(gcf,'Position',[ 1654         781         560         420])
    hold on
    
    
    c=c+1;
end




%% relationship between cortex and STN

figure(533)

% absolute value relationship
figure(534)
mdl = fitglm((xx2),yy);
hold on
plot(xx2,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2).*xx2,'k','LineWidth',3);
textt2 = sprintf('p-val=%4.3f', mdl.Coefficients.pValue(2));
text(0,2,textt2);
[rho,pval_spearman] = corr(xx2',yy','Type','Spearman');
xlabel('Absolute Distance From Central Sulcus (mm)');
ylabel('STN Center of Voltage (mm)');
set(gca, 'XDir','reverse')
[rho,pval_pearson] = corr(xx2',yy','Type','Pearson');
%textt = sprintf('linear model p=%4.3f, Pearson p=%4.3f, Spearman p=%4.3f \n', mdl.Coefficients.pValue(2),pval_pearson, pval_spearman)
%text(20,0.5,textt)


% quadratic model
% I used the Matlab curve fitting tool, polynomial order 2, weights = aoi = [ 6     4     3    22    43    44    45    40]
figure(533)
p1 =     0.001383 ;%(6.546e-05, 0.0004828)
p2 =   -0.001527 ;%(0.002119, 0.01024)
p3 =     -0.3973 ;% (0.3098, 0.5716)
 
x_hat = -35:0.1:32;
y_hat = p1.*x_hat.^2 + p2.*x_hat + p3;
plot(-x_hat,y_hat,'k','LineWidth',3)

xlabel('Distance From Central Sulcus (mm)');
ylabel('STN Center of Voltage (mm)');
set(gca, 'XDir','reverse')

%}






figure(1001)
cd /Users/ajorge/Desktop/hyperdirectfigs
filenamee = sprintf('allcenters.eps');
saveas(gcf,filenamee,'epsc')

figure(533)
cd /Users/ajorge/Desktop/hyperdirectfigs
filenamee = sprintf('quadratic.eps');
saveas(gcf,filenamee,'epsc')
set(gcf,'Position',[ 1395         286         336         285])

figure(534)
cd /Users/ajorge/Desktop/hyperdirectfigs
filenamee = sprintf('linear.eps');
saveas(gcf,filenamee,'epsc')
set(gca, 'XDir','normal')
set(gcf,'Position',[  1735         290        336         285])























%% plot center of ECOG stim in CORTEX


figure(301)
%load MNI brain
cd('/Users/ajorge/ajorge/data/ecog_stnpathways')
load('cortex_MNI.mat')
%plot MNI brain
atlas=2; % specify Desikan-Killiany atlas
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
hold on

centerr = nanmean(xe)
plot3(centerr(1), centerr(2), centerr(3),'*','MarkerSize',40 )







































































