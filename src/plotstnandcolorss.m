function [outt, vector_tail, vector_head] = plotstnandcolorss(plotfactor,brodmannoi,y_loc,y_loc_n_sig,y_loc_n_all,colorss,EPoi,iq1,med,iq3,ulocs,avg_aux, partt, method)









%% 
if partt == 1
    figure(plotfactor+brodmannoi+partt)
elseif partt == 2
    figure(433)
    hold on
elseif partt == 3
    figure(500+brodmannoi)
elseif partt == 4
    figure(1000)
end






%% plot stn
%subplot(2,1,1)
cd('/Users/ajorge/software2/lead/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL Minimal (Ewert 2017)')
load('atlas_index.mat') 
sidee = 2;
patch('vertices',atlases.fv{1,sidee}.vertices,'faces', atlases.fv{1,sidee}.faces,'facecolor',atlases.colormap(1,:),'FaceAlpha',0.05,'EdgeColor','none'); hold on
patch('vertices',atlases.fv{2,sidee}.vertices,'faces', atlases.fv{2,sidee}.faces,'facecolor',atlases.colormap(2,:),'FaceAlpha',0.05,'EdgeColor','none')
patch('vertices',atlases.fv{3,sidee}.vertices,'faces', atlases.fv{3,sidee}.faces,'facecolor',atlases.colormap(3,:),'FaceAlpha',0.05,'EdgeColor','none')







%% plot STN with little bubbles

%subplot(5,1,[1,2,3,4])
for i = 1:length(y_loc)
    
    noplot=0;
    
    
    %all the same size
    sizee = 25;
    
    %size depends on number of ecog contacts
    %{
    if y_loc_n_all(i)<=10;  sizee =  5; end
    if y_loc_n_all(i)<=20; sizee = 10; end
    if y_loc_n_all(i)<=30; sizee = 15; end
    if y_loc_n_all(i)> 40;  sizee = 25; end
    %}
    
    
    if method == 4
        %pick color QUARTIILE METHOD
        if y_loc(i) <=iq1
            colorsss(i,:) = colorss(1,:);
            classs(i) = 1;
        elseif y_loc(i) > iq1 && y_loc(i) <= med
            colorsss(i,:) = colorss(2,:);
            classs(i) = 2;
        elseif y_loc(i) > med && y_loc(i) <=iq3
            colorsss(i,:) = colorss(3,:);
            classs(i) = 3;
        elseif y_loc(i)>iq3
            colorsss(i,:) = colorss(4,:);
            classs(i) = 4;
        elseif isnan(y_loc(i))
            %plot3(ulocs(i,1), ulocs(i,2), ulocs(i,3), 'o','MarkerSize',sizee2,'MarkerEdgeColor','k')
            hold on
            noplot = 1;
        else
            colorsss(i) = [0 0 0];
            noplot = 1;
        end
        %}

    elseif method == 3
        %pick color TERTILE METHOD
        if y_loc(i) <=iq1
            colorsss(i,:) = colorss(1,:);
            classs(i) = 1;
        elseif y_loc(i) > iq1 && y_loc(i) <= med
            colorsss(i,:) = colorss(2,:);
             classs(i) = 2;
        elseif y_loc(i) > med 
            colorsss(i,:) = colorss(3,:);
             classs(i) = 3;
        elseif isnan(y_loc(i))
            hold on
            noplot = 1;
        else
            colorsss(i) = [0 0 0];
            noplot = 1;
        end
        %}
        
    end
    
    
    
    

    
    
   
    if noplot==0
        %subplot(2,1,1)
        if partt == 1
            plot3(ulocs(i,1), ulocs(i,2), ulocs(i,3), 'o','MarkerSize', sizee,'MarkerFaceColor',colorsss(i,:),'MarkerEdgeColor',colorsss(i,:))
            hold on
        end
        
        if partt == 3 && classs(i) == 4
            plot3(ulocs(i,1), ulocs(i,2), ulocs(i,3), 'o','MarkerSize', sizee,'MarkerFaceColor',colorsss(i,:),'MarkerEdgeColor',colorsss(i,:))
            hold on
        else
            plot3(ulocs(i,1), ulocs(i,2), ulocs(i,3), 'o','MarkerSize', sizee,'MarkerEdgeColor',[0.5 0.5 0.5])
        end
        
    end
    
    %plot all available channels on cortical area of interest
    %plot3(ulocs(i,1), ulocs(i,2), ulocs(i,3), 'o','MarkerSize',y_loc_n_all(i),'MarkerEdgeColor','k')
    
end

%grid











%% find axis with LDA
clc
mdl = fitcdiscr(ulocs,classs');

n_bigballs = size(mdl.ClassNames,1);

%plot 4 balls (red, orange, light blue, blue)
if partt == 1 || partt == 2
    for i = 1:n_bigballs
        %plot3(mdl.Mu(i,1), mdl.Mu(i,2), mdl.Mu(i,3), 'o','MarkerSize', 50,'MarkerFaceColor',colorss(i,:),'MarkerEdgeColor',colorss(i,:)) 
        hold on
    end
end

if partt == 3
    plot3(mdl.Mu(n_bigballs,1), mdl.Mu(n_bigballs,2), mdl.Mu(n_bigballs,3),'o','MarkerSize',50,'MarkerFaceColor',avg_aux.cortical_colors(brodmannoi,:),'LineWidth',1,'MarkerEdgeColor',avg_aux.cortical_colors(brodmannoi,:)) 
    %{
    plot3(mdl.Mu(n_bigballs,1), mdl.Mu(n_bigballs,2), mdl.Mu(n_bigballs,3),'.','MarkerSize',100,'Color',avg_aux.cortical_colors(brodmannoi,:)) 
    plot3(mdl.Mu(n_bigballs,1), mdl.Mu(n_bigballs,2), mdl.Mu(n_bigballs,3),'o','MarkerSize',50,'Color',avg_aux.cortical_colors(brodmannoi,:),'LineWidth',3) 
    plot3(mdl.Mu(n_bigballs,1), mdl.Mu(n_bigballs,2), mdl.Mu(n_bigballs,3),'o','MarkerSize',75,'Color',avg_aux.cortical_colors(brodmannoi,:),'LineWidth',3) 
    %}
end

outt.bigredbubble(brodmannoi,:) = [mdl.Mu(n_bigballs,1), mdl.Mu(n_bigballs,2), mdl.Mu(n_bigballs,3)];

% dynamic vector (changes according to LDA for specific cortical area)
vector_head = [mdl.Mu(n_bigballs,1),mdl.Mu(n_bigballs,2),mdl.Mu(n_bigballs,3)];
vector_tail = [mdl.Mu(1,1),mdl.Mu(1,2),mdl.Mu(1,3)];

% static vector (M1) (one vector to reflect them all)
%vector_tail = [-13.5485  -14.5225   -5.8487];
%vector_head = [ -12.3479  -13.8637   -6.1459];


if partt == 2
    plot3(vector_tail(1), vector_tail(2), vector_tail(3), 'r*', 'MarkerSize',8)
    plot3(vector_head(1), vector_head(2), vector_head(3), 'g*', 'MarkerSize',8)

    arrowhead1 = 1.5; %2.5;
    arrowhead2 = 1.5; %4.0;
    arrow3(vector_tail, vector_head,'k',arrowhead1,arrowhead2)
end

%save these values
outt.stn_m_center = vector_tail;
outt.stn_a_center = vector_head;










%% find axis with center of STN motor and STN assoc DEPRECRATED

%{
stn_m_center = mean(atlases.fv{2,sidee}.vertices);
plot3(stn_m_center(1), stn_m_center(2), stn_m_center(3), 'r*', 'MarkerSize',8)
stn_a_center = mean(atlases.fv{3,sidee}.vertices);
plot3(stn_a_center(1), stn_a_center(2), stn_a_center(3), 'g*', 'MarkerSize',8)


%redefine stm centers optimized for the x and y planes (STN LONG AXIS)
%stn_m_center = [-13.03, -15.12, -4.849];
%stn_a_center = [-9.912, -11.61, -6.009];

% STN axis from motor center EP to premotor cnete EP
%stn_m_center = [-13.03, -15.12, -4.849];
%stn_a_center = [-13.7, -13.33, -3.776];

% STN premotor optimized vector
%stn_m_center = [-13.4 -17.4 -6.2];
%stn_a_center = [-14.3 -10.6 -3.5];

% STN motor optimized vector
%stn_m_center = [-14.3 -17 -6.6];
%stn_a_center = [-11.5 -10.6 -3.5];

%stn_m_center = ([-13.4 -17.4 -6.2] + [-14.3 -17 -6.6])./2;
%stn_a_center = ([-14.3 -10.6 -3.5] + [-11.5 -10.6 -3.5])./2;

if partt == 2
    plot3(vector_tail(1), vector_tail(2), vector_tail(3), 'r*', 'MarkerSize',8)
    plot3(vector_head(1), vector_head(2), vector_head(3), 'g*', 'MarkerSize',8)

    arrowhead1 = 1; %2.5;
    arrowhead2 = 1; %4.0;
    arrow3(vector_tail, vector_head,'k',arrowhead1,arrowhead2)
end

%save these values
outt.stn_m_center = vector_tail;
outt.stn_a_center = vector_head;
%}


















%% plot center of voltage 

    

% (leave 20% out)
%{
for i = 1:10
    a = 1;
    b = length(y_loc);
    r = round((b-a).*rand(round(length(y_loc)*0.80),1) + a);
    centerofvolt = nansum(y_loc(r)'.*ulocs(r,:))./nansum(y_loc(r)');
    %plot3(centerofvolt(1), centerofvolt(2), centerofvolt(3),'.k','MarkerSize',22,'Color',avg_aux.cortical_colors(brodmannoi,:))
    centerofvolt_save(:,:,:,i,brodmannoi) = centerofvolt;
end
%}



%plot center of voltage (concentric circles)
centerofvolt = nansum(y_loc'.*ulocs)./nansum(y_loc);
outt.centerofvolt = centerofvolt;

if partt == 2
    plot3(centerofvolt(1), centerofvolt(2), centerofvolt(3),'.','MarkerSize',100,'Color',avg_aux.cortical_colors(brodmannoi,:)) 
elseif partt == 1  
    %{
    plot3(centerofvolt(1), centerofvolt(2), centerofvolt(3),'.','MarkerSize',100,'Color',avg_aux.cortical_colors(brodmannoi,:)) 
    plot3(centerofvolt(1), centerofvolt(2), centerofvolt(3),'o','MarkerSize',50,'Color',avg_aux.cortical_colors(brodmannoi,:),'LineWidth',3) 
    plot3(centerofvolt(1), centerofvolt(2), centerofvolt(3),'o','MarkerSize',75,'Color',avg_aux.cortical_colors(brodmannoi,:),'LineWidth',3) 
    %}
end



sig = 2*nanstd((y_loc'.*ulocs)./nansum(y_loc'));
r = mean(sig);
th = 0:pi/50:2*pi;
xunit = sig(1) * cos(th) + centerofvolt(1);
yunit = sig(2) * sin(th) + centerofvolt(2);
zunit = sig(3) * cos(th) + centerofvolt(3);
%h = plot3(xunit, yunit, zunit,'LineWidth',5,'Color',avg_aux.cortical_colors(brodmannoi,:));
%ellipsoid(centerofvolt(1),centerofvolt(2),centerofvolt(3),sig(1),sig(2),sig(3))

%plot center of locations only
%centerofvolt = nanmean(ulocs);
%plot3(centerofvolt(1), centerofvolt(2), centerofvolt(3),'o','MarkerSize',20,'Color',avg_aux.cortical_colors(brodmannoi,:) )

%sig = nanstd((y_loc(r)'.*ulocs(r,:))./nansum(y_loc(r)'));
%lot3(centerofvolt(1), centerofvolt(2), centerofvolt(3),'.','MarkerSize',40,'Color',avg_aux.cortical_colors(brodmannoi,:) )





view([0 90])   %axial view
%view([-90 0]) %coronal view


set(gcf,'Position',[  690           1        1193        1200])
xlabel('x mni (mm)','FontSize',40)
ylabel('y mni (mm)','FontSize',40)
a = get(gca,'XTickLabel');
set(gca,'fontsize',40)

% limitis to highlight all STN
xlim([-16.5 -6])
ylim([-20   -8])

%lim to highlight motor STN
%xlim([-16.1 -10])
%ylim([-19 -12])

box off
%set(gca,'visible','off')


cd /Users/ajorge/Desktop/hyperdirectfigs
filenamee = sprintf('STN_%2.0f',partt*100+brodmannoi);
saveas(gcf,filenamee,'png')











