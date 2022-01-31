function distances = STN2lipskiaxis( vector_tail, vector_head, avg_aux, ulocsi, brodmannoi, plotornot)
%
%
% ajorge
%
%








%% code used pre 190927 

% distance of the points on the "STN axis"

%{
%make a line between m_center and a_center
if plotornot == 1
    hold on
    figure(105)
end

% make line to connect between 2 coordinates
startp = stn_m_center(1:2);
endp   = stn_a_center(1:2);
pts    = 1000;
m=(endp(2)-startp(2))/(endp(1)-startp(1)); %gradient 
xx=linspace(startp(1)-10,endp(1)+10,pts);
yy=m*(xx-startp(1))+startp(2);

startp = stn_m_center(2:3);
endp   = stn_a_center(2:3);
m=(endp(2)-startp(2))/(endp(1)-startp(1)); %gradient 
zz=m*(xx-startp(1))+startp(2);

if plotornot == 1
    plot3(xx,yy,zz)
    hold on
end

%find closest point on the STN axis to the Stimuation coordinate
distt = sqrt(sum(([xx' yy' zz'] -   dummy  ).^2,2));
[~, distt_min] = min(distt);

if plotornot == 1
    plot3(xx(distt_min)+randn(1)/100, yy(distt_min)+randn(1)/100, zz(distt_min)+randn(1)/100,'o','MarkerSize',20,'MarkerFaceColor',avg_aux.cortical_colors(brodmannoi,:))
end
hold on

coord_save = [xx(distt_min), yy(distt_min), zz(distt_min)];

distances = sqrt(sum(((stn_m_center-coord_save).^2),2));
%}


















%% code used after 190927 

if plotornot == 1
    figure(107)
    hold on
elseif plotornot == 2
    figure(106)   
    hold on
elseif plotornot == 4
    figure(1001)
end





%% plot STN 

if plotornot == 1 || plotornot == 6
    cd('/Users/ajorge/software2/lead/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL Minimal (Ewert 2017)')
    load('atlas_index.mat') 
    sidee = 2;
    patch('vertices',atlases.fv{1,sidee}.vertices,'faces', atlases.fv{1,sidee}.faces,'facecolor',atlases.colormap(1,:),'FaceAlpha',0.05,'EdgeColor','none'); hold on
    patch('vertices',atlases.fv{2,sidee}.vertices,'faces', atlases.fv{2,sidee}.faces,'facecolor',atlases.colormap(2,:),'FaceAlpha',0.05,'EdgeColor','none')
    patch('vertices',atlases.fv{3,sidee}.vertices,'faces', atlases.fv{3,sidee}.faces,'facecolor',atlases.colormap(3,:),'FaceAlpha',0.05,'EdgeColor','none')


end






%% find 1D line
% find equation of a line in 3D space (ref: https://brilliant.org/wiki/3d-coordinate-geometry-equation-of-a-line/):
% consider a line which passes thru P = (x1 y1 z1) 
% and has direction D = (l,m,n).
% Let X = (x,y,z) be a random point on the line. 
% Then the vector PX is parallel to D, thus:
% PX = td ... t = (x-x1)/l = (y-y1)/m = (z-z1)/n

P = vector_tail;
Q = vector_head;

D = P - Q;

x = -18:0.0001:-6;

if plotornot == 4
    x = -13.5:0.001:-10;
end

y = ((x - P(1))/D(1))*D(2) + P(2);
z = ((y - P(2))/D(2))*D(3) + P(3);

if plotornot == 1 ||  plotornot == 2 || plotornot == 5
    plot3(x,y,z,'.k'); hold on
    plot3(P(1),P(2),P(3),'sr','MarkerSize',10)
    hold on
    plot3(Q(1),Q(2),Q(3),'g^','MarkerSize',10)
    hold on
end

if plotornot == 4
    %plot3(P(1),P(2),P(3),'k.','MarkerSize',10)
    %plot3(Q(1),Q(2),Q(3),'k.','MarkerSize',10)
    %plot3(x,y,z,'.k'); hold on
    QQ = [x(1) y(1) z(1)];
    PP = [x(end) y(end) z(end)];
    arrow3(QQ,PP,'k',3,3)
    %arrow3(vtail',vhead','k',arrowhead1,arrowhead2)
end



%% define 0 in our 1D line as the location of the projection of STN_motor onto this 1D line


stn_motor_center = [-12.6901  -15.0986   -7.1118];
stn_assc_center  = [-10.4126  -11.7141   -7.6281];

%stn_motor_center = [-13.99 -14.74 -5.4];
%stn_assc_center  = [-12.97 -12.57 -5.08];



if plotornot == 1 ||  plotornot == 2
    plot3(stn_motor_center(1), stn_motor_center(2), stn_motor_center(3),'om')
    text(stn_motor_center(1), stn_motor_center(2), stn_motor_center(3),'origin')
end


clear distt
distt = sqrt(sum(      ([x' y' z'] - stn_motor_center).^2        ,2));
[~, line1d_zero] = min(distt);

if plotornot == 1 ||  plotornot == 2 || plotornot == 9
    plot3(x(line1d_zero),y(line1d_zero),z(line1d_zero),'om','MarkerSize',15)
    plot3(x(line1d_zero),y(line1d_zero),z(line1d_zero),'om','MarkerSize',10)
    plot3(x(line1d_zero),y(line1d_zero),z(line1d_zero),'om','MarkerSize',5)
    text( x(line1d_zero),y(line1d_zero),z(line1d_zero),'origin')

    hold on
    axis equal
end






%% project points onto 1D line

% original point
if plotornot == 1 ||  plotornot == 2 || plotornot > 5
    plot3(ulocsi(1),ulocsi(2),ulocsi(3),'ks','MarkerSize',10,'markerfacecolor',avg_aux.cortical_colors(brodmannoi,:))
    hold on
end

%find point on 1D line that is the closest to ulocsi
clear distt
distt = sqrt(sum(      ([x' y' z'] - ulocsi).^2        ,2));
[~, inx] = min(distt);

%final point
if plotornot == 1 ||  plotornot == 2 || plotornot > 5
    plot3(x(inx),y(inx),z(inx),'ko','MarkerSize',7,'markerfacecolor',avg_aux.cortical_colors(brodmannoi,:))
end

%connecting dots
x2 = [ulocsi(1) x(inx) ];
y2 = [ulocsi(2) y(inx) ];
z2 = [ulocsi(3) z(inx) ];

plot3(x2, y2, z2,'k.-')
hold on

%distances = (inx - line1d_zero).*0.01;
%distances = sqrt(sum(([x(inx),y(inx),z(inx)] - stn_motor_center).^2,2));
distances = sqrt(sum(([x(inx),y(inx),z(inx)] - [x(line1d_zero),y(line1d_zero),z(line1d_zero)]).^2,2));




% doing this before 2021 02 06
%determine if behind or forward of origin
%{
if inx >= line1d_zero
elseif inx < line1d_zero
    distances = -distances; 
end
%}


% doing this after 2021 02 06
if inx < line1d_zero
elseif inx >= line1d_zero
    distances = -distances; 
end




