Subject = 'DBS3002';
Session = 3;

Val = squeeze(mean(mean(mean_zscore_Cue(251:1750, findfreq(75, 115, fq), :),1),2));

Val = squeeze(mean(mean(mean_zscore_SpOnset(251:1750, findfreq(75, 115, fq), :),1),2));

% Val = squeeze(mean(mean(mean_zscore_SpOnset(251:1750, findfreq(75, 115, fq), :)-mean_zscore_Cue(251:1750, findfreq(75, 115, fq), :),1),2));


ValueRange(1)=min(min(Val));
ValueRange(2)=max(max(Val));
ValueRange = [-max(abs(ValueRange)), max(abs(ValueRange))];



Val = (Val - ValueRange(1))/(ValueRange(2)-ValueRange(1));


[ a, ~ ] = project2verts( CortElecLoc, cortex.vert );
a = cell2mat(a);
figure;
colormap jet;

d = 2.4;
de = 1.2;
tau = 1.2;
cm = colormap;


V = zeros(length(cortex.vert),1);
V_color = 1*ones(length(cortex.vert),3);

for v = 1:length(a)
    aeconn = find(pdist2(cortex.vert(a(v),:), cortex.vert)<=de);
    aconn = find(pdist2(cortex.vert(a(v),:), cortex.vert)>de & ...
        pdist2(cortex.vert(a(v),:), cortex.vert)<=d);
    dconn = pdist2(cortex.vert(a(v),:), cortex.vert(aconn,:));
%     V(aeconn) = V(aeconn)+Val(v);
%     V(aconn) = V(aconn)+Val(v)*exp(-(dconn-de)/tau)';
    V(aeconn) = arrayfun(@(x) max(Val(v), x), V(aeconn));
    V(aconn) = arrayfun(@max, Val(v)*exp(-(dconn-de)/tau)', V(aconn));
end

for v = find(V~=0)'
    V_color(v,:) = getColor(V(v), [0 1], cm);
end


Hp = patch('vertices',cortex.vert,'faces',cortex.tri),'FaceVertexCData', V_color,'edgecolor','none','FaceColor','interp',...
    'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5);
axis equal
camlight('headlight','infinite');
axis off;
colormap jet

elec = reshape(cell2mat(CortElecLoc),3,length(CortElecLoc))';
cm = colormap;
%frange = [min(Val), max(Val)];
for e=1:length(CortElecLoc)
    hold on; plot3(elec(e,1), elec(e,2), elec(e,3), 'o', 'markersize', 6, 'color', getColor( Val(e), [0 1], cm ));
end


set(gca,'CameraPosition',DispCamPos.cp,...
    'CameraTarget',DispCamPos.ct,...
    'CameraViewAngle',DispCamPos.cva,...
    'CameraUpVector',DispCamPos.uv);

DispCamPos.cp = campos;
DispCamPos.cva = camva;
DispCamPos.ct = camtarget;
DispCamPos.uv = camup;

saveas(gcf, [Subject, '_Session', num2str(Session), '_CueGamma.fig'], 'fig');
saveas(gcf, [Subject, '_Session', num2str(Session), '_CueGamma.pdf'], 'pdf');

saveas(gcf, [Subject, '_Session', num2str(Session), '_SpOnsetGamma.fig'], 'fig');
saveas(gcf, [Subject, '_Session', num2str(Session), '_SpOnsetGamma.pdf'], 'pdf');