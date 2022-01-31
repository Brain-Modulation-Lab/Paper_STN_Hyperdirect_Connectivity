function voltage_vs_radiusfrommotorSTN(ulocs,plotfactor,avg_aux,c6,brodmannoi)

stn_motor_center = [-13.45 -15.06 -7.10];

clear normm
for i = 1:length(y_loc)
    normm(i) = norm(ulocs(i,:) - stn_motor_center);
end



figure(1000+plotfactor)
subplot(3,7,c6+7)
c6=c6+1;
plot(normm,y_loc,'o','MarkerSize',4,'MarkerFaceColor',avg_aux.cortical_colors(brodmannoi,:) ,'MarkerEdgeColor',avg_aux.cortical_colors(brodmannoi,:))
xlabel('Distance to Center STN_m_o_t_o_r (mm)')
ylabel('Cortical Evoked Potential (\muV)')


mdl = fitglm(normm,y_loc);
hold on
xt     = [min(normm) max(normm)];
yt     = [mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*min(normm)      mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*max(normm)  ];
plot( [xt(1) xt(2)]   , [yt(1) yt(2)] , '-','LineWidth',3,'Color',avg_aux.cortical_colors(brodmannoi,:))

strr = sprintf('pval=%4.3f',mdl.Coefficients.pValue(2));
text(4,5,strr)

%ylim([0 5])
xlim([0 8])

set(gcf,'Position',[     1685         253         371         478])


cd('/Users/ajorge/Desktop/')
filenamee = sprintf('fit_%2.0f',brodmannoi)
saveas(gcf,filenamee,'epsc')


























