function plothistovoltdistribution(y_loc,y_loc_n_sig,y_loc_n_all,colorss,brodmannoi,plotfactor2,c6,EPoi,iq1,med,iq3,method)




%y_loc = y_loc./y_loc_n_all;




%% interquantile method

if method == 4

    close(figure(44))
    figure(44)
    h1 = histogram(y_loc,15);
    figure(45)
    %subplot(1,2,plotfactor2)
    b = bar(h1.BinEdges(1:(end-1)), h1.Values, 0.9);
    xlabell = sprintf('Cortical EP_%1.0f Amplitude (uV)',EPoi);
    xlabel(xlabell )
    ylabel('STN Stim Site Count')
    b.FaceColor = 'flat';

    for i = 1:length(h1.BinEdges(1:(end-1)))
        if h1.BinEdges(i) < iq1
            b.CData(i,:) = [colorss(1,:)];
        elseif h1.BinEdges(i) <= med
            b.CData(i,:) = [colorss(2,:)];
        elseif h1.BinEdges(i) <= iq3
            b.CData(i,:) = [colorss(3,:)];
        elseif h1.BinEdges(i) > iq3
            b.CData(i,:) = [colorss(4,:)];

        end
    end
    set(gcf,'Position',[  171   482   216   154])
    %ylim([0 10])
    close(44)


    cd /Users/ajorge/Desktop/hyperdirectfigs
    filenamee = sprintf('barplot_%2.0f',brodmannoi);
    saveas(gcf,filenamee,'epsc')



elseif method == 3

    close(figure(44))
    figure(44)
    h1 = histogram(y_loc,20);
    figure(45)
    %subplot(1,2,plotfactor2)
    b = bar(h1.BinEdges(1:(end-1)), h1.Values, 0.9);
    xlabell = sprintf('Cortical EP_%1.0f Amplitude (uV)',EPoi);
    xlabel(xlabell )
    ylabel('STN Stim Site Count')
    b.FaceColor = 'flat';

    for i = 1:length(h1.BinEdges(1:(end-1)))
        
        if h1.BinEdges(i) < iq1
            b.CData(i,:) = [colorss(1,:)];
        elseif h1.BinEdges(i) <= med
            b.CData(i,:) = [colorss(2,:)];
        else
            b.CData(i,:) = [colorss(3,:)];
        end
        
    end
    set(gcf,'Position',[  171   482   216   154])
    ylim([0 12])
    close(44)


    cd /Users/ajorge/Desktop/hyperdirectfigs
    filenamee = sprintf('barplot_%2.0f',brodmannoi);
    saveas(gcf,filenamee,'epsc')

end



