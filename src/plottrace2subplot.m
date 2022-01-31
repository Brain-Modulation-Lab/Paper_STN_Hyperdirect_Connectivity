function plottrace2subplot(brodmann, x, y, EPt, EPv) 

figure(301)
switch brodmann
    case 6
        subplot(4,2,1)
        title('premotor')
    case 4
        subplot(4,2,3)
         title('m1')
    case 3
        subplot(4,2,5)
         title('s1')
    case 22
        subplot(4,2,7) 
         title('stg')

    case 45
        subplot(4,2,2)
        title('pars t.')
    case 44
        subplot(4,2,4)
         title('pars o.')
    case 43
        subplot(4,2,6)
         title('subcentral')


end

plot(x,y)
hold on
ylim([-7 7])
plot(EPt, EPv,'+r')



