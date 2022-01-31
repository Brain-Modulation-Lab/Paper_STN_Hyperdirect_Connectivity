function plot_rnd_individual_avg(avg,avg_aux)


clc

%pick a random subsample of significant peaks to plot
howmany = 100;
totall = 1:length(avg.data);
a = 1;
b = length(avg.data);
r = round((b-a).*rand(howmany,1) + a);

epx=2; %assuming that the second peak is the "3ms peak"
figure(340+epx)


for gg = 1:howmany %r%1:length(avg.data)
    
    trace = r(gg);
    
    if isempty(avg.EPv{trace}) == 0
        if isnan(avg.EPv{trace}(epx))==0 

            switch avg.brodmann(trace)
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
                case 40
                    subplot(4,2,8)
                     title('supramarginal')
            end

            %only plot if min/max voltages between 2ms-5ms are >0.1 and <10
            i8 = find(avg.time{trace}>0.002 & avg.time{trace}<0.006);
            
            if  max(abs(avg.data{trace}(i8)))<10     
                figure(251)
                plot(avg.time{trace}, (avg.data{trace}))
                hold on
                ylim([-7 7])
                plot(avg.EPt{trace}(epx), avg.EPv{trace}(epx),'+r')
                xlim([0 0.010])
            end
            
            
        end
    end
    
end