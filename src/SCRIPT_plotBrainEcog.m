% SCRIPT plot brain and ECoG coordinates
clear all
close all
clc


%% plot brain and ECoG coordinates (MNI space)
clc

close all

for ptoi = 3016%[ 3005 3008 3010 3011 3012 3016 3017 3018 3019 3023 3024 3026 3027 3028 3029 3030  3032]
    
    %load coordinates
    try
        cfg=[];
        cfg.ptoi = ptoi;
        values = ones(1,126);
        cfg.plotbrain=1;
        cfg.enumerate=1; %number the electrodes
        plot_ecog2cortex_v2(cfg,values)
        
        if ptoi == 30333
            cfg.plotbrain=1;
            plot_ecog2cortex_v2(cfg,values)
        end
        
        
     
        
    catch
        ptoi    
    end
    

end












%% plot brain and ECoG coordinates (native space)
clc

for pt = 3027%[ 3006 3008 3011  3016 3017 3018 3019 3020 3022 3023 3024 3025 3026 3027 3028 3029 3030 3031 3032]
    
    %load coordinates
    

    try
        cfg=[];
        cfg.pt = pt;
        cfg.hot = 1;
        cfg.pos1 = [536   228   927   717];
        cfg.mni = 1;  
        cfg.plotcortex=1;
        
        if pt == 3031
            cfg.last=1;
        else
            cfg.last=0;
        end
        
        
        values = ones(1,126);
        plot_ecog2cortex(cfg,values)
        
    catch
        pt    
    end
    

end












