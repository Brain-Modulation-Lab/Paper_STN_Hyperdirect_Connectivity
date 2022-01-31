function [tract_was_chosen, stim_amp, stim_freq] = chosen_tract_during_sx(pt,queried_tract)

% look up table to know which tract was the chosen one (according to my
% word document that highlights which tract was localized (crammond's notes) and accoriding to
% the research stim set number assigned by me) i.e. the tract is numbered in 
% chronological order when the research stim was done

% ajorge dec 2018



%% look up table for chosen tract 

chosen_tract = -5;

switch pt
    
    case 3005
        chosen_tract = [3 4];
        
    case 3006
        chosen_tract = 1;
    
    case 3008
        chosen_tract = 2;
             
    case 3010
        chosen_tract = [2];
         
    case 3011
        chosen_tract = [1 4]; 
        
    case 3012
        chosen_tract = 1;
        
    case 3016
         chosen_tract = [3 5];
        
    case 3017
        chosen_tract = 2;
        
    case 3018
        chosen_tract = 4;
        
    case 3019
        chosen_tract = 2;
    
    case 3020 
        chosen_tract = 3;

    case 3022   
        chosen_tract = 1:3;
        
    case 3023
        chosen_tract = [1 2 3];
        
    case 3024
        chosen_tract = [1 4];
        
    case 3025
        chosen_tract = [1];
          
    case 3026
        chosen_tract = [1];
        
    case 3027
        chosen_tract = [1];
        
    case 3028
        chosen_tract = [1 2];
    
    case 3029
        chosen_tract = [1 2];
    
    case 3030
        chosen_tract = [3 4];
    
    case 3031
        chosen_tract = [2];
        
    case 3032
        chosen_tract = 1;
        
        
    otherwise
        warning('not a valid patient on the lookup table')
        chosen_tract = 999;
        
        
end
         

if any(queried_tract == chosen_tract) %this allows for multiple sets to be designated as "seletected tract" during surgery
    tract_was_chosen = 1;
else
    tract_was_chosen = 0;
end

if chosen_tract == -5
    tract_was_chosen = NaN;
end




%% look up table for stim amplitude and frequency

settt = queried_tract;

if pt == 3005
    
    if settt == 1 || settt == 3 || settt == 5
        stim_amp = 1;
        stim_freq = 1;
    elseif settt == 2 || settt == 4 || settt == 6
        stim_amp = 2;
        stim_freq = 1;     
    end
    
elseif pt == 3006
    stim_amp = 1;
    stim_freq = 1;
     
elseif pt == 3008
    stim_amp = 1;
    stim_freq = 1;
        
elseif pt == 3010
    stim_amp = 1;
    stim_freq = 1;
    
elseif pt == 3011
    stim_amp = 1;
    stim_freq = 1;
    
elseif pt == 3012
    stim_amp = 3;
    stim_freq = 1;
    
elseif pt == 3016
    stim_amp = 3;
    stim_freq = 1;
    
elseif pt == 3017
    stim_amp = 3;
    stim_freq = 1;
    
elseif pt == 3018
    stim_amp = 3;
    stim_freq = 1;
    
elseif pt == 3019
    stim_amp = 3;
    stim_freq = 1;
    
elseif pt == 3020
    stim_amp = 3;
    stim_freq = 1;
    
elseif pt == 3022
    if settt == 1 || settt == 4 || settt == 7
        stim_amp = 1;
    elseif settt == 2 || settt == 5 || settt == 8
        stim_amp = 2;
    elseif settt == 3 || settt == 6 || settt == 9
        stim_amp = 3;
    end
    
elseif pt == 3023
    stim_freq = 1;
    
    if settt == 1 || settt == 4 || settt == 7
        stim_amp = 1;
    elseif settt == 2 || settt == 5 || settt == 8
        stim_amp = 2;
    elseif settt == 3 || settt == 6 || settt == 9
        stim_amp = 3;
    end
   
    
elseif pt == 3024
    stim_amp = 3;
    stim_freq = 1;
    
elseif pt == 3025
    stim_amp = 3;
    stim_freq = 1;
    
    
elseif pt == 3026
    stim_amp = 3;
    stim_freq = 1;
    
elseif pt == 3027
    stim_amp = 3;
    stim_freq = 1;
    
elseif pt == 3028
    stim_amp = 3;
    
    if settt == 1 || settt == 3 || settt == 5 || settt == 7 || settt == 9 || settt == 11
        stim_freq = 1;
    elseif settt == 2 || settt == 4 || settt == 6 || settt == 8 || settt == 10 || settt == 12
        stim_freq = 10;
    end
    
elseif pt == 3029
    stim_amp = 3;
    
    if settt == 1 || settt == 3 || settt == 6 || settt == 7 || settt == 9 || settt == 11
        stim_freq = 1;
    elseif settt == 2 || settt == 4 || settt == 5 || settt == 8 || settt == 10 || settt == 12
        stim_freq = 10;
    end
    
    
elseif pt == 3030
     stim_amp = 3;
    
    if settt == 1 || settt == 3 || settt == 5 || settt == 7 || settt == 9 || settt == 11
        stim_freq = 1;
    elseif settt == 2 || settt == 4 || settt == 6 || settt == 8 || settt == 10 || settt == 12
        stim_freq = 10;
    end
    
    
elseif pt == 3031
    stim_amp = 3;
    stim_freq = 10;
    
elseif pt == 3032    
    stim_amp = 3;
    stim_freq = 10;

    
else
    
    stim_amp = nan;
    stim_freq = nan;
    
end



























