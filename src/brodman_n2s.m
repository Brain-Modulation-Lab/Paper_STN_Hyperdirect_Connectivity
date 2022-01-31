function titlee = brodman_n2s(brodmann_n)



switch brodmann_n
    
    case 22
        titlee = 'STG';
    case 4
        titlee = 'M1';
    case 6
        titlee = 'Premotor';
    case 3 
        titlee = 'S1';
    case 43 
        titlee = 'Subcentral';
    case 44 
        titlee = 'Pars O.';
    case 45 
        titlee = 'Pars T.';
    case 40 
        titlee = 'Supramarginal';
        
    otherwise
        titlee = '';
        
end
        
        
        



