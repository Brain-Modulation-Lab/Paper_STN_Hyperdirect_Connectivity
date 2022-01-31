%% SCRIPT make brodmann tables



%load patients
ptoi = load_patient_list_STNhyperdirect(1);

%initialize the lookup table
ecog_brodmann_lookupt = nan*ones(3100,126);


list_failed_pts=0;
for pt = ptoi
        
    try
        %get location
        [array1, array2] = ecog_cortexlocs_auto(pt,0);

        ecog_brodmann_lookupt(pt,[1:63 64 65:127]) = [array1(1:63)' nan array2(1:63)'];
    catch
        list_failed_pts = [list_failed_pts pt];
       
    end
    
end


list_failed_pts

cd('/Users/ajorge/ajorge/src/stn_pathways/')
save('ecog_brodmann_lookupt','ecog_brodmann_lookupt')




