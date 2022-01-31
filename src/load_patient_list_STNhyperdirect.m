function ptoi = load_patient_list_STNhyperdirect(option)

%% latest patient list for STN hyperdirect pathway project


% option = 1; gives all patients of interest (better to obtain all preprocessing files that we can)
if option == 1
    ptoi = [3004 3005 3006 3008 3010 3011 3012 3016 3017 3018 3019 3020 3022 3023 3024 3025 3026 3027 3028 3029 3030 3031 3032 4068 4079 4081 4085 4088];
end




% option = 2; gives final patients (not all patients will meet all requirements to process (e.g. pts with noisy ECOG recordings)
if option == 2
    ptoi = [ 3005  3008 3010 3011 3012 3016 3017 3018 3019    3023 3024  3026 3027 3028 3029 3030 3032          4068 4079 4081 4085 4088];
end



