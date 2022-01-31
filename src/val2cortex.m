function vert_val = val2cortex(subject_list, val_mat_list, average_method)

%% inputs:
% subject_list         = {'DBSXXXX', 'DBSXXXX', ...};
% val_mat_list         = {1XN, 1XM, 1XZ,...};
% average_method       = 'mean' or 'median' (default mean);
%% output:
% vert_val         = M X 1 matrix, M is the number of brain vertices

%%%%% can be changed: dd, df, decay function, and average method

if ~exist('average_method')
    average_method = 'mean';
end

%tepuy characteristics (dd>df)
dd   =  8.0;                % diameter of influence (from de to d decay)   (units of distance)
df   =  2.0;               % diameter of full weight                      (units of distance)    

% if onely ECOG electrode footpath wanted, dd-df=0 and df=0.5

% get all the cortelecloc for all the subject in subject_list
clearvars elec_dir
CortElecLoc_total = [];
for subject_id = 1:length(subject_list)
    subject_name = subject_list{subject_id};
    

    % direct to CortElecLoc file
    
    %original
    %{
    elec_dir{subject_id} = dir(['/Volumes/Nexus/DBS/' subject_name '/Anatomy/Freesurfer/*/Electrode_Locations/MNI_on_cortex/CortElecLocL.mat']);
    if isempty(elec_dir{subject_id})
        elec_dir{subject_id} = dir(['/Volumes/Nexus/DBS/' subject_name '/Anatomy/Freesurfer/Electrode_Locations/MNI_on_cortex/CortElecLocL.mat']);
        if isempty(elec_dir{subject_id})
            error(['No CortElecLoc file found for ' subject_name]);
        end
    end
    %}
    %local folder option:
    elec_dir{subject_id} = dir(['/Users/ajorge/ajorge/data/ecog_stnpathways/freesurfer_preop_ElecLoc/' subject_name '_preop_ElecLoc/MNI_on_cortex/CortElecLocL.mat']);
    if isempty(elec_dir{subject_id})
        %elec_dir{subject_id} = dir(['/Volumes/Nexus/DBS/' subject_name '/Anatomy/Freesurfer/Electrode_Locations/MNI_on_cortex/CortElecLocL.mat']);
        if isempty(elec_dir{subject_id})
            error(['No CortElecLoc file found for ' subject_name]);
            pause
        end
    end
    

    if length(elec_dir{subject_id}) ~= 1
        error(['More than one files found for ' subject_name]);
    end
    
    % load in CortElecLoc for this subject
    load([elec_dir{subject_id}.folder filesep elec_dir{subject_id}.name], 'CortElecLoc');
    
    length(CortElecLoc);
    
    % stack them in CortElecLoc_total
    CortElecLoc_total = [CortElecLoc_total, CortElecLoc];
end

% get the total value array also
try
    val_mat_total = cell2mat(val_mat_list);
catch
    val_mat_total = (val_mat_list);
end

% load in MNI brain BS1
load('/Users/ajorge/ajorge/data/ecog_stnpathways/cortex_MNI.mat', 'BS1');


% Next step is for each contact to find its corresponding vertices
[~, vert_idx] = min(pdist2(cell2mat(CortElecLoc_total'), BS1.Vertices),[],2);





% project them onto the cortex
Cort_elec_coords = BS1.Vertices(vert_idx,:);
%%%%%%%%%%%%%% to be added%%%%%%%
assert(length(Cort_elec_coords) == length(val_mat_total), 'Number of contacts not match with number of values');
%remove nans from the list
idx2rm = isnan(val_mat_total);
val_mat_total(idx2rm) = [];
Cort_elec_coords(idx2rm,:) = [];
%%%%%%%%%%%%%%%% 






V = nan(length(BS1.Vertices),1);

for which_vert = 1:length(V)
    %which_vert;
    
    dist_list = pdist2(BS1.Vertices(which_vert,:), Cort_elec_coords);
    
    if min(dist_list) > dd
        V(which_vert) = NaN;
    else
        
        % find contacts that enclose the vertex (DF)
        full_elec_idx = find(dist_list <= df);
        full_weight_list = ones(1,length(full_elec_idx));
        full_val_list = val_mat_total(full_elec_idx);
        
        % find contracts that are enclosed by the decay (DD) distance
        decay_elec_idx = find(dist_list > df ...
            & dist_list <= dd);
        decay_dist_list = dist_list(decay_elec_idx);
        decay_val_list = val_mat_total(decay_elec_idx);
        % assign weigths for decay contacts 
        decay_weight_list = 1 - ((decay_dist_list - 0.5) * (2/3)); % linear decay from 1 to 0
        
        % mean/median DF vertices with DD vertices 
        if strcmp(average_method, 'mean')
            V(which_vert) = ([full_val_list, decay_val_list] * [full_weight_list, decay_weight_list]') / sum([full_weight_list, decay_weight_list]);
        elseif strcmp(average_method, 'median')
            V(which_vert) = weightedMedian([full_val_list, decay_val_list],[full_weight_list, decay_weight_list]);
        end
    end
end

vert_val = V;

        
        
            
        
        
    
    
 




%%%%%%%%%%%%% weighted median function
    function wMed = weightedMedian(D,W)

    % ----------------------------------------------------------------------
    % Function for calculating the weighted median 
    % Sven Haase
    %
    % For n numbers x_1,...,x_n with positive weights w_1,...,w_n, 
    % (sum of all weights equal to one) the weighted median is defined as
    % the element x_k, such that:
    %           --                        --
    %           )   w_i  <= 1/2   and     )   w_i <= 1/2
    %           --                        --
    %        x_i < x_k                 x_i > x_k
    %
    %
    % Input:    D ... matrix of observed values
    %           W ... matrix of weights, W = ( w_ij )
    % Output:   wMed ... weighted median                   
    % ----------------------------------------------------------------------


    if nargin ~= 2
        error('weightedMedian:wrongNumberOfArguments', ...
          'Wrong number of arguments.');
    end

    if size(D) ~= size(W)
        error('weightedMedian:wrongMatrixDimension', ...
          'The dimensions of the input-matrices must match.');
    end

    % normalize the weights, such that: sum ( w_ij ) = 1
    % (sum of all weights equal to one)

    WSum = sum(W(:));
    W = W / WSum;

    % (line by line) transformation of the input-matrices to line-vectors
    d = reshape(D',1,[]);   
    w = reshape(W',1,[]);  

    % sort the vectors
    A = [d' w'];
    ASort = sortrows(A,1);

    dSort = ASort(:,1)';
    wSort = ASort(:,2)';

    sumVec = [];    % vector for cumulative sums of the weights
    for i = 1:length(wSort)
        sumVec(i) = sum(wSort(1:i));
    end

    wMed = [];      
    j = 0;         

    while isempty(wMed)
        j = j + 1;
        if sumVec(j) >= 0.5
            wMed = dSort(j);    % value of the weighted median
        end
    end


    % final test to exclude errors in calculation
    if ( sum(wSort(1:j-1)) > 0.5 ) & ( sum(wSort(j+1:length(wSort))) > 0.5 )
         error('weightedMedian:unknownError', ...
          'The weighted median could not be calculated.');
    end
    end
%%%%%%%%%%%%%%%
end



















