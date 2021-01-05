function region_matrices = create_PT_responsemat_fordprime( data_table ,...
         PT_protocol_struct , brain_regions , WINDOW_SIZE , PRE_STIM_MSEC ) 

    % The function CREATE_PT_RESPONSEMAT_FORDPRIME receives a data table,
    % the pure tone protocol file, a cell of brain regions, the size of the
    % units response window and the time of baseline activity pre stimulus
    % presentation.
    % The function outputs a cell "region_matrices" containing trial 
    % matrices for each brain region.
    % Each matrix is a 3D matrix with dimensions in the following sizes:
    % N x S x rep
    % N - num of units in a regions
    % S - num of stimuli
    % rep - number of trials/repetitions of each stimulus
    
    % The function essentially takes data from the units raster plot and 
    % rearranges it to fit the input format of the function "findDprime" 


    n_stims = length( PT_protocol_struct.freqs ) ;
    n_reps = PT_protocol_struct.n_reps ;
    n_attens = length( PT_protocol_struct.atten ) ;

    trial_vecs = n_reps + 1 : n_reps* n_attens : n_stims* n_reps *n_attens ; 
    trial_vecs = repmat( [0 : n_reps-1]' , 1 , n_stims ) + trial_vecs ; 

    response_matrix = zeros( size( data_table , 1 ) , n_stims , n_reps ) ;
    region_matrices = cell( size( brain_regions , 2 ) , 1 ) ;

    for kk = 1 : size( data_table , 1 )

        if data_table.is_excited(kk,1) || isnan( data_table.is_inhibited( kk , 1 ) ) || (data_table.is_excited(kk,1)==0 && data_table.is_inhibited(kk,1)==0 ) 

            win_onset = data_table.max_resp_wind( kk ,  1 ) ;        

        else

            win_onset = data_table.min_resp_wind( kk , 1 ) ;

        end   

        raster = data_table.PT_raster{kk , 1 } ;
        for nn = 1 : n_stims

            stim_trials = trial_vecs( : , nn ) ;
            stim_trials = sum( raster( stim_trials , PRE_STIM_MSEC + win_onset + 1 : PRE_STIM_MSEC + win_onset + WINDOW_SIZE ) , 2 ) ; 
            response_matrix( kk , nn , : ) = stim_trials ;

        end

    end
    
    for mm = 1 : size( brain_regions , 2 )
        
        region_indices = find( data_table.acronym == brain_regions{ 1 , mm } ) ;
        region_matrices{mm,1} = response_matrix( region_indices , : , : ) ;
        
    end
    
end

