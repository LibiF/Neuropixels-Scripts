function [ logFM_region_matrices , linFM_region_matrices ] = create_FM_responsemat_fordprime( data_table ,...
                FM_protocol_struct , brain_regions , WINDOW_SIZE , PRE_STIM_MSEC ) 

    % The function CREATE_FM_RESPONSE_MAT_FORDPRIME receives the data
    % table, the FM protocol struct, a list of brain regions and parameters
    % such as the response window size (WINDOW_SIZE) and the pre stimulus
    % baseline time (PRE_STIM_MSEC). 
    % The function creates response matrices for linear and logarithmic FM
    % responses that will fit the "findDprime" function.
    % The output matrices are 3 dimensional with the following sizes:
    % 1 - size of units per region
    % 2 - number of stimuli per protocol
    % 3 - number of trial repetitions per stimulus
            
    n_stims = length( FM_protocol_struct.oct_speeds ).*2 ;
    n_reps = FM_protocol_struct.n_reps ;
    n_attens = length( FM_protocol_struct.atten ) ;

    FM_trial_vecs = 1 : n_reps : n_stims* n_reps *n_attens ; 
    FM_trial_vecs = repmat( [ 0 : n_reps-1]' , 1 , n_stims ) + FM_trial_vecs ; 
         
    logFM_region_matrices = cell( size( brain_regions , 2 ) , 1 ) ;
    linFM_region_matrices = cell( size( brain_regions , 2  ) , 1 ) ;

    for kk = 1 : size( data_table , 1 )
            
        lin_win_onset = data_table.best_lin_win{kk,1} ;        
        log_win_onset = data_table.best_log_win{kk,1} ; 

        lin_raster = data_table.lin_responses{kk , 1 } ;
        log_raster = data_table.log_responses{kk , 1 } ;

        for nn = 1 : n_stims

            stim_trials = FM_trial_vecs( : , nn ) ;
            lin_stim_trials = sum( lin_raster( stim_trials , PRE_STIM_MSEC + lin_win_onset(nn) + 1 : PRE_STIM_MSEC + lin_win_onset(nn) + WINDOW_SIZE ) , 2 ) ; 
            linFM_response_matrix( kk , nn , : ) = lin_stim_trials ;

            log_stim_trials = sum( log_raster( stim_trials , PRE_STIM_MSEC + log_win_onset(nn) + 1 : PRE_STIM_MSEC + log_win_onset(nn) + WINDOW_SIZE ) , 2 ) ; 
            logFM_response_matrix( kk , nn , : ) = log_stim_trials ;

        end

    end
       
   for mm = 1 : size( brain_regions , 2 )
        
        region_indices = find( data_table.acronym == brain_regions{ 1 , mm } ) ;
        logFM_region_matrices{mm,1} = logFM_response_matrix( region_indices , : , : ) ;
        linFM_region_matrices{mm,1} = linFM_response_matrix( region_indices , : , : ) ;
        
    end
    
end

