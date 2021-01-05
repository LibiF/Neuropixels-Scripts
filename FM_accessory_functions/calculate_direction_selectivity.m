function [ total_selectivity_index_log ,  total_selectivity_index_lin , ...
           total_selectivity_ind_allstims ] = calculate_direction_selectivity( data_table ,...
           protocol_struct , PRE_STIM_MSEC , WINDOW_SIZE )

    % The function CALCULATE_DIRECTION_SELECTIVITY receives the data table,
    % a FM protocol struct, the pre stimulus baseline time (PRE_STIM_MSEC)
    % and the size of the defined response window (WINDOW_SIZE) and
    % calculates the degree of direction selectivity for all units. The
    % function outputs several metrics of direction selectivity:
    % total_selectivity_index_log - selectivity index for logarithmic FMs
    % total_selectivity_index_lin - selectivity index for linear FMs
    % total_selectivity_ind_allstims - selectivity index for all FMs
    % grouped
            
    MSEC_IN_SEC = 1000 ;
    
    % Extract parameters from FM protocol
    n_reps = protocol_struct.n_reps ;
    n_stims_per_prot = 2 * length( protocol_struct.oct_speeds ) ;
    ind_vec = [ 1 : n_reps : n_reps * n_stims_per_prot ] ;
    stim_times = round( [ protocol_struct.time_of_stim , fliplr( protocol_struct.time_of_stim ) ] .* MSEC_IN_SEC )' ; 

    % Prepare empty data structures
    total_selectivity_index_log = zeros( size( data_table , 1 ) ,1 ) ;
    total_selectivity_index_lin = zeros( size( data_table , 1 ) ,1 ) ;
    total_selectivity_ind_allstims = zeros( size( data_table , 1 ) ,1 ) ;
  
    % Calculate the DSIs
    for kk = 1 : size( data_table , 1 )

        log_resp_mat = data_table.log_responses{ kk , 1 } ;        
        lin_resp_mat = data_table.lin_responses{ kk , 1 } ;
        
        best_lin_win = data_table.best_lin_win{ kk , 1 } ;
        best_log_win = data_table.best_log_win{ kk , 1 } ;
        
        activity_sum_log = zeros( n_stims_per_prot , 1 ) ;
        activity_sum_lin = zeros( n_stims_per_prot , 1 ) ;
        
        for mm = 1 : n_stims_per_prot
    
            log_trials = log_resp_mat( [ ind_vec(mm) : ind_vec(mm) + n_reps - 1 ] , : ) ;
            lin_trials = lin_resp_mat( [ ind_vec(mm) : ind_vec(mm) + n_reps - 1 ] , : ) ;

            activity_sum_log( mm , 1 ) = sum( sum( log_trials( : , ( PRE_STIM_MSEC + best_log_win( mm )  ) : ( PRE_STIM_MSEC + best_log_win( mm ) + WINDOW_SIZE  ) ) ) ) ./WINDOW_SIZE .* MSEC_IN_SEC;
            activity_sum_lin( mm , 1 ) = sum( sum( lin_trials( : , ( PRE_STIM_MSEC + best_lin_win( mm )  ) : ( PRE_STIM_MSEC + best_lin_win( mm ) + WINDOW_SIZE  ) ) ) )  ./WINDOW_SIZE .* MSEC_IN_SEC ;
            
        end    
         
        direct_selectivity_log_total = ( sum( activity_sum_log( end - 0.5* n_stims_per_prot +1 : end  ) ) - sum( activity_sum_log( 1 : 0.5* n_stims_per_prot ) ) ) ./ sum( activity_sum_log ) ; 
        direct_selectivity_lin_total = ( sum( activity_sum_lin( end - 0.5* n_stims_per_prot +1 : end  ) ) - sum( activity_sum_lin( 1 : 0.5* n_stims_per_prot ) ) ) ./ sum( activity_sum_lin ) ; 

        per_stim_selectivity_log = reshape( activity_sum_log , 0.5* n_stims_per_prot  , 2 ) ;
        per_stim_selectivity_log( : , 2 ) = flipud( per_stim_selectivity_log( : , 2 ) ) ;
        per_stim_selectivity_log = ( per_stim_selectivity_log( : , 2 ) - per_stim_selectivity_log( : ,1 ) ) ./ sum( per_stim_selectivity_log , 2 ) ;
        
        per_stim_selectivity_lin = reshape( activity_sum_lin , 0.5* n_stims_per_prot  , 2 ) ;
        per_stim_selectivity_lin( : , 2 ) = flipud( per_stim_selectivity_lin( : , 2 ) ) ;
        per_stim_selectivity_lin = ( per_stim_selectivity_lin( : , 2 ) - per_stim_selectivity_lin( : , 1 ) ) ./ sum( per_stim_selectivity_lin , 2 ) ;
       
        total_selectivity_index_log( kk , 1) = direct_selectivity_log_total  ;
        total_selectivity_index_lin( kk , 1) = direct_selectivity_lin_total ;
    
        all_rising_FMs = [ activity_sum_lin(  end - 0.5* n_stims_per_prot +1 : end ) ; activity_sum_log(  end - 0.5* n_stims_per_prot +1 : end ) ] ; 
        all_falling_FMs = [ activity_sum_lin( 1 : 0.5* n_stims_per_prot ) ; activity_sum_log(   1 : 0.5* n_stims_per_prot ) ] ; 
        
        total_DSI_allstims = ( sum( all_rising_FMs) - sum(all_falling_FMs) ) ./ ( sum( all_rising_FMs) + sum(all_falling_FMs) ) ;
        total_selectivity_ind_allstims( kk , 1 ) = total_DSI_allstims ;
        
    end   
    
end

