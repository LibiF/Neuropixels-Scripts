function [ min_lin_windows_mat , min_log_windows_mat , lin_signif_mat ,...
                log_signif_mat ] = find_min_window_in_FMresponses( data_table ,...
                protocol_struct , window_size , post_stim_window , PRE_STIM_MSEC )
    
    % The function FIND_MIN_WINDOW_IN_FMRESPONSES receives a data table,
    % an FM protocol file, and a desired window size and finds the window
    % with minimal response for each FM presented. 
    % The parameter "post_stim_window" - sets the latest point for the
    % response window to end (stimulus-end + post_stim_window)
    % The parameter "PRE_STIM_MSEC" - sets the baseline time of activity
    % before the FM stimuli were presented.
    % The function outputs:
    % min_lin_windows_mat - a cell size nx1 where n is the number of units
    % in the data_table. In each row, there is a vector with k values of
    % the onsets for the minimal response windows (most suppressive windows) 
    % to all linear FMs (k=number of linear FMs presented. In my protocol 
    % this is 10).
    % lin_signif_mat - A cell containing a boolean vector in each row
    % marking whether the suppressive response window was significantly 
    % different from baseline.
    % 2 analougous variables are outputted for logarithmic FMs
            
    MSEC_IN_SEC = 1000 ;
    
    n_reps = protocol_struct.n_reps ;
    n_stims_per_prot = 2 * length( protocol_struct.oct_speeds ) ;
    ind_vec = [ 1 : n_reps : n_reps * n_stims_per_prot ] ;
    stim_times = round( [ protocol_struct.time_of_stim , fliplr( protocol_struct.time_of_stim ) ] .* MSEC_IN_SEC )' ; 
    
    min_lin_windows_mat = cell( size( data_table , 1 ) , 1 ) ;
    min_log_windows_mat = cell( size( data_table , 1 ) , 1 ) ;
    
    lin_signif_mat = cell( size( data_table , 1 ) , 1 ) ;
    log_signif_mat = cell( size( data_table , 1 ) , 1 ) ;
    
    for kk = 1 : size( data_table , 1 )
        
        log_resp_mat = data_table.log_responses{ kk , 1 } ;        
        lin_resp_mat = data_table.lin_responses{ kk , 1 } ;
        
        min_lin_unit_mat = zeros( 1 , n_stims_per_prot ) ;
        signif_lin_unit_mat = zeros( 1 , n_stims_per_prot ) ;
        min_log_unit_mat = zeros( 1 , n_stims_per_prot ) ;
        signif_log_unit_mat = zeros( 1 , n_stims_per_prot ) ;
        
        for mm = 1 : n_stims_per_prot
            
            lin_trials = lin_resp_mat( [ ind_vec(mm) : ind_vec(mm) + n_reps - 1 ] , : ) ;
            log_trials = log_resp_mat( [ ind_vec(mm) : ind_vec(mm) + n_reps - 1 ] , : ) ;
            
            possible_windows = stim_times( mm ) + post_stim_window - window_size ;
            sum_lin = zeros( possible_windows , 1 ) ;
            sum_log = zeros( possible_windows , 1 ) ;
            
            for ww = 1 : possible_windows
            
                sum_lin( ww , 1 ) = sum( sum( lin_trials( : , PRE_STIM_MSEC + ww : PRE_STIM_MSEC + ww + window_size - 1 ) ) ) ;
                sum_log( ww , 1 ) = sum( sum( log_trials( : , PRE_STIM_MSEC + ww : PRE_STIM_MSEC + ww + window_size - 1 ) ) ) ;
  
            end
            
            [ ~ , min_lin_ind ] =  min( sum_lin ) ;
            [ ~ , min_log_ind ] =  min( sum_log ) ;
                        
            spont_lin_trials = sum( lin_trials( : , [ PRE_STIM_MSEC - window_size + 1 : PRE_STIM_MSEC ] ) , 2 ) ;
            spont_log_trials = sum( log_trials( : , [ PRE_STIM_MSEC - window_size + 1 : PRE_STIM_MSEC ] ) , 2 ) ;

            min_win_lin_trials = sum( lin_trials( : , [ ( PRE_STIM_MSEC + min_lin_ind ) : ( PRE_STIM_MSEC + min_lin_ind + window_size-1 ) ] ) , 2 ) ;
            min_win_log_trials = sum( log_trials( : , [ ( PRE_STIM_MSEC + min_log_ind ) : ( PRE_STIM_MSEC + min_log_ind + window_size -1 ) ] ) , 2 ) ;
            
            [ ~ , h_lin ] = ranksum(  min_win_lin_trials , spont_lin_trials , 'tail' , 'left' ) ;  
            [ ~ , h_log ] = ranksum(  min_win_log_trials , spont_log_trials , 'tail' , 'left' ) ;  

            min_lin_unit_mat( 1 , mm ) = min_lin_ind ;
            signif_lin_unit_mat(1 , mm ) = h_lin ;
            min_log_unit_mat( 1, mm ) = min_log_ind ;
            signif_log_unit_mat( 1 , mm ) = h_log ;            

        end   
        
        lin_signif_mat{ kk , 1 } = signif_lin_unit_mat ;
        log_signif_mat{ kk , 1 } = signif_log_unit_mat ; 
        
        min_lin_windows_mat{ kk , 1 } = min_lin_unit_mat ;            
        min_log_windows_mat{ kk , 1 } = min_log_unit_mat ;
        
    end   

end

