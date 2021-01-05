function [ best_lin_windows_mat , best_log_windows_mat , lin_signif_mat ,log_signif_mat ,...
           peak_resp_lin_mat , peak_resp_log_mat ] = find_best_window_in_FMresponses( data_table ,...
           protocol_struct , window_size , post_stim_window , PRE_STIM_MSEC )

    % The function FIND_BEST_WINDOW_IN_FMRESPONSES receives a data table,
    % an FM protocol file, and a desired window size and finds the window
    % with maximal response for each FM presented. 
    % The parameter "post_stim_window" - sets the latest point for the
    % response window to end (stimulus-end + post_stim_window)
    % The parameter "PRE_STIM_MSEC" - sets the baseline time of activity
    % before the FM stimuli were presented.
    % The function outputs:
    % best_lin_windows_mat - a cell size nx1 where n is the number of units
    % in the data_table. In each row, there is a vector with k values of
    % the onsets for the best excitatory response windows to all linear FMs
    % (k=number of linear FMs presented. In my protocol this is 10).
    % lin_signif_mat - A cell containing a boolean vector in each row
    % marking whether the excitatory response window was statistically
    % significantly different from baseline.
    % peak_response_lin_mat - the time of maximal response during the
    % optimal response window.
    % 3 analougous variables are outputted for logarithmic FMs
       
    MSEC_IN_SEC = 1000 ;
    
    n_reps = protocol_struct.n_reps ;
    n_stims_per_prot = 2 * length( protocol_struct.oct_speeds ) ;
    ind_vec = [ 1 : n_reps : n_reps * n_stims_per_prot ] ;
    stim_times = round( [ protocol_struct.time_of_stim , fliplr( protocol_struct.time_of_stim ) ] .* MSEC_IN_SEC )' ; 
    
    best_lin_windows_mat = cell( size( data_table , 1 ) , 1 ) ;
    best_log_windows_mat = cell( size( data_table , 1 ) , 1 ) ;
    
    lin_signif_mat = cell( size( data_table , 1 ) , 1 ) ;
    log_signif_mat = cell( size( data_table , 1 ) , 1 ) ;
    
    peak_resp_lin_mat = cell( size( data_table , 1 ) , 1 ) ;
    peak_resp_log_mat = cell( size( data_table , 1 ) , 1 ) ;
    
    for kk = 1 : size( data_table , 1 )
        
        log_resp_mat = data_table.log_responses{ kk , 1 } ;        
        lin_resp_mat = data_table.lin_responses{ kk , 1 } ;
        
        best_lin_unit_mat = zeros( 1 , n_stims_per_prot ) ;
        signif_lin_unit_mat = zeros( 1 , n_stims_per_prot ) ;
        best_log_unit_mat = zeros( 1 , n_stims_per_prot ) ;
        signif_log_unit_mat = zeros( 1 , n_stims_per_prot ) ;
        lat_topeak_lin_unit_mat = zeros( 1 , n_stims_per_prot ) ;
        lat_topeak_log_unit_mat = zeros( 1 , n_stims_per_prot ) ;
               
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
            
            [ ~ , max_lin_ind ] =  max( sum_lin ) ;
            [ ~ , max_log_ind ] =  max( sum_log ) ;
            
            spont_lin_trials = sum( lin_trials( : , [ PRE_STIM_MSEC - window_size + 1 : PRE_STIM_MSEC ] ) , 2 ) ;
            spont_log_trials = sum( log_trials( : , [ PRE_STIM_MSEC - window_size + 1 : PRE_STIM_MSEC ] ) , 2 ) ;

            max_win_lin_trials = sum( lin_trials( : , [ ( PRE_STIM_MSEC + max_lin_ind ) : ( PRE_STIM_MSEC + max_lin_ind + window_size-1 ) ] ) , 2 ) ;
            max_win_log_trials = sum( log_trials( : , [ ( PRE_STIM_MSEC + max_log_ind ) : ( PRE_STIM_MSEC + max_log_ind + window_size -1 ) ] ) , 2 ) ;
            
            [ ~ , lat_in_maxwin_lin_trials ] = max( sum( lin_trials( : , [ ( PRE_STIM_MSEC + max_lin_ind ) : ( PRE_STIM_MSEC + max_lin_ind + window_size-1 ) ] ) , 1) ) ;
            [ ~ , lat_in_maxwin_log_trials ]= max( sum( log_trials( : , [ ( PRE_STIM_MSEC + max_log_ind ) : ( PRE_STIM_MSEC + max_log_ind + window_size -1 ) ] ) , 1 ) ) ;
            lat_in_maxwin_lin_trials =  lat_in_maxwin_lin_trials + max_lin_ind ;
            lat_in_maxwin_log_trials = lat_in_maxwin_log_trials + max_log_ind ;
            
            [ ~ , h_lin ] = ranksum(  max_win_lin_trials , spont_lin_trials , 'tail' , 'right' ) ;  
            [ ~ , h_log ] = ranksum(  max_win_log_trials , spont_log_trials , 'tail' , 'right' ) ;  

            best_lin_unit_mat( 1 , mm ) = max_lin_ind ;
            signif_lin_unit_mat(1 , mm ) = h_lin ;
            best_log_unit_mat( 1, mm ) = max_log_ind ;
            signif_log_unit_mat( 1 , mm ) = h_log ;    
            lat_topeak_lin_unit_mat( 1 , mm ) =  lat_in_maxwin_lin_trials ;
            lat_topeak_log_unit_mat( 1 , mm ) = lat_in_maxwin_log_trials ;

        end   
        
        lin_signif_mat{ kk , 1 } = signif_lin_unit_mat ;
        log_signif_mat{ kk , 1 } = signif_log_unit_mat ; 
        
        best_lin_windows_mat{ kk , 1 } = best_lin_unit_mat ;            
        best_log_windows_mat{ kk , 1 } = best_log_unit_mat ;
        
        peak_resp_lin_mat{ kk , 1 } = lat_topeak_lin_unit_mat ;
        peak_resp_log_mat{ kk , 1 } = lat_topeak_log_unit_mat ;
        
    end   

end

