function [ S_stat_lin , S_stat_log , S_stat_total ] = calculate_FM_lifetimesparsness( data_table ,...
           protocol_struct , PRE_STIM_MSEC , WINDOW_SIZE )

    % The function CALCULATE_FM_LIFETIMESPARSENESS receives a data table,
    % an FM protocol struct, a desired window size and time of baseline
    % activity per trial (PRE_STIM_MSEC), and calculates the lifetime 
    % sparseness of FM responses.
    % The calculation is done over linear and logarithmic FMs separately
    % and also for all FM stimuli altogether.
    
    MSEC_IN_SEC = 1000 ;

    n_reps = protocol_struct.n_reps ;
    n_stims_per_prot = 2 * length( protocol_struct.oct_speeds ) ;
    total_stims = 2 * n_stims_per_prot ;
    ind_vec = [ 1 : n_reps : n_reps * n_stims_per_prot ] ;

    S_stat_lin = zeros( size( data_table , 1 ) , 1 ) ;
    S_stat_log = zeros( size( data_table , 1 ) , 1 ) ;
    S_stat_total = zeros( size( data_table , 1 ) , 1 ) ;

for kk = 1 : size( data_table , 1 )
        
        log_resp_mat = data_table.log_responses{ kk , 1 } ;        
        lin_resp_mat = data_table.lin_responses{ kk , 1 } ;
        
        lin_win_times = data_table.best_lin_win{ kk ,1 } ;
        log_win_times = data_table.best_log_win{ kk , 1 } ;
        
        lin_FRs_in_wins = zeros( n_stims_per_prot , 1 ) ;
        log_FRs_in_wins = zeros( n_stims_per_prot , 1 ) ;

        for mm = 1 : n_stims_per_prot
            
            lin_trials = lin_resp_mat( [ ind_vec(mm) : ind_vec(mm) + n_reps - 1 ] , : ) ;
            log_trials = log_resp_mat( [ ind_vec(mm) : ind_vec(mm) + n_reps - 1 ] , : ) ;
            
            lin_FRs_in_wins( mm , 1 ) = sum( sum( lin_trials( : , PRE_STIM_MSEC + lin_win_times( mm ) : PRE_STIM_MSEC + lin_win_times( mm ) + WINDOW_SIZE ) ) ) ./ WINDOW_SIZE .* MSEC_IN_SEC;
            log_FRs_in_wins( mm , 1 ) = sum( sum( log_trials( : , PRE_STIM_MSEC + log_win_times( mm ) : PRE_STIM_MSEC + log_win_times( mm ) + WINDOW_SIZE ) ) ) ./ WINDOW_SIZE .* MSEC_IN_SEC ;

        end

        S_stat_lin(kk,1) = ( 1 - ( ( sum( lin_FRs_in_wins ) / n_stims_per_prot )^2 / ( sum( lin_FRs_in_wins.^2 ) / n_stims_per_prot) ) )/( 1 - ( 1 / n_stims_per_prot ) )  ;
        S_stat_log(kk,1) = ( 1 - ( ( sum( log_FRs_in_wins ) / n_stims_per_prot )^2 / ( sum( log_FRs_in_wins.^2 ) / n_stims_per_prot) ) )/( 1 - ( 1 / n_stims_per_prot ) )  ;
        
        all_FRs = [ lin_FRs_in_wins ; log_FRs_in_wins ] ;
        S_stat_total( kk , 1 ) =  ( 1 - ( ( sum( all_FRs ) / total_stims )^2 / ( sum( all_FRs.^2 ) / total_stims) ) )/( 1 - ( 1 /total_stims ) )  ;
        
end

