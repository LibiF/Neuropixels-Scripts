function [ FRA , FRA_significance_mat ] = make_FRA_from_suppressive_window( sorted_RP ,...
           protocol_struct , evoked_window_vec , spontaneous_window_vec )

    % The function MAKE_FRA_FROM_SUPPRESSIVE_WINDOW receives a raster plot
    % (sorted_RP), a pure tone protocol file (protocol_struct) and times of
    % both spontaneous and evoked response windows and calculates an FRA
    % for suppressed units. It also calculates the significance of
    % responses in terms of suppression (meaning using a left-tailed test)

    MSEC_IN_SEC = 1000 ;
    
    freqs_times_attens = length( protocol_struct.freqs) * length( protocol_struct.atten ) ;
    n_reps = protocol_struct.n_reps ;
    FRA_significance_mat = zeros( length( protocol_struct.atten ) , length( protocol_struct.freqs) ) ;
    
    for mm = 1 : freqs_times_attens
    
        evoked_firing_rates( mm ) = ( mean( sum( sorted_RP( ( mm * n_reps -( n_reps -1 ) ) : mm * n_reps , evoked_window_vec ) , 2 ) )/ length( evoked_window_vec ) ) * MSEC_IN_SEC ; 
        freq_trials = sum( sorted_RP( mm * n_reps -( n_reps -1 ) : mm * n_reps , evoked_window_vec ) , 2 ) ;
        baseline_trials = sum( sorted_RP( mm * n_reps -( n_reps -1 ) : mm * n_reps , spontaneous_window_vec ) , 2 )  ;
        [ ~ , significane_rez( mm ) ] = ranksum( freq_trials ,baseline_trials , 'Tail' , 'left' ) ;
        spontaneous_FR_pertrials = ( mean( sum( sorted_RP( ( mm * n_reps -( n_reps -1 ) ) : mm * n_reps , spontaneous_window_vec ) , 2 ) )/ length( spontaneous_window_vec ) ) * MSEC_IN_SEC ; 
        evoked_firing_rates( mm ) = evoked_firing_rates( mm ) - spontaneous_FR_pertrials ;
        
    end
    
    FRA = reshape( evoked_firing_rates , length( protocol_struct.atten ) , length( protocol_struct.freqs) ) ;
    FRA_significance_mat = reshape( significane_rez , length( protocol_struct.atten ) , length( protocol_struct.freqs) ) ;
    
end

