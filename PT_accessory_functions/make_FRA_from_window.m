function [ FRA , FRA_significance_mat ] = make_FRA_from_window( sorted_RP , protocol_struct , spont_rate , evoked_window_vec , spontaneous_window_vec )
    
    % A function that calculates the FRA while subtracting the same
    % spontaneous firing rate ("spont_rate") from all responses

    MSEC_IN_SEC = 1000 ;
    
    freqs_times_attens = length( protocol_struct.freqs) * length( protocol_struct.atten ) ;
    n_reps = protocol_struct.n_reps ;
    FRA_significance_mat = zeros( length( protocol_struct.atten ) , length( protocol_struct.freqs) ) ;
    
    for mm = 1 : freqs_times_attens
    
        evoked_firing_rates( mm ) = ( mean( sum( sorted_RP( ( mm * n_reps -( n_reps -1 ) ) : mm * n_reps , evoked_window_vec ) , 2 ) )/ length( evoked_window_vec ) ) * MSEC_IN_SEC ; 
        freq_trials = sum( sorted_RP( mm * n_reps -( n_reps -1 ) : mm * n_reps , evoked_window_vec ) , 2 ) ;
        baseline_trials = sum( sorted_RP( mm * n_reps -( n_reps -1 ) : mm * n_reps , spontaneous_window_vec ) , 2 )  ;
        [ ~ , significane_rez( mm ) ] = ranksum( freq_trials ,baseline_trials , 'Tail' , 'right' ) ;
        spontaneous_FR_pertrials = ( mean( sum( sorted_RP( ( mm * n_reps -( n_reps -1 ) ) : mm * n_reps , spontaneous_window_vec ) , 2 ) )/ length( spontaneous_window_vec ) ) * MSEC_IN_SEC ; 
        
%         evoked_firing_rates( mm ) = evoked_firing_rates( mm ) - spontaneous_FR_pertrials ;
        
    end
    
    evoked_firing_rates = evoked_firing_rates - spont_rate ;
    FRA = reshape( evoked_firing_rates , length( protocol_struct.atten ) , length( protocol_struct.freqs) ) ;
    FRA_significance_mat = reshape( significane_rez , length( protocol_struct.atten ) , length( protocol_struct.freqs) ) ;
    
end

