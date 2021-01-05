function [ FRA , FRA_significance_mat ] = make_FRA_from_window2( sorted_RP ,...
           protocol_struct , spont_rate , evoked_window_vec , spontaneous_window_vec )
    
    % The function MAKE_FRA_FROM_WINDOW2 receives a raster plot
    % (sorted_RP), a pure tone protocol file (protocol_struct) , a vector
    % with all time point of the evoked window (evoked_window_vec) and a
    % vector of the spontaneous time window time points(
    % spontaneous_window_vec) and calculates the units FRA and tests for
    % the significance of the auditory responses for each (frequency x
    % attenuation) pairs
    %
    % The variable spont_rate is not used in this version of the function

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
        evoked_firing_rates( mm ) = evoked_firing_rates( mm ) - spontaneous_FR_pertrials ;
        
    end
    
    FRA = reshape( evoked_firing_rates , length( protocol_struct.atten ) , length( protocol_struct.freqs) ) ;
    FRA_significance_mat = reshape( significane_rez , length( protocol_struct.atten ) , length( protocol_struct.freqs) ) ;
    
end

