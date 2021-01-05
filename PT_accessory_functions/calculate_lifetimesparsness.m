function S_stat = calculate_lifetimesparsness( data_table , PRE_STIM_MSEC ,...
         window_size ,  protocol_struct , check_excitation )

    % The function CALCULATE_LIFETIMESPARSENESS calculates lifetime
    % sparseness for pure tone responses in excited units.
    % Responses are calcaulted for each unit based on its optimal response
    % window with no spontaneous firing rate subtraction (to avoid cases of
    % "negative firing rates").

    % Extract relevant parameters
    n_reps = protocol_struct.n_reps ;
    n_attens = length( protocol_struct.atten ) ;
    
    % Create empty data structure
    S_stat = zeros( size( data_table , 1 ) , 1 ) ;
    
    for kk = 1 : size( data_table , 1 ) 
        
        raster = data_table.PT_raster{ kk , 1 } ; 
        raster = raster( : , PRE_STIM_MSEC : end ) ;
        
        if check_excitation
            
            is_excitatory = data_table.is_excited( kk , 1 ) ;
        
        else
            
            is_excitatory = true ;
            
        end
        
        if is_excitatory
            
            % Set spontaneous FR to zero so that FRA output will be actual
            % firing rates in the optimal window
            FAKE_SPONT = 0 ;                                                                                                   
            evoked_win_start = data_table.max_resp_wind( kk , 1 ) ;
            evoked_win_vec = [ evoked_win_start : evoked_win_start + window_size - 1 ] ;
            spont_win_vec = [ PRE_STIM_MSEC - window_size + 1 : PRE_STIM_MSEC ] ;

            % In this function I calculate the FRA without subtracting
            % spontaneous firing rate (set it to 0) because I didn't want
            % to have "negative firing rates" that would effect the
            % sparseness results, but rather just take raw number of spikes
            FRA = make_FRA_from_window( raster , protocol_struct , FAKE_SPONT , evoked_win_vec , spont_win_vec ) ;
            % Extract row with best atten
            [ ~ , BF_ind ] = max( max( FRA ) ) ; 
            [ ~ , Best_Atten_ind ] = max( FRA( : , BF_ind ) ) ;
            Best_Atten_FRs = FRA( Best_Atten_ind , : ) ;
            n_freqs = length( Best_Atten_FRs ) ; 
            
            % Calculate S_stat
            S_stat(kk,1) = ( 1 - ( ( sum( Best_Atten_FRs ) / n_freqs )^2 / ( sum( Best_Atten_FRs.^2 ) /n_freqs) ) )/( 1 - ( 1 / n_freqs ) )  ;
            
        end
        
    end
            
end

