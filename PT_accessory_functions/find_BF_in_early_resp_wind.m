function [FRA_cell , FRA_signif_cell , BF_mat, lat_to_BFfirstspike ] = find_BF_in_early_resp_wind( data_table ,...
         PRE_STIM_MSEC , protocol_struct , window_size , check_excitation )

    % The function FIND_BF_IN_EARLY_RESP_WIND is a function that finds the
    % best frequency for  units with an early excitatory response window. 
    % The function receives the data table, the protocol_struct, the size
    % of the excitatory window and the pre stimulus baseline time 
    % (PRE_STIM_MSEC). The parameter "check_excitation" is set to 1 if we 
    % want to check whether the unit is significantly excitated or 0 if we 
    % don't (and we calculate the FRA anyway).
    % The function outputs:
    % FRA_cell - FRAs for all the (excited) units in the data 
    % FRA_signif_cell - A boolean matrix summarizing the statistical
    % significance of responses for all (frequency x attenuation)
    % combinations 
    % BF_mat - matrix summarizing the BFs of all (excited) units
    % lat_to_BFfirstspike - the latency to the first spike in response to
    % BF trials. This is a nx4 matrix with n=num. of units and each column
    % containing a different statistics of the latency to 1st spike:
    % Col 1 - median latency
    % Col 2 - std latency 
    % Col 3 - CV = std/mean
    % Col 4 - num of trials with a spike

    % Extract protocol parameters
    n_reps = protocol_struct.n_reps ;
    n_attens = length( protocol_struct.atten ) ;
    
    % Create empty data structures 
    FRA_cell = cell( size( data_table , 1 ) , 1 ) ;
    FRA_signif_cell = cell( size( data_table , 1 ) , 1 ) ;
    BF_mat = zeros( size( data_table , 1 ) , 1 ) ;
    lat_to_BFfirstspike = zeros( size( data_table , 1 ) , 4 ) ;
    
    for kk = 1 : size( data_table , 1 ) 
        
        raster = data_table.PT_raster{ kk , 1 } ; 
        
        if check_excitation
            
            is_excitatory = data_table.is_early_excited( kk , 1 ) ;
            if isnan( is_excitatory )
        
                is_excitatory = false ;
                
            end
            
        else
            
            is_excitatory = true ;
            
        end
        
        if is_excitatory
            
            % Extract spontaneous and evoked windows
            spont_FR = data_table.spont_FR( kk , 1 ) ;
            evoked_win_start = data_table.max_early_resp_wind( kk , 1 ) ;
            evoked_win_vec = [ evoked_win_start : evoked_win_start + window_size ] + PRE_STIM_MSEC ;
            spontaneous_window_vec = [ PRE_STIM_MSEC - window_size : PRE_STIM_MSEC ] ; 
            
            % Calculate FRAs
            [ FRA  , FRA_significance ] = make_FRA_from_window2( raster , protocol_struct , spont_FR , evoked_win_vec , spontaneous_window_vec ) ;
            FRA_cell{ kk ,1 } = FRA ;
            FRA_signif_cell{ kk , 1 } = FRA_significance ; 
            
            % Find BFs
            [ ~ , BF_ind ] = max( max( FRA ) ) ;
            [ ~ , atten_ind ] = max( FRA( : , BF_ind ) ) ;
            BF = protocol_struct.freqs( BF_ind ) ; 
            BF_mat( kk , 1 ) = BF ; 
            
            % Find BF trials and extract latency to 1st spike in all trials
            % which contain at least 1 spike
            freq_rows_in_rast =  [ (BF_ind - 1 )*n_reps*n_attens + 1 : BF_ind*n_reps*n_attens ] ;
            freq_rows_in_rast = freq_rows_in_rast( (atten_ind-1)*n_reps+1 : atten_ind*n_reps ) ;
            freq_trials = raster( freq_rows_in_rast , : ) ;
            freq_trials = freq_trials( : , PRE_STIM_MSEC + 1 : end ) ;
            freq_trials = freq_trials( sum( freq_trials , 2 ) > 0 , : ) ;
            [ ~ , first_spikes ] = max(freq_trials,[],2) ; 
            if ~isempty( first_spikes )

                % Calculate latency summary statistics
                lat_to_BFfirstspike( kk , 1 ) = median( first_spikes ) ; 
                lat_to_BFfirstspike( kk , 2 ) = std( first_spikes ) ; 
                lat_to_BFfirstspike( kk , 3 ) =  std( first_spikes ) ./ mean( first_spikes ) ; 
                lat_to_BFfirstspike( kk , 4 ) = length( first_spikes ) ;
                
            else

                lat_to_BFfirstspike( kk , : ) = nan ;
                
            end 
        
        end
        
    end

end

