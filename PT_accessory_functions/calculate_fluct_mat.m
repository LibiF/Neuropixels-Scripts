function fluctuation_cell = calculate_fluct_mat( data_table , protocol_example , binned_rast_column )

    % The function CALCULATE_FLUCT_MAT takes creates fluctuations matrices
    % showing the trial-by-trial variability in responses to the different
    % stimuli.
    % The function takes binned raster plots, and calculates the mean
    % response to each stimulus (freqxency x attenuation combination) and
    % subtracts it from each trial to obtain "fluctuation vectors".
    % Essentially, each trials in the fluctuation matrix contains the delta
    % of the response in each time bin and the mean response in this time
    % point for the given stimulus 
    
    fluctuation_cell = cell( size( data_table , 1 ) , 1 ) ;
    n_reps = protocol_example.n_reps ;
    
    for ii = 1 : size( fluctuation_cell , 1 )
        
        raster = data_table{ ii , binned_rast_column }{1,1} ;
        fluct_mat = [] ;
        index_vec = 1 : n_reps : size( raster , 1 ) ;
        
        for kk = 1 : length( index_vec ) 
            
            stim_indexes = [ index_vec( kk ) : 1 : ( index_vec( kk ) + n_reps - 1 ) ] ;
            stim_vecs = raster( stim_indexes , : ) ;
            stim_mean = mean( raster( stim_indexes , : ) , 1 ) ;
            fluct_vecs = stim_vecs - stim_mean ;
            fluct_mat = [ fluct_mat ; fluct_vecs ] ;
            
        end    
        
        fluctuation_cell{ ii , 1 } = fluct_mat ;
        
    end    

end
