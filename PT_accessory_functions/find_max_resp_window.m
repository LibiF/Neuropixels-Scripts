function [ best_window_mat , is_excited ] = find_max_resp_window( data_table , PRE_STIM_MSEC ,...
          window_size , response_end_lim )

    % The function FIND_MAX_RESP_WINDOW looks for a temporal window of size
    % "window_size" in which excitation post stimulus is maximal.
    % The window is searched in the limited time window between stimulus
    % onset and up to "response_end_lim".
    % After the window is detected, the significance of the excitatory
    % response is tested via ttest comparing to baseline activity
    %
    % The function outputs:
    % best_window_mat - Matrix of the onset time of the maximal excitatory
    % window for all units in the data_table
    % is_excited - Bollean vector size nx1 showing the results of the 
    % performed ttest
      
    best_window_mat = zeros( size( data_table , 1 ) , 1 ) ;  
    is_excited = zeros( size( data_table , 1 ) , 1 ) ; 
    for kk = 1 : size( data_table , 1 ) 

        PSTH = data_table.PT_PSTH{ kk ,1 } ;
        raster = data_table.PT_raster{ kk , 1 } ;
        PSTH = PSTH( PRE_STIM_MSEC + 1 : PRE_STIM_MSEC + response_end_lim ) ;
        response_sum = zeros( 1 , response_end_lim - window_size ) ;
        for mm = 1 : response_end_lim - window_size

            response_sum( 1 , mm ) = sum( PSTH( mm : mm + window_size ) ) ; 
            
        end
        
        [ ~ , max_window ] = max( response_sum ) ;
        best_window_mat( kk , 1 ) = max_window ;
        trials_sum = sum( raster( : , PRE_STIM_MSEC +  max_window : PRE_STIM_MSEC + max_window + window_size ) , 2 ) ;
        baseline_sum = sum( raster( : , PRE_STIM_MSEC - window_size : PRE_STIM_MSEC )  , 2 ) ;
        [ is_excited( kk ,1 ) , ~ ] = ttest2( trials_sum , baseline_sum , 'Tail' ,'right' ) ;
        
    end
    
end

