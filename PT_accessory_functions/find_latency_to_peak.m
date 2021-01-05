function [ latency_to_peak ] = find_latency_to_peak( data_table , PRE_STIM_MSEC , SMOOTH_WINDOW ) 

    % This function finds the latency to the peak of each units PSTH
    % PSTH is smoothes using a window size "SMOOTH_WINDOW"
    
    Fs = 1000 ;
    latency_to_peak = nan( size( data_table , 1 ) , 1 ) ;

    for kk = 1 : size( data_table , 1 )

        if data_table.is_excited( kk , 1 ) == 1
            
            unit_PSTH = data_table.PT_PSTH{ kk , 1 } ;
            smooth_PSTH = smooth( unit_PSTH , SMOOTH_WINDOW ) ;   
            smooth_PSTH = smooth_PSTH( PRE_STIM_MSEC + 1 : end ) ;
            [ ~ , latency_to_peak( kk  ,1 ) ] = max( smooth_PSTH(1:end-1) ) ;
            
        end
        
    end
    
end

