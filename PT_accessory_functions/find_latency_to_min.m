function [ latency_to_min ] = find_latency_to_min( data_table , PRE_STIM_MSEC , SMOOTH_WINDOW ) 

    % This function calculates the latency to the minima of each units
    % PSTH. PSTHs are smoothed using a window of size "SMOOTH_WINDOW".

    Fs = 1000 ;
    latency_to_min = nan( size( data_table , 1 ) , 1 ) ;

    for kk = 1 : size( data_table , 1 )

        if data_table.is_excited( kk , 1 ) == 0 && data_table.is_inhibited( kk , 1 ) == 1
            
            unit_PSTH = data_table.PT_PSTH{ kk , 1 } ;
            smooth_PSTH = smooth( unit_PSTH , SMOOTH_WINDOW ) ;
            smooth_PSTH = smooth_PSTH( PRE_STIM_MSEC + 1 : end ) ;
            [ ~ , latency_to_min( kk  ,1 ) ] = min( smooth_PSTH(1:end-1) ) ;
            
        end
        
    end
    
end

