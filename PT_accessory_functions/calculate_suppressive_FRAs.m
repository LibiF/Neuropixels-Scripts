function [FRA_cell , FRA_signif_cell ] = calculate_suppressive_FRAs( data_table ,...
          PRE_STIM_MSEC , protocol_struct , window_size )

    % The function CALCULATE_SUPPRESSIVE_FRAs receives a data table, the
    % pure tone protocol struct, the time of baseline activity pre stimulus
    % (PRE_STIM_MSEC) and the window size and calculates an FRA for units
    % which are suppressed 
    
    % Extract parameters from protocol
    n_reps = protocol_struct.n_reps ;
    n_attens = length( protocol_struct.atten ) ;
 
    % Create empty data structures
    FRA_cell = cell( size( data_table , 1 ) , 1 ) ;
    FRA_signif_cell = cell( size( data_table , 1 ) , 1 ) ;
    
    for kk = 1 : size( data_table , 1 ) 
        
        raster = data_table.PT_raster{ kk , 1 } ; 
        is_suppressed = data_table.is_inhibited( kk , 1 ) ;

        
        if is_suppressed > 0 
            
            suppressed_win_start = data_table.min_resp_wind( kk , 1 ) ;
            suppresed_win_vec = [ suppressed_win_start : suppressed_win_start + window_size ] + PRE_STIM_MSEC ;
            spontaneous_window_vec = [ PRE_STIM_MSEC - window_size : PRE_STIM_MSEC ] ; 
            [ FRA  , FRA_significance ] = make_FRA_from_suppressive_window( raster , protocol_struct , suppresed_win_vec , spontaneous_window_vec ) ;
            FRA_cell{ kk ,1 } = FRA ;
            FRA_signif_cell{ kk , 1 } = FRA_significance ; 
            
        end
        
    end

end

