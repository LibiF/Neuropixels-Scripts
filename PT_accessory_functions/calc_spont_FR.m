function [spont_FR_mat , spont_FR_std_mat ] = calc_spont_FR(data_table , PRE_STIM_MSEC )

    % This function calculates the spontaneous firing rate in Hz based on
    % baseline activity pre pure tone trials.
    % PRE_STIM_MSEC - the number of baseline ms pre stimulus onset within
    % each trial

    MSEC_IN_SEC = 1000 ;
    spont_FR_mat = zeros( size( data_table , 1 ) , 1 ) ;
    spont_FR_std_mat = zeros( size( data_table , 1 ) , 1 ) ;
    for mm = 1 : size( data_table , 1 )
        
        raster = data_table.PT_raster{mm , 1} ;
        if ~isempty( raster )
        
            spont_FR = sum( raster( : , 1 : PRE_STIM_MSEC ) , 2 ) .* ( MSEC_IN_SEC/PRE_STIM_MSEC) ;
            spont_FR_mat( mm ) = mean( spont_FR ) ;
            spont_FR_std_mat( mm ) = std( spont_FR ) ; 
       
        end
       
    end

end

