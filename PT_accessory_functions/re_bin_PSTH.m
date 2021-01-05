function PSTH_binned = re_bin_PSTH( data_table , old_bin , new_bin , smooth_bin )

    % Re-bin the PSTH from time bins in size "old_bin" to size "new_bin"
    % Smooth the PSTH with a window size "SMOOTH_BIN"
    
    PSTH_binned = cell( size( data_table , 1 ) , 1 ) ; 
    for ii = 1 : size( data_table , 1 )

        PSTH = data_table.PT_PSTH{ ii , 1 } ;
        PSTH_new = re_bin_data( old_bin , new_bin , PSTH ) ;
        PSTH_new = smooth( PSTH_new , smooth_bin )' ;
        PSTH_binned{ ii , 1 } = PSTH_new ;

    end    

end

