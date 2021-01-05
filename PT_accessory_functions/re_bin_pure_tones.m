function data_table = re_bin_pure_tones( data_table , old_bin , new_bin )

    % This function re-bins raster plots in temporal bins size "new_bin"
    % insted of "old_bin".
    % The binned rasters are added as the last column of the data_table
    % suuplied
    
    rasters_binned = cell( size( data_table , 1 ) , 2 ) ;
    
for ii = 1 : size( data_table , 1 )
    
    rast = data_table.PT_raster{ ii , 1 } ;
    
    for kk = 1 : size( rast , 1 )
    
       rast_new( kk , : ) = re_bin_data( old_bin , new_bin , rast( kk , : ) ) ;

    end    
     
    rasters_binned{ ii , 1 } = rast_new ;
    rasters_binned{ ii , 2 } = new_bin ;
    
end    

data_table.binned_PT_rast = rasters_binned ;

end

