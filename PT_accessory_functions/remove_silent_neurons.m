function data_table = remove_silent_neurons( data_table , thresh )

    % This function removes from the data table all units with less than
    % "thresh" spikes during entire pure tone protocol

    exclude_list = [] ;
    for kk = 1 : size( data_table , 1 )
        
        sum_spikes = sum( sum( data_table.PT_raster{kk,1} ) ) ;
        if sum_spikes < thresh
           
            exclude_list = [ exclude_list , kk ] ;
            
        end    
        
    end
    
    data_table( exclude_list , : ) = [] ;
    
end

