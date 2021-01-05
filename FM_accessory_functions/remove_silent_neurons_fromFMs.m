function data_table = remove_silent_neurons_fromFMs( data_table , thresh , columns)

    % The function REMOVE_SILENT_NEURONS_FROMFMS receives a data table and
    % removes from it all units with missing FM protocol or less than
    % "thresh" spikes fired during the FM protocol
    
    exclude_list = [] ;
    for kk = 1 : size( data_table , 1 )
        
        unit_column_cell = cell( length( columns ) , 1 ) ;
        for mm = 1 : length( columns ) 
            
            unit_column_cell{ mm , 1 } = data_table{ kk , columns(mm)}{1,1}  ;
                        
        end
        
        test_sum = cellfun( @(x) sum( sum( x ) ) , unit_column_cell ) ;
        if any( test_sum < thresh )
           
            exclude_list = [ exclude_list , kk ] ;
            
        end    
        
    end
    
    data_table( exclude_list , : ) = [] ;
    
end
