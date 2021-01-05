function evoked_FR_mat = calc_BFevoked_FRs( data_table )

    % The function CALC_BFEVOKED_FRS receives a data table and extracts
    % from the FRAs of excited units their evoked FRs in response to each
    % units BF
    
    MSEC_IN_SEC = 1000 ;
    
    evoked_FR_mat = zeros( size( data_table , 1 ) , 1 ) ;
    
    for nn = 1 : size( data_table , 1 )
        
        is_excited = data_table.is_excited( nn , 1 ) ;
        if is_excited
            
            FRA = data_table.FRA{ nn , 1 } ;
            BF_evoked_FR = max( max( FRA ) ) ;
            evoked_FR_mat( nn , 1 ) = BF_evoked_FR ;
        
        else
            
             evoked_FR_mat( nn , 1 ) = nan ;
            
        end
        
    end
 
end
