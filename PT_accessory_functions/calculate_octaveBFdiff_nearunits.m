function diff_sum = calculate_octaveBFdiff_nearunits( data_table , channel_radius , my_colormap )
    
    % The function CALCULAT_OCTAVEBFDIFF_NEARUNITS find the BF diffence
    % between units recorded from contacts =< channel_radius apart from one
    % another. The function discriminates between two groups of such unit
    % pairs - units recorded from the same brain region and unit pairs
    % recorded from different brain regions.
    % The function plots a summary scatter showing BF distances in octave
    % scale for the two compared groups.

    data_table = data_table( data_table.BF > 0 , : ) ;
    rec_dates = unique( data_table.rec_date ) ;
    
    % Create empty data structures
    BF_diff_same_reg = [] ;
    BF_diff_diff_reg = [] ; 
    
    for mm = 1 : size( rec_dates )

        date_table = data_table( data_table.rec_date == rec_dates(mm , : ) , : ) ;
        % Find all possible unit pairs
        unit_combinations =nchoosek( 1: size(date_table ,1 ) , 2 ) ;
        % Find channels of all unit pairs
        combination_channels = zeros( size( unit_combinations , 1 ) , 2 ) ;
        combination_channels( : , 1 ) = date_table.clus_channel( unit_combinations(: ,1) ) ;
        combination_channels( : , 2 ) = date_table.clus_channel( unit_combinations(: ,2) ) ;
        combination_regions = zeros( size( unit_combinations , 1 ) , 2 ) ;
        % Find brain region affiliation for all unit pairs
        combination_regions( : , 1 ) = grp2idx( date_table.acronym( unit_combinations( : , 1 ) , : ) ) ;
        combination_regions( : , 2 ) = grp2idx( date_table.acronym( unit_combinations( : , 2 ) , : ) ) ;
        
        % Define all "near by" pairs
        is_nearby = abs( combination_channels(:,2) - combination_channels(:,1) ) <= channel_radius ;
        nearby_pairs = unit_combinations( is_nearby , : ) ;
        nearby_channels = combination_channels( is_nearby , : ) ;
        nearby_regions = combination_regions( is_nearby , : ) ; 
        
        % Divide pairs to same/different brain regions
        is_same_reg = nearby_regions(:,1) == nearby_regions(:,2) ; 
        nearby_pairs_same_reg = nearby_pairs( is_same_reg , : ) ;
        nearby_pairs_diff_reg = nearby_pairs( not( is_same_reg) , : ) ;
        
        nearby_channels_same_reg = nearby_channels( is_same_reg , : ) ;
        nearby_channels_diff_reg = nearby_channels( not( is_same_reg ) , : ) ;
        
        % Calculate BF differences for unit pairs from the same region
        for same = 1 : size( nearby_pairs_same_reg , 1 )
            
            unit_num1 = nearby_pairs_same_reg( same , 1 ) ;
            unit_num2 = nearby_pairs_same_reg( same , 2 ) ;
            
            BF1 = date_table.BF( unit_num1 ) ;
            BF2 = date_table.BF( unit_num2 ) ;
            
            largeBF = max( BF1 , BF2 ) ;
            smallBF = min( BF1 , BF2 ) ;
            
            channel1 = date_table.clus_channel( unit_num1 ) ;
            channel2 = date_table.clus_channel( unit_num2 ) ;
            
            delta_BF = log2( largeBF/ smallBF ) ;
            delta_channs = abs( channel2 - channel1 ) ;
            BF_diff_same_reg = [ BF_diff_same_reg ; delta_BF ] ;
          
        end
        
        % Calculate BF differences for pairs from different regions
        for diff = 1 : size( nearby_pairs_diff_reg , 1 )
            
            unit_num1 = nearby_pairs_diff_reg( diff , 1 ) ;
            unit_num2 = nearby_pairs_diff_reg( diff , 2 ) ;
            
            BF1 = date_table.BF( unit_num1 ) ;
            BF2 = date_table.BF( unit_num2 ) ;
            
            largeBF = max( BF1 , BF2 ) ;
            smallBF = min( BF1 , BF2 ) ;
            
            channel1 = date_table.clus_channel( unit_num1 ) ;
            channel2 = date_table.clus_channel( unit_num2 ) ;
            
            delta_BF =  log2( largeBF/ smallBF ) ;
            delta_channs = abs( channel2 - channel1 ) ;
            BF_diff_diff_reg = [ BF_diff_diff_reg ; delta_BF ] ;

            
        end    

    end    
        
    diff_sum{1,1} = BF_diff_same_reg ;
    diff_sum{2,1} = BF_diff_diff_reg ;
    
    % Plot summary scatter plot
    title_text = 'BF differences of nearby units' ; 
    create_scatterbp_forcells( diff_sum(:,1) , {'Within reg' , 'Across reg'} , my_colormap , title_text )
    
    
end

