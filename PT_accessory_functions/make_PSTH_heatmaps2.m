function make_PSTH_heatmaps2( data_table  , brain_regions , PRE_STIM_BINS , bin_size , time_frame )

    % The function plots PSTH heatmaps ordered by latencies to peak (top
    % part) and latencies to minima (bottom part) for all brain_regions

    RED = [	174, 96, 88 ] ./ 256 ; 
    BLUE = [ 1,85,151] ./ 256 ; 

    h_fig = figure() ; 
    
    for kk = 1 : size( brain_regions , 2 )

        reg_units = data_table( data_table.acronym == brain_regions{ 1 , kk } , : ) ;
        excited_units = reg_units( reg_units.is_excited == 1 , : ) ;
        excited_PSTH_mat = cell2mat( excited_units.PSTH_binned ) ;
        inhibited_units = reg_units( ( reg_units.is_excited == 0 & reg_units.is_inhibited == 1 ) , : ) ;
        inhibited_PSTH_mat = cell2mat( inhibited_units.PSTH_binned ) ;
        excited_PSTH_mat = excited_PSTH_mat( : , PRE_STIM_BINS + 1 : end ) ;
        inhibited_PSTH_mat =  inhibited_PSTH_mat( : , PRE_STIM_BINS + 1 : end ) ;

        normalized_excited_PSTH = excited_PSTH_mat ./ max( excited_PSTH_mat , [] , 2 ) ;
        [ ~ , max_ind ] = max( normalized_excited_PSTH , [] , 2 ) ;
        [ sorted_max_exc , sorted_order_exc ] = sortrows( max_ind ) ;
        normalized_excited_PSTH = normalized_excited_PSTH( sorted_order_exc , : ) ;
        
        normalized_inhibited_PSTH = inhibited_PSTH_mat./ max( inhibited_PSTH_mat , [] , 2 ) ;
        [ ~ , min_ind ] = min( normalized_inhibited_PSTH , [] , 2 ) ;
        [ sorted_min_inh , sorted_order_inh ] = sortrows( min_ind ) ;
        normalized_inhibited_PSTH = normalized_inhibited_PSTH( flipud( sorted_order_inh ) , : ) ;
        
        total_PSTH = [ normalized_excited_PSTH ; normalized_inhibited_PSTH ] ;
        
        h_ax = subplot( 1 , size( brain_regions , 2 ) , kk ) ;
        imagesc( h_ax , total_PSTH ) ;
        xticks( h_ax , [ 1 : bin_size : size( total_PSTH , 2 ) ] ) ;
%         xticks( h_ax , [] ) ;
        xticklabels( h_ax , [ time_frame(1) : 100 : time_frame(end) ]' ) ;
%         yticks( h_ax, [] ) ; 
        xtickangle( h_ax, 90 ) ;
        xlabel( h_ax, 'Time [msec]' ) ;
        ylabel( h_ax, 'Unit num' ) ;
        title( h_ax , brain_regions( 1 , kk ) ) ;
        hold( h_ax , 'on' ) ;
        plot( h_ax , sorted_max_exc , 1 : length( sorted_max_exc ) , 'Color' ,  RED , 'LineWidth' , 2  , 'LineStyle' ,'--' ) ;
        hold( h_ax , 'on' ) ;
        plot( h_ax , flipud(sorted_min_inh ) , [ ( length( sorted_max_exc ) + 1 ): length( total_PSTH ) ] , 'Color' , BLUE , 'LineWidth' , 2  , 'LineStyle' ,'--') ;
        set (h_ax,'Ydir','reverse')
        ylim(h_ax, [1 , length( total_PSTH )])
%         set(h_ax, 'box','off' )

    end    
        
end