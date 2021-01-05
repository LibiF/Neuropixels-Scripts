function [ signif_FMs_per_region ] = plot_significance_CDF_per_region_allstimuli( data_table ,...
           brain_regions , my_colormap )
       
    % The function PLOT_SIGNIFICANCE_CDF_PER_REGION_ALLSTIMULI calculates 
    % the number of FMs (both linear and logarithmic) evoking a sigificant 
    % excitatory response. The function outputs a cell with size nx1 where
    % n is the number of brain regions analyzed. Each inset in the cell 
    % contains a vector of the number of significant responses for each 
    % unit in this region.
    % The function plots a CDF of the results.   

    h_fig = figure() ;
    h_ax = axes( 'Parent' , h_fig ) ;
    
    signif_FMs_per_region = cell( size( brain_regions , 2 ) , 1 ) ;

    for kk = 1 : size( brain_regions , 2 )
        
        reg_units = data_table( data_table.acronym == brain_regions{ 1 , kk } , : ) ;
                
        log_signifs = cell2mat( reg_units.log_sig_mat ) ;
        lin_signifs = cell2mat( reg_units.lin_sig_mat ) ;
        total_signifs = [ log_signifs , lin_signifs ] ; 

        n_log_stims = size( log_signifs , 2 ) ;                
        n_lin_stims = size( lin_signifs , 2 ) ;
        n_total_stims = n_log_stims + n_lin_stims ;
        all_signifs = sum ( total_signifs , 2 ) ;

        signif_FMs_per_region{ kk ,1 } = all_signifs ;
        
        histogram( h_ax , all_signifs , 'Normalization' , 'cdf' , 'DisplayStyle' , 'stairs' , 'EdgeColor' , my_colormap( kk , : ) ,'LineWidth' , 2.5 , 'BinWidth' , 1 ) ;
        hold( h_ax , 'on' ) ;
        
    end
    
    xlim( h_ax , [ 0 , n_total_stims ] ) ;     
    ylim( h_ax, [ 0 , 1 ] ) ; 
    legend( h_ax , brain_regions , 'Location' ,'northwest' ) ;
    title( h_ax , 'Significant excitatory FM responses per unit' ) ;
    ylabel( h_ax , 'Probability' ) ;
    xlabel( h_ax , 'Significant responses per unit' ) ;
        
end

