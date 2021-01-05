function [ log_signif_syll_per_region ,  lin_signif_syll_per_region ] = plot_significance_CDF_per_region( data_table ,...
            brain_regions , my_colormap , resp_mode )

    % The function PLOT_SIGNIFICANCE_CDF_PER_REGION calculates the number
    % of linear and logarithmic FMs evoking a sigificant response. The
    % response type is set by the parameter "resp_mode":
    % resp_mode = 1 - check excitatory responses
    % resp_mode = 0 - check suppressive responses
    % resp_mode = 2 - check either excitatory/inhibitory
    % The function outputs cells with size nx1 where n is the number of
    % brain regions analyzed. Each inset in the cell contains a vector of
    % the number of significant responses for each unit in this region.
    % The function plots a CDF of the results.

    h_fig = figure() ;
    h_ax = axes( 'Parent' , h_fig ) ;
    
    h_fig2 = figure() ;
    h_ax2 = axes( 'Parent' , h_fig2 ) ;
    
    log_signif_syll_per_region = cell( size( brain_regions , 2 ) , 1 ) ;
    lin_signif_syll_per_region = cell( size( brain_regions , 2 ) , 1 ) ;
    
    for kk = 1 : size( brain_regions , 2 )
        
        reg_units = data_table( data_table.acronym == brain_regions{ 1 , kk } , : ) ;
        
        if resp_mode == 1 
        
            log_signifs = cell2mat( reg_units.log_sig_mat ) ;
            lin_signifs = cell2mat( reg_units.lin_sig_mat ) ;
        
        elseif resp_mode == 0
            
            log_signifs = cell2mat( reg_units.log_sig_inh_mat ) ;
            lin_signifs = cell2mat( reg_units.lin_sig_inh_mat ) ;
        
        elseif resp_mode == 2
            log_exc = cell2mat( reg_units.log_sig_mat ) ;
            log_inh = cell2mat( reg_units.log_sig_inh_mat ) ;
            log_signifs = ( log_exc + log_inh ) > 0 ;
            
            lin_exc = cell2mat( reg_units.lin_sig_mat ) ;
            lin_inh = cell2mat( reg_units.lin_sig_inh_mat ) ;
            lin_signifs = ( lin_exc + lin_inh ) > 0 ;
        
        end    
        n_log_stims = size( log_signifs , 2 ) ;
        log_signifs = sum ( log_signifs , 2 ) ;
        
        n_lin_stims = size( lin_signifs , 2 ) ;
        lin_signifs = sum ( lin_signifs , 2 ) ;
        
        log_signif_syll_per_region{ kk ,1 } = log_signifs ;
        lin_signif_syll_per_region{ kk , 1 } = lin_signifs ; 
        
        histogram( h_ax , log_signifs , 'Normalization' , 'cdf' , 'DisplayStyle' , 'stairs' , 'EdgeColor' , my_colormap( kk , : ) ,'LineWidth' , 2.5 , 'BinWidth' , 1 ) ;
        hold( h_ax , 'on' ) ;
        
        histogram( h_ax2 , lin_signifs , 'Normalization' , 'cdf' , 'DisplayStyle' , 'stairs' , 'EdgeColor' , my_colormap( kk , : ) ,'LineWidth' , 2.5 , 'BinWidth' , 1 ) ;
        hold( h_ax2 , 'on' ) ;
        
    end
    
        xlim( h_ax , [ 0 , n_log_stims] ) ;
        xlim( h_ax2 , [ 0 , n_lin_stims ] ) ;
        
        legend( h_ax , brain_regions , 'Location' ,'northwest' ) ;
        legend( h_ax2 , brain_regions , 'Location' ,'northwest' ) ; 
        
        if resp_mode == 1 
            
            title( h_ax , 'Num. of significant excitatory log-FM responses' ) ;
            title( h_ax2 , 'Num. of significant excitatory lin-FM responses' ) ;
            
        elseif resp_mode == 0
            
            title( h_ax , 'Num. of significant inhibitory log-FM responses' ) ;
            title( h_ax2 , 'Num. of significant inhibitory lin-FM responses' ) ;
            
        elseif resp_mode == 2
            
            title( h_ax , 'Num. of significant log-FM responses - E & I ' ) ;
            title( h_ax2 , 'Num. of significant lin-FM responses - E & I ' ) ;
            
        end    
        
end

