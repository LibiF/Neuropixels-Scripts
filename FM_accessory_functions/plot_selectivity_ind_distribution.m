function [total_selec_ind_log , total_selec_ind_lin ,...
          total_DSI_allstims ] = plot_selectivity_ind_distribution( data_table ,...
          brain_regions , my_colormap )

    % The function PLOT_SELECTIVITY_IND_DISTRIBUTION receives the data
    % table after DSIs were calculated and plots the DSI distributions of
    % different brain regions.
    % It also outputs the distributions taking into consideration only FM
    % responsive units

    total_selec_ind_log = cell( size( brain_regions , 2 ) , 1 ) ;
    total_selec_ind_lin = cell( size( brain_regions , 2 ) , 1 ) ;
    total_DSI_allstims =  cell( size( brain_regions , 2 ) , 1 ) ;
    
    h_fig = figure() ;
    h_ax = axes('Parent' , h_fig ) ;
    
    h_fig2 = figure() ;
    h_ax2 = axes('Parent' , h_fig2 ) ;
    
    h_fig3 = figure() ;
    h_ax3 = axes('Parent' , h_fig3 ) ;

    for kk = 1 : size( brain_regions , 2 )
        
        reg_units = data_table( data_table.acronym == brain_regions{ 1 , kk } , : ) ;
        active_reg_units = reg_units( reg_units.is_excited_lin == 1 | reg_units.is_inhibited_lin == 1 | reg_units.is_excited_log == 1 | reg_units.is_inhibited_log == 1 , : ) ;
        log_reg_units = reg_units( reg_units.is_excited_log == 1 | reg_units.is_inhibited_log == 1 , : ) ;
        lin_reg_units = reg_units( reg_units.is_excited_lin == 1 | reg_units.is_inhibited_lin == 1 , : ) ;
        
        total_selec_ind_log{ kk , 1 } = log_reg_units.total_selectivity_index_log ;
        total_selec_ind_lin{ kk , 1 } = lin_reg_units.total_selectivity_index_lin ;
        total_DSI_allstims{ kk , 1 } = active_reg_units.total_DSI ;

        histogram( h_ax ,  log_reg_units.total_selectivity_index_log , 'EdgeColor' , my_colormap( kk , : ) , 'Normalization' ,  'CDF' , 'BinWidth' , 0.1 , 'DisplayStyle' , 'stairs' ) ;
        hold( h_ax , 'on' );
        
        histogram( h_ax2 ,  lin_reg_units.total_selectivity_index_lin , 'EdgeColor' , my_colormap( kk , : ) , 'Normalization' ,  'CDF' , 'BinWidth' , 0.1, 'DisplayStyle' , 'stairs' ) ;
        hold( h_ax2 , 'on' );
        
        histogram( h_ax3 ,  active_reg_units.total_DSI , 'EdgeColor' , my_colormap( kk , : ) , 'Normalization' , 'CDF' , 'BinWidth' , 0.1  , 'DisplayStyle' , 'stairs') ;
        hold( h_ax3 , 'on' );
        
    end    

    title( h_ax , 'Total log S.I' ) ;
    title( h_ax2 , 'Total lin S.I' ) ;
    title( h_ax3 , 'DSI - all FMs' ) ;
    
    xlim( h_ax , [ -1 , 1 ] ) ;
    xlim( h_ax2 , [ -1 , 1 ] ) ;
    xlim( h_ax3 , [ -1 , 1 ] ) ;
    
end

