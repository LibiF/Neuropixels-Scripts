function hist_counts_cell = plot_noisecorrelation_distribution2( corr_cell , shuffle_cell , brain_regions , suptitle_text , my_colormap, y_limits)

    % The function plots the noise correlation distribution alongside its
    % shuffled distribution
    
    h_fig = figure(); 
    n_brain_regions = length( brain_regions ) ;
    hist_counts_cell = cell( n_brain_regions  , 2 ) ;
    for kk = 1 : n_brain_regions
        
        total_corrs = corr_cell{ kk, 1 }  ;
        total_shuffs = shuffle_cell{ kk , 1} ;
        
        h_ax = subplot( 1 , n_brain_regions  , kk ) ; 
        histogram( h_ax , total_corrs , 'Normalization' , 'probability' , 'FaceColor' , my_colormap(kk, : ) , 'BinWidth' , 0.01) ; 
        hold( h_ax , 'on' ) ;
        histogram( h_ax , total_shuffs , 'Normalization' , 'probability' , 'EdgeColor' , 'r' , 'LineWidth' , 1.5 ,'BinWidth' , 0.01 , 'DisplayStyle' , 'stairs') ; 
        hist_counts =  histcounts(total_corrs , [ -1 : 0.05 : 1 ] ) ; 
        hist_counts = hist_counts ./ sum( hist_counts ) ; 
        hist_counts_cell{ kk , 1 } = hist_counts ; 
        
        xlim( h_ax , [ -1 , 1 ] ) ;
        ylim( h_ax , y_limits ) ;
        xlabel( h_ax , [ 'Noise Correlations' ] , 'FontSize' , 12) ;
        title( h_ax , brain_regions{ 1 , kk } ) ;

        ax_xticks = get(h_ax ,'XTickLabel');  
        set(h_ax,'XTickLabel',ax_xticks,'fontsize',10)
        set(h_ax,'XTickLabelMode','auto')
        
        ax_yticks = get(h_ax ,'YTickLabel');  
        set(h_ax,'YTickLabel',ax_yticks,'fontsize',10)
        set(h_ax,'YTickLabelMode','auto')
        
        
    end       
    
end
