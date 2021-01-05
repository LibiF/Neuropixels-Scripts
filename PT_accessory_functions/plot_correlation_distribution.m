function hist_counts_cell = plot_correlation_distribution( corr_cell , brain_regions , suptitle_text , my_colormap)

    % This function plots the correlation distributions found in corr_cell
    % alongside its shuffled correlation distribution.
    % The functions assumes the input to be a cell similar to the one
    % outputted by "calculate_sig_corrs3", where the cell contains n rows 
    % and 4 columns:
    % n - corresponds to the number of cross region combinations
    % Columns:
    % 1 - correlations between all units
    % 2 - correlations between simultaneously recorded units
    % 3 - shuffled distribution for all units
    % 4 - shuffled distribution for simultaneously recorded pairs
    
    h_fig = figure(); 
    n_brain_regions = length( brain_regions ) ;
    hist_counts_cell = cell( n_brain_regions  , 2 ) ;
    for kk = 1 : n_brain_regions
        
        total_corrs = corr_cell{ kk, 1 }  ; 
        by_date_corrs = corr_cell{ kk , 2 } ;
        total_shuffled_corrs = corr_cell { kk , 3 } ; 
        total_bydate_shuffled_corrs = corr_cell { kk , 4 } ; 

        h_ax = subplot( 2 , n_brain_regions  , kk ) ; 
        histogram( h_ax , total_corrs , 'Normalization' , 'probability' , 'FaceColor' , my_colormap(kk, : ) , 'BinWidth' , 0.05) ; 
        hold( h_ax ,'on' ) ;
        histogram( h_ax , total_shuffled_corrs , 'Normalization' , 'probability' , 'EdgeColor' , 'r' , 'LineWidth' , 1.5 ,'BinWidth' , 0.05 , 'DisplayStyle' , 'stairs') ; 
        hist_counts =  histcounts(total_corrs , [ -1 : 0.05 : 1 ] ) ; 
        hist_counts = hist_counts ./ sum( hist_counts ) ; 
        hist_counts_cell{ kk , 1 } = hist_counts ; 
        
        xlim( h_ax , [ -1 , 1 ] ) ;
        ylim( h_ax , [ 0 , 0.4 ] ) ;
        xlabel( h_ax , [ 'Signal Correlations' ] ) ;
        title( h_ax , brain_regions{ 1 , kk } ) ;
        h_ax2 = subplot( 2 , n_brain_regions  , kk+n_brain_regions ) ;
        histogram( h_ax2 , by_date_corrs , 'Normalization' , 'probability' , 'FaceColor' , my_colormap(kk, : ) , 'BinWidth' , 0.05 ) ; 
        hold( h_ax2 ,'on' ) ;
        histogram( h_ax2 , total_bydate_shuffled_corrs , 'Normalization' , 'probability' , 'EdgeColor' , 'r' , 'LineWidth' , 1.5 ,'BinWidth' , 0.05 , 'DisplayStyle' , 'stairs') ; 
        xlim( h_ax2 ,  [ -1 , 1 ] ) ;
        xlabel( h_ax2 , 'S. C. by date' )
        ylim( h_ax2 , [ 0 , 0.4 ] ) ;
%         suptitle( suptitle_text ) ;
        
        hist_counts_by_date =  histcounts( by_date_corrs , [ -1 : 0.05 : 1 ] ) ; 
        hist_counts_by_date = hist_counts_by_date ./ sum( hist_counts_by_date ) ; 
        hist_counts_cell{ kk , 2 } = hist_counts_by_date ; 

    end       
    
end
