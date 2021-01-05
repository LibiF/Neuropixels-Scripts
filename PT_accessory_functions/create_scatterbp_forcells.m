function create_scatterbp_forcells( data_cell , brain_regions , my_colormap , title_text , data_col )

    % This function receives a data_cell and plots a scatter plot
    % representing it also printing summary statistics of the populations
    
    SCALE_RANDOM = 0.05 ;
    SCATTER_SIZE = 18 ;
    LINE_DIMS = 0.2 ;

    if nargin < 5
       
        data_col = 1 ;
        
    end    
    
    h_fig = figure() ;
    h_ax = axes( 'Parent' , h_fig ) ; 
    for kk = 1 : size( data_cell , 1 )
       
        data = data_cell{ kk , 1 } ;
        x_positions = kk + randn( size( data , 1 ) ,1 ).* SCALE_RANDOM ;    
        scatter( h_ax , x_positions , data, SCATTER_SIZE , my_colormap( kk, : ) ) ;
        hold( h_ax, 'on' ) ;
        line( h_ax , [ kk - LINE_DIMS , kk + LINE_DIMS ] , [nanmedian(data), nanmedian(data)] , 'Color' , my_colormap( kk , : ) , 'LineWidth' , 3) ;
        hold( h_ax, 'on' ) ;

        try
        
            disp( ['Region - ' brain_regions{1,kk}  ] )
            disp( ['Mean '  num2str(nanmean( data_cell{ kk , 1 } ) ) ] ) ;
            disp( ['Std '  num2str(nanstd(data_cell{ kk , 1 } ) ) ] ) ;
            disp( ['Median '  num2str(nanmedian( data_cell{ kk , 1 } ) ) ] ) ;
            disp( [ 'Quantiles (25%,75%) ' num2str( quantile( data_cell{ kk , 1 }(:,data_col),[ 0.25 , 0.75 ] ) ) ] ) ;
            disp( [ 'Num of units' num2str( size( data_cell{ kk ,1 } , 1 ) ) ] ) ;
        end
        
    end

    xticks( h_ax , 1 : length( brain_regions ) ) ;
    xticklabels( h_ax , char( brain_regions ) ) ;
    h_ax.FontSize = 12 ;    
    
    title( h_ax , title_text ) ;



end

