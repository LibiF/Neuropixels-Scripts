function create_BFbp_from_cells_log_scale( window_cell , brain_regions ,...
         my_colormap , title_text , protocol_struct , data_col)

    % This function creates a logarithmically scale box-plot of BF
    % distributions
    
    kHZ_TO_HERTZ = 1000 ;
    freqs = protocol_struct.freqs .*  kHZ_TO_HERTZ ;
    tick_dist = logspace( log( freqs(1) )/log(10) , log( freqs(end) )/log(10) , 10 ) ;
    tick_labels = round( (tick_dist(1:end)' ./ kHZ_TO_HERTZ ) .* 100 ) ./ 100  ;

    if nargin < 6
       
        data_col = 1 ;
        
    end    
    
    data = [] ;
    groups = [] ;
    for kk = 1 : size( window_cell , 1 )
       
        data = [ data ; window_cell{ kk , 1 } ];
        groups = [ groups ; kk.*ones( size( window_cell{kk,1} , 1 ) , 1 ) ] ;

%         disp( nanmedian( window_cell{ kk , 1 } ) ) ;
%         disp( iqr( window_cell{ kk , 1} ) ) ;
        disp( ['Region - ' brain_regions{1,kk}  ] )
        disp( ['Mean '  num2str(mean( window_cell{ kk , 1 } ) ) ] ) ;
        disp( ['Std '  num2str(std(window_cell{ kk , 1 } ) ) ] ) ;
        disp( ['Median '  num2str(median( window_cell{ kk , 1 } ) ) ] ) ;
        disp( [ 'Quantiles (25%,75%) ' num2str( quantile( window_cell{ kk , 1 },[ 0.25 , 0.75 ] ) ) ] ) ;
        disp( [ 'Num of units' num2str( size( window_cell{ kk ,1 } , 1 ) ) ] ) ;
        
    end
    
    h_fig = figure() ;
    h_ax = axes( 'Parent' , h_fig ) ; 
    h_bp = boxplot( h_ax , log( data(:,data_col).*  kHZ_TO_HERTZ)./log(10), groups , 'Labels' , brain_regions , 'symbol','+' ) ;
    h_boxes = findobj(h_ax,'Tag','Box'); 
    for j=1:length(h_boxes) 

        patch(get(h_boxes(j),'XData'),get(h_boxes(j),'YData'),my_colormap(length(h_boxes) + 1 - j,:),'FaceAlpha',.5);

    end 
    
    title( h_ax , title_text ) ;
    yticks( h_ax , log(tick_dist)./log(10) ) ;
    yticklabels( h_ax ,  tick_labels ) ;
    ylabel( h_ax , 'BF [kHz]' , 'FontSize' , 14 ) ;
    h_ax.FontSize = 12 ;
    
end