function create_bp_from_cells( data_cell , brain_regions , my_colormap , title_text , data_col)

    % The function CREATE_BP_FROM_CELLS receives a data_cell and plots the
    % summary of the data as boxplots such that each box corresponds to a
    % single brain region.
    % If the data cell contains more than 1 column, the parameter data_col
    % is used to extract the relevant data
    
    if nargin < 5
       
        data_col = 1 ;
        
    end 
    
    data = [] ;
    groups = [] ;
    
    for kk = 1 : size( data_cell , 1 )
       
        data = [ data ; data_cell{ kk , 1 } ];
        groups = [ groups ; kk.*ones( size( data_cell{kk,1} , 1 ) , 1 ) ] ;

        try
        
            disp( ['Region - ' brain_regions{1,kk}  ] )
            disp( ['Mean '  num2str(nanmean( data_cell{ kk , 1 } ) ) ] ) ;
            disp( ['Std '  num2str(nanstd(data_cell{ kk , 1 } ) ) ] ) ;
            disp( ['Median '  num2str(nanmedian( data_cell{ kk , 1 } ) ) ] ) ;
            disp( [ 'Quantiles (25%,75%) ' num2str( quantile( data_cell{ kk , 1 }(:,data_col),[ 0.25 , 0.75 ] ) ) ] ) ;
            disp( [ 'Num of units' num2str( size( data_cell{ kk ,1 } , 1 ) ) ] ) ;
        
        end
        
    end
    
    h_fig = figure() ;
    h_ax = axes( 'Parent' , h_fig ) ; 
    h_bp = boxplot( h_ax , data(:,data_col) , groups , 'Labels' , brain_regions , 'symbol','+' ) ;
    h_boxes = findobj(h_ax,'Tag','Box'); 
    for j=1:length(h_boxes) 

        patch(get(h_boxes(j),'XData'),get(h_boxes(j),'YData'),my_colormap(length(h_boxes) + 1 - j,:),'FaceAlpha',.5);

    end 
    title( h_ax , title_text ) ;

end

