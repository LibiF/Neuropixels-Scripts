function [ ttest_scores , signif_mat_cell ] = compare_speedselectivity_via_ttest( data_cell ,...
           brain_regions , title_text ) 

    % The function COMPARE_SPEEDSELECTIVITY_VIA_TTEST receives a data cell
    % with n-rows = number of brain regions compared and 1 column. In
    % every inset of the data cell we find the matrix with rows
    % corresponding the number of units in every region and columns
    % correspond to responses to FMs of different speeds/velocities.
    % The function statistically compares between population responses of
    % different regions (output variable "ttest_scores") and between the
    % population response to different speeds within regions (output
    % variable signif_mat_cell).
    
    % Combinations of cross-region comparisons
    combins = nchoosek( 1 : size(data_cell , 1 ) , 2 ) ; 
    
    ttest_scores = zeros( size( combins , 1 ) , size( data_cell(1,1) , 2 ) ) ; 
    
    for cc = 1 : size( combins , 1 ) 
        
        columns = combins( cc , : ) ;
        pop1 = data_cell{ columns( 1 ) , 1} ;
        
        pop2 = data_cell{ columns( 2 ) , 1 } ; 
        
        for tt = 1 : size( pop1 , 2 ) 
        
            [ ~ , ttest_scores( cc , tt ) ] = ttest2( pop1( : , tt ) , pop2( : , tt ) ) ;
                                   
        end    
            
        h_fig = figure() ;
        h_ax = axes( 'Parent' , h_fig ) ;
        imagesc( h_ax , ttest_scores( cc , : ) ) ;
        
        % These following rows can be used to make the color bar log-scaled
%         caxis( h_ax , [ 1e-10 , 0.1] ) ;
%         decimals = (-10:1:-1) ;
%         ticks = 10.^decimals ;
%         cb = colorbar('Ticks',ticks, 'TickLabels',num2str( (log(ticks)./log(10))' ) ) ;
%         cb = colorbar
%         set( h_ax, 'ColorScale' ,'log' ) ;
%         cb = colorbar  ;
%         cb.Ticks = log(ticks)./log(10) ;
%         cb.TickLabels = decimals ;
        title( h_ax , [ brain_regions( columns(1) ) ' vs. ' brain_regions( columns( 2 ) ) ] ) ;
        colormap( bone ) 
        
    end
    
    signif_mat_cell = cell( size( data_cell , 1 ) , 1 ) ;
    
    for kk =1 : size( data_cell , 1 )
    
        region_data = data_cell{ kk , 1} ;
        signif_mat = 0.5.*ones( size( region_data , 2 ) ) ;
        for ii = 1 : size( region_data, 2 ) - 1
        
            data_point = region_data( : , ii ) ;
            
            for jj = ii + 1 : size( region_data , 2 )
            
                data_point2 = region_data( : , jj ) ;
                [ ~ , signif_mat( ii , jj ) ] = ttest2( data_point , data_point2 ) ;
                
            end    
                
        end
        
        signif_mat = signif_mat + signif_mat'  - 0.5;
        signif_mat_cell{ kk , 1 } = signif_mat ; 
        h_fig = figure() ;
        h_ax = axes( 'Parent' , h_fig ) ;
        imagesc( h_ax , signif_mat ) ;
        title( h_ax , [ title_text , brain_regions( kk ) ] ) ;
        colormap( bone ) 
        caxis( h_ax , [ 0 , 0.05] ) ;
        
    end
    
end

