function ttest_scores = compare_PSTHs_via_ttest( data_cell , brain_regions ) 

    % The function COMPARE_PSTHS_VIA_TTEST recieves a data cell each
    % contatining the population response in k brain regions (num.
    % corresponding to the length of "brain_regions") for
    % n=size(data_cell,1) types of stimuli.
    % The function statistically compares the population response between
    % the each two brain regions and plot a heat-map of the p-values over
    % time
    
    combins = nchoosek( 1 : size(data_cell , 2 ) , 2 ) ; 
    
    ttest_scores = zeros( size( combins , 1 ) , size( data_cell(1,1) , 2 ) ) ; 
    
    for cc = 1 : size( combins , 1 ) 
        
        columns = combins( cc , : ) ;
        pop1 = data_cell( 1 : end , columns( 1 ) ) ;
        pop1 = cell2mat( pop1 ) ;
        
        pop2 = data_cell( 1: end , columns( 2 ) ) ; 
        pop2 = cell2mat( pop2 ) ; 
        
        for tt = 1 : size( pop1 , 2 ) 
        
            [ ~ , ttest_scores( cc , tt ) ] = ttest2( pop1( : , tt ) , pop2( : , tt ) ) ;
                                   
        end    
            
        h_fig = figure() ;
        h_ax = axes( 'Parent' , h_fig ) ;
        imagesc( h_ax , ttest_scores( cc , : ) ) ;
        caxis( h_ax , [ 1e-10 , 0.1] ) ;
        decimals = (-10:1:-1) ;
        ticks = 10.^decimals ;
%         cb = colorbar('Ticks',ticks, 'TickLabels',num2str( (log(ticks)./log(10))' ) ) ;
%         cb = colorbar
        set( h_ax, 'ColorScale' ,'log' ) ;
        cb = colorbar  ;
%         cb.Ticks = log(ticks)./log(10) ;
%         cb.TickLabels = decimals ;
        title( h_ax , [ brain_regions( columns(1) ) ' vs. ' brain_regions( columns( 2 ) ) ] ) ;
        colormap( bone ) 
        
    end
    
end

