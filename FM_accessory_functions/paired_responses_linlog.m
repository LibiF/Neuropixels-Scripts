function difference_vecs = paired_responses_linlog( data_table , brain_regions , my_colormap )

    % The function PAIREd_RESPONSES_LINLOG receives the data table, a list
    % of brain regions and a colormap.
    % For each brain region, the function goes over all units and finds the
    % difference in the number of linear FMs and logarithmic FM stimuli 
    % evoking a significant excitatory response. 
    % The function plots a scatter plot with the summary of results and
    % output "differene_vecs" - a cell the size of "brain_region" with each
    % inset containing a vector of the delta between linear and logarithmic
    % FM significant responses for each unit in a given brain region.

    symbols = [ 's-' ; 'o-' ; '^-' ; '*-' ] ;

    h_fig = figure() ;
    h_ax = axes( 'Parent' , h_fig ) ;
    
    difference_vecs = cell( size( brain_regions , 2 ) ,1 ) ;
    
    for mm = 1 : size( brain_regions , 2 )
        
        reg_units = data_table( data_table.acronym == brain_regions{ 1 , mm } , : ) ;
        signif_lins = sum( cell2mat( reg_units.lin_sig_mat ) , 2 ) ;
        signif_logs = sum( cell2mat( reg_units.log_sig_mat ) , 2 ) ; 
        locs = randn( size( signif_logs , 1 ) , 1 ) .* 0.05 + mm *2 ;
        
        for aa = 1 : size( signif_lins , 1 ) 
            
            plot( h_ax , [ locs(aa) , locs(aa) + 1 ] , [ signif_lins(aa) , signif_logs(aa) ] , symbols( mm , : ) , 'Color' , my_colormap( mm , : ) ) ;
            hold( h_ax , 'on' ) ;
            
        end
        
        [ h , p ] = ttest( signif_lins , signif_logs )  
        
        difference_vecs{ mm , 1 } = signif_lins - signif_logs ; 
        
    end    

    xticks( h_ax , [ 1 : size( brain_regions , 2 ) ].*2 + 0.5 ) ;
    xticklabels( h_ax, brain_regions ) ;
    ylabel( h_ax , '# Number of stimuli' ) ;
    title( h_ax , 'Linear vs. Logarithmic FMs significant response num.' ) ; 
    h_ax.FontSize = 12 ; 
    
end

