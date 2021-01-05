function bandwidth_cell = calculate_plot_bandwidth( data_table , brain_regions ,...
        FP_num , protocol_struct , my_colormap , if_print , display_mode ) 

    % The function CALCULATE_PLOT_BANDWIDTH takes the data_table and
    % calculates based on each unit's FRA significance the number of
    % stimuli for which it shows a significant auditory response.
    % The number of expected false positives (FP_num) is subtracted from
    % this quantity to obtain the true number of "significant stimuli".
    % If the number of significant stimuli is smaller than what is expected
    % from FPs, it is rounded to 0.
    % IF_PRINT - controls the display of statistics
    % DISPLAY_MODE - controls the presentation of response histograms

    total_stimuli = length( protocol_struct.atten ) .* length( protocol_struct.freqs ) ;

    bandwidth_cell = cell( size(brain_regions , 2 ) ,1 ) ;
    n_brain_regions = size(brain_regions , 2 ) ; 
    h_fig = figure(); 

    for kk = 1 : length( brain_regions)
        
        reg_table = data_table( data_table.acronym == brain_regions{ 1,kk } , : ) ;
        FRA_significance = reg_table.FRA_signif ;
        FRA_signif_indices = cellfun( @(x) isempty(x) , FRA_significance ) ;
        FRA_significance = FRA_significance( not( FRA_signif_indices ) ) ;
        num_signif_responses = cellfun( @(x) sum( sum( x ) ) , FRA_significance ) ;
        num_signif_responses = num_signif_responses - FP_num ;
        num_signif_responses( num_signif_responses < 0 ) = 0 ;
        
        h_ax = subplot( 1 , n_brain_regions  , kk ) ; 
        if display_mode == 1 
            
            num_signif_responses = num_signif_responses ./ total_stimuli ;
            histogram( h_ax , num_signif_responses , 'Normalization' , 'probability' , 'FaceColor' , my_colormap(kk, : ) , 'BinWidth' , 0.1) ; 
            xlim( h_ax , [0 , 1 ] ) ;
            xlabel( h_ax , 'Frac. signif. responses' ) ;
            
        else
            
            histogram( h_ax , num_signif_responses , 'Normalization' , 'probability' , 'FaceColor' , my_colormap(kk, : ) , 'BinWidth' , 5 ) ; 
            xlim( h_ax , [0 , total_stimuli - 20 ] ) ;
            xlabel( h_ax , 'Num. signif. responses' ) ;

        end
        
        bandwidth_cell{ kk , 1 } = num_signif_responses ;
        ylim( h_ax , [0 , 1 ] ) ;
        title( h_ax , brain_regions{ 1 , kk } ) ;
        
        if if_print
            
            disp( ['Region - ' brain_regions{1,kk}  ] )
            disp( ['Mean significant fraction '  num2str(mean( num_signif_responses ) ) ] ) ;
            disp( ['Median significant fraction '  num2str(median( num_signif_responses ) ) ] ) ;

        end    
        
    end       

end

