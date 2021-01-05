function [ ttest_scores , time_to_signif ] = find_time_to_baseline_ttest( data_cell ,...
           brain_regions , PRE_STIM_MSEC , pval_thresh , time_delay) 
    
    % The function FIND_TIME_TO_BASELINE_TTEST receives a data cell
    % contating population responses to different stimuli within different
    % regions.
    % For each region, it statistically compares each time point in the
    % population vector to the mean population activity during baseline
    % (defined as time 0:PRE_STIM_MSEC).
    % After the statistics are computed, the function detects the first 
    % time points in which the population activity deviates from baseline 
    % and the time point in which the population activity returns to
    % baseline (both values returned in "time_to_signif").
    % Note that first time point to return to baseline is detected
    % following a specified temporal delay (time_delay) to avoid premature
    % detections (for ex. when there is post stimulus suppression which is
    % part of the population response).
    % Statiscs are corrected for multiple comparisons by looking for a
    % p-value corrected for the number of compared time bins (pval_thresh)


    ttest_scores = zeros( size( data_cell , 2 ) , size( data_cell(1,1) , 2 ) ) ; 
    times_to_signif = zeros ( size( data_cell , 2 )  , 2 ) ;
    
    for cc = 1 : size( data_cell , 2 ) 
        
        pop = data_cell( 1 : end , cc ) ;
        pop = cell2mat( pop ) ;

        baseline_vec = mean( pop( : , 1: PRE_STIM_MSEC ) , 2 ) ;
        
        for tt = 1 : size( pop , 2 ) 
        
            [ ~ , ttest_scores( cc , tt ) ] = ttest2( pop( : , tt ) , baseline_vec ) ;
                                   
        end    

        time_from_basline = find( ttest_scores(cc,:) < pval_thresh , 1 ) ;
        time_to_baseline = ttest_scores(cc,:) > pval_thresh ;
        [ ~ , time_to_baseline ] = max( time_to_baseline(:,time_delay+PRE_STIM_MSEC:end) , [] , 2 );
        time_to_baseline  = time_to_baseline + time_delay + PRE_STIM_MSEC ; 
        time_to_signif( cc , 1 ) = time_from_basline ;
        time_to_signif( cc, 2 ) = time_to_baseline ; 
        
        h_fig = figure() ;
        h_ax = axes( 'Parent' , h_fig ) ;
        imagesc( h_ax , ttest_scores( cc , : ) ) ;
        caxis( h_ax , [ 1e-10 , 0.1] ) ;
        decimals = (-10:1:-1) ;
        ticks = 10.^decimals ;
        set( h_ax, 'ColorScale' ,'log' ) ;
        cb = colorbar  ;
        title( h_ax , [ brain_regions(cc)] ) ;
        h_ax.XTickLabel = {[]} ;
        h_ax.YTick = [] ;
        colormap( bone ) 
        hold( h_ax, 'on' )
        line( h_ax , [ time_from_basline , time_from_basline ] , [ 0.5 , 1.5 ] , 'Color' ,'r' , 'LineWidth' , 2 ) ;
        hold( h_ax ,'on' )
        line( h_ax , [time_to_baseline , time_to_baseline ] , [ 0.5 , 1.5 ] , 'Color' ,'r' , 'LineWidth' , 2 ) ; 
        
        
    end
    
    time_to_signif = time_to_signif - PRE_STIM_MSEC ;
    
end

