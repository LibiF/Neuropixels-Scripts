function plot_evoked_FRs_distribution( evoked_FR_cell , brain_regions , suptitle_text , my_colormap)

    % This function plots the evoked FR distribution in regular and
    % log-scale
    
    h_fig = figure(); 
    n_brain_regions = length( brain_regions ) ;
    for kk = 1 : n_brain_regions
        
        reg_spont_FR = evoked_FR_cell{ kk, 1 }  ;        
        
        h_ax = subplot( 2 , n_brain_regions  , kk ) ; 
        histogram( h_ax , reg_spont_FR , 'Normalization' , 'probability' , 'FaceColor' , my_colormap(kk, : ) , 'BinWidth' , 5) ; 
        xlim( h_ax , [ 0 , 140 ] ) ;
        ylim( h_ax , [ 0 , 0.35 ] ) ;
        xlabel( h_ax , [ 'Evoked FR [Hz]' ] ) ;
        title( h_ax , brain_regions{ 1 , kk } ) ;
        h_ax2 = subplot( 2 , n_brain_regions  , kk+n_brain_regions ) ;
        histogram( h_ax2 , log( reg_spont_FR(:) ) , 'Normalization' , 'probability' , 'FaceColor' , my_colormap(kk, : ) ) ; 
        xlim( h_ax2 , [0 , 6 ] ) ;
        xlabel( h_ax2 , 'log (Evoked FR)' )
        ylim( h_ax2 , [ 0 , 0.4 ] ) ;
%         suptitle( suptitle_text ) ;
        
    end       
    


end

