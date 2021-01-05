function spont_FR_cell = plot_spont_FR( data_table , brain_regions , IF_PRINT , my_colormap )

    % The function PLOT_SPONT_FR plots the spontaneous firing rate 
    % distributions of all brain region in regular and log scale
    % This function should be used after spontaneous firing rates were
    % already calculated and assigned to the table column "spont_FR"
    
    % The boolean variable IF_PRINT is used to indicate whether the user
    % wants the summary statistics to be printed to the MATLAB cmd window
    
    spont_FR_cell = cell( size(brain_regions , 2 ) ,1 ) ;
    n_brain_regions = size(brain_regions , 2 ) ; 
    h_fig = figure(); 

    for kk = 1 : length( brain_regions)
        
        reg_spont_FR = data_table.spont_FR( data_table.acronym == brain_regions{ 1,kk } , : ) ;
        spont_FR_cell{ kk ,1 } = reg_spont_FR ;
        
        if IF_PRINT
            
            disp( ['Region - ' brain_regions{1,kk}  ] )
            disp( ['Mean spont FR '  num2str(mean( spont_FR_cell{ kk , 1 } ) ) ] ) ;
            disp( ['Std spont FR '  num2str(std( spont_FR_cell{ kk , 1 } ) ) ] ) ;
            disp( ['Median spont FR '  num2str(median( spont_FR_cell{ kk , 1 } ) ) ] ) ;
            disp( [ 'Quantiles (25%,75%) spont FR ' num2str( quantile(spont_FR_cell{ kk , 1 },[ 0.25 , 0.75 ] ) ) ] ) ;
            disp( [ 'Num of units' num2str( size( spont_FR_cell{ kk ,1 } , 1 ) ) ] ) ;
            
        end    
        
        h_ax = subplot( 2 , n_brain_regions  , kk ) ; 
        histogram( h_ax , reg_spont_FR , 'Normalization' , 'probability' , 'FaceColor' , my_colormap(kk, : ) , 'BinWidth' , 5) ; 
        xlim( h_ax , [0 , 50 ] ) ;
        ylim( h_ax , [0 , 1 ] ) ;
        xlabel( h_ax , 'Spontaneous FR' ) ;
        title( h_ax , brain_regions{ 1 , kk } ) ;
        h_ax2 = subplot( 2 , n_brain_regions  , kk+n_brain_regions ) ;
        histogram( h_ax2 , log( reg_spont_FR(:) ) , 'Normalization' , 'probability' , 'FaceColor' , my_colormap(kk, : ) ) ; 
        xlim( h_ax2 , [- 5 , 5 ] ) ;
        xlabel( h_ax2 , 'log(Spontaneous FR)' )
        
    end       
    
end

