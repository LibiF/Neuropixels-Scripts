function BF_by_date_cell = plot_BF_bymouse( data_table , rec_dates , brain_regions , by_date_colormap , protocol_struct )

    % This function plots the BF distribution of units color-coded
    % according to session they were recorded in
    
    kHZ_TO_HERTZ = 1000 ;
    freqs = protocol_struct.freqs .*  kHZ_TO_HERTZ ;
    tick_dist = logspace( log( freqs(1) )/log(10) , log( freqs(end) )/log(10) , 10 ) ;
    tick_labels = round( (tick_dist(1:end)' ./ kHZ_TO_HERTZ ) .* 100 ) ./ 100  ;
    legend_cell = cell( size(rec_dates, 1 ) , 1 ) ;
    
    n_regions = size( brain_regions ,2 ) ;
    x_loc_vec = 1 : n_regions ;
    
    h_fig = figure() ;
    h_ax = axes( 'Parent' , h_fig ) ;
    
    BF_by_date_cell = cell( n_regions , 1 ) ;
   for nn = 1 : n_regions 
   
       reg_by_date_cell = cell( 1 , size( rec_dates , 1 )  ) ;
       region_units = data_table( data_table.acronym == brain_regions{ 1 , nn } , : ) ;
                   
       for kk = 1 : size( rec_dates , 1 )
        
            date_table = region_units( region_units.rec_date == rec_dates( kk , : ) , : ) ;
            reg_BFs = date_table.BF ;
            if ~isempty( reg_BFs ) 
                
                random_loc = randn( size( reg_BFs , 1 ) , 1 ).*0.06 ; 
                scatter( h_ax , nn + random_loc , log( reg_BFs .*  kHZ_TO_HERTZ)./log(10)  , 18 , by_date_colormap( kk , : ) , 'filled' ) ;
                reg_by_date_cell{ 1 , kk } = reg_BFs ;

            end
            
            hold( h_ax , 'on' ) ; 
            legend_cell{kk,1} = [num2str(kk) ] ;

       end
        
       BF_by_date_cell{ nn , 1 } = reg_by_date_cell ; 
        
   end    
        
   title( h_ax , 'BF distribution by recording' ) ;
   xticks( h_ax , 1 : n_regions ) ;
   xticklabels( h_ax , brain_regions ) ;
   yticks( h_ax , log(tick_dist)./log(10) ) ;
   yticklabels( h_ax ,  tick_labels ) ;
   ylabel( h_ax , 'BF [kHz]' , 'FontSize' , 14 ) ;
   h_ax.FontSize = 12 ;
   legend( h_ax, cellstr(legend_cell) ,'Location' , 'northeastoutside') ;

end

