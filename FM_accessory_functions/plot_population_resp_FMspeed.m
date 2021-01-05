function curves_cell = plot_population_resp_FMspeed( data_table ,...
         brain_regions , is_normalized, data_col , data_sub_col ,...
         my_colormap , title_text )

    % The function PLOT_POPULATION_RESP_FMSPEED receives the data table
    % after the mean response for different FM speeds/velocities have
    % already been calculated. The function receives also a data column
    % (data_col) and sub column and looks for data within these columns to
    % plot the FM speed selectivity.
    % The function also receives the parameter "is_normalized" determines
    % whether the population curves will be summed up after each unit is
    % normalized to it's maximal response or not.
    % The function uses the function "plot_areaerrorbar" by 
    % Victor Martinez-Cagigal that can be found on
    % https://www.mathworks.com/matlabcentral/fileexchange/58262-shaded-area-error-bar-plot
    
    % Define parameters and data structures
    MSEC_IN_SEC = 1000 ;    
    options.handle = figure() ;
    h_ax= axes( 'Parent' , options.handle ) ; 
    new_legend = cell( size( brain_regions , 2 )*2 , 1 ) ;
    curves_cell = cell( size( brain_regions , 2) , 1 ) ;
    
    for kk = 1 : size( brain_regions , 2 )

        new_legend{ 2*(kk-1) + 1 , 1 } = brain_regions{ 1,  kk } ;
        new_legend{ 2*(kk-1) + 2 , 1 } = 'SEM' ;
        
        wanted_table = data_table( data_table.acronym == brain_regions{ 1 , kk } , : ) ;
        curves = wanted_table{ [ 1:end ] , data_col }([1 :end] , data_sub_col ) ;
        curves = cell2mat(curves) ;
       
        if is_normalized
              
            normalized_curves = curves ./ max( curves , [] , 2 ) ;

        else
            
            normalized_curves = curves ;
            
        end   
             
        curves_cell{ kk ,1 } = normalized_curves ; 
        mean_curves = nanmean( normalized_curves ) ; 

        gca() ;
        options.color_line = my_colormap(kk,:) ;
        options.error = 'sem' ;
        options.color_area = my_colormap(kk,:) ; 
        options.alpha = 0.5 ; 
        options.line_width  = 2.5 ;
        options.x_axis = 1 : length( mean_curves ) ; 
        plot_areaerrorbar( normalized_curves, options ) ; 
        hold( h_ax , 'on' ) ; 
        
    end
    
    title( h_ax, title_text ) ;
     
end

