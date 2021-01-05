function h_fig = create_FM_fig2( log_FM, lin_FM, oct_speeds, lin_speeds , n_reps , PRE_STIM_TIME ,...
               cluster_num , log_PSTH_unit , lin_PSTH_unit , FM_protocol , my_colormap , color_ind )

    % The function CREATE_FM_FIG2 creates a figure of rasters of responses 
    % to FM sweeps
    
    MSEC_IN_SEC = 1000 ;
    SCALE_FACTOR = 5 ;

    GRAY_RGB = [ 36 , 36 , 36 ] ./255 ;
    RED = [ 213 , 0 , 50 ] ./255 ;
    RED = [ RED , 0.3 ] ;
    
    stim_times = FM_protocol.time_of_stim ;
    stim_times = round( MSEC_IN_SEC.* [ stim_times , fliplr( stim_times ) ] ) ;
    
    LOG_COLOR = [0,0,0] ./ 255 ;
    LIN_COLOR = [0,0,0] ./ 255  ;

    log_PSTH_unit = log_PSTH_unit ./ (max( log_PSTH_unit , [] , 2) ) .*SCALE_FACTOR ;
    lin_PSTH_unit =  lin_PSTH_unit ./ (max( lin_PSTH_unit , [] , 2) ) .*SCALE_FACTOR ;
    
    log_PSTH_unit(isnan(log_PSTH_unit) ) = 0 ;
    lin_PSTH_unit(isnan(lin_PSTH_unit) ) = 0 ;
    
    t_PSTH = linspace( 1 , size(log_FM , 2 ) , size(log_PSTH_unit , 2 ) ) ;
    
    
    h_fig = figure() ;
    h_ax = subplot('Position' , [ 0.1 , 0.15 , 0.8 , 0.30 ] ) ;
    h_ax2 = subplot('Position' , [  0.1 , 0.6 , 0.8 , 0.30 ] ) ;
    
    new_log_mat = [ ] ;
    new_lin_mat = [ ] ;
    
    oct_speed_labels = round(100 .* [-oct_speeds' ;flipud(oct_speeds')] ) ./ 100 ;
    oct_speed_labels = num2str( oct_speed_labels ) ;
    
    lin_speed_labels = round(100 .* [-lin_speeds' ;flipud(lin_speeds')] ) ./ 100 ;
    lin_speed_labels = num2str( lin_speed_labels ) ;
        
    for ii = 1 : size( log_FM , 1 ) / n_reps
        
        new_log_mat = [ new_log_mat ; log_FM((ii-1)*n_reps + 1 : ii*n_reps , : ) ] ;
        new_log_mat = [ new_log_mat ; zeros(1, size(log_FM ,2 ) ) ] ;
 
        
        new_lin_mat = [ new_lin_mat ; lin_FM((ii-1)*n_reps + 1 : ii*n_reps , : ) ] ;
        new_lin_mat = [ new_lin_mat ; zeros(1, size(lin_FM ,2 ) ) ] ;
        
 
        if ii > 1
            
            plot( h_ax , 1: size(log_FM ,2 ) , (n_reps+1)*(ii-1).*ones( 1 , size(log_FM ,2 )) , 'LineStyle' , '--' ,'Color' , GRAY_RGB ) ;
            hold( h_ax , 'on' ) ;
            plot( h_ax2 , 1: size(lin_FM ,2 ) , (n_reps+1)*(ii-1).*ones( 1 , size(lin_FM ,2 )) , 'LineStyle' , '--' ,'Color' , GRAY_RGB ) ;
            hold( h_ax2 , 'on' ) ;
            
        end
        
    end    
     
    [I1 , J1 ] = find( new_log_mat ) ;
    plot( h_ax, J1, I1, '.' , 'MarkerSize' , 6 , 'Color' , LOG_COLOR ) ; 
    hold( h_ax, 'on' ) ;

    plot( h_ax,  [ PRE_STIM_TIME + 1 , PRE_STIM_TIME + 1 ] , [ 1 , (n_reps+1)*ii ] , '-' , 'Color' ,[ my_colormap(color_ind , : ) , 0.3 ] , 'LineWidth' , 1.5 ) ;
    title(h_ax , [ 'Log FMs, Cluster num. ' num2str(cluster_num)] , 'FontSize' , 16, 'FontName' , 'Arial') ;
    xlabel( h_ax, 'Time [msec]' , 'FontSize' , 12 , 'FontName' , 'Arial') ;
%     ylabel( h_ax, 'Sweep speed [Oct/sec]' , 'FontSize' , 16 , 'FontWeight' , 'bold' ) ;
    yticks( h_ax, (n_reps+1)/2:n_reps+1:size( new_log_mat,1) ) ;
    yticklabels( h_ax, oct_speed_labels) ;
    xticks( h_ax, PRE_STIM_TIME: 500 : size(  new_log_mat ,2 ) ) ;
    xticklabels( h_ax, strrep(cellstr(num2str( ( 0 : 500 : size(  new_log_mat, 2) )' )),' ','' )) ;
    ylim( h_ax, [0 , size( new_log_mat , 1 ) ] ) ;
    xlim( h_ax, [0 , size( new_log_mat , 2 ) ] ) ; 
    h_ax.FontSize = 10 ;
    h_ax.FontName = 'Arial' ; 
    
    
    [I2 , J2 ] = find( new_lin_mat ) ;
    plot( h_ax2, J2, I2, '.' , 'MarkerSize' , 6 ,'Color', LIN_COLOR) ;
    hold( h_ax2 , 'on' ) ; 
    plot( h_ax2,  [ PRE_STIM_TIME + 1 , PRE_STIM_TIME + 1 ] , [ 1 , (n_reps+1)*ii ], '-' , 'Color' ,[ my_colormap(color_ind , : ) , 0.3 ] , 'LineWidth' , 1.5 ) ;
    title(h_ax2 , ['Lin FMs, Cluster num. ' num2str(cluster_num) ], 'FontSize' , 16 , 'FontName' , 'Arial') ;
    xlabel( h_ax2, 'Time [msec]' , 'FontSize' , 12 , 'FontName' , 'Arial') ;
%     ylabel( h_ax2, 'Sweep speed [kHz/sec]' , 'FontSize' , 16 , 'FontWeight' , 'bold' ) ;
    yticks( h_ax2, (n_reps+1)/2:n_reps+1:size( new_lin_mat,1) ) ;
    yticklabels( h_ax2, lin_speed_labels ) ;
    xticks( h_ax2, PRE_STIM_TIME: 500 : size(  new_lin_mat ,2 ) ) ;
    xticklabels( h_ax2, strrep(cellstr(num2str( ( 0 : 500 : size(  new_lin_mat, 2) )' )),' ','' )) ;
    ylim( h_ax2, [0 , size( new_lin_mat , 1 ) ] ) ;
    xlim( h_ax2, [0 , size( new_lin_mat , 2 ) ] ) ; 
    h_ax2.FontSize = 10 ;
    h_ax2.FontName = 'Arial' ;
    
    for kk = 1 : size( log_FM , 1 ) / n_reps
    
%         plot( h_ax, PRE_STIM_TIME + [ stim_times( kk ) , stim_times( kk ) ] , [ (n_reps+1)*(kk-1) + 1 , (n_reps+1)*kk ] , '-' , 'Color' , LIGHT_GRAY , 'LineWidth' , 1.5 ) ;
%         hold( h_ax , 'on' ) ;
        rectangle(h_ax, 'Position',[ PRE_STIM_TIME+1 ,  (n_reps+1)*(kk-1) ,  stim_times( kk ) ,  (n_reps+1) ],'FaceColor', [ my_colormap(color_ind , : ) , 0.3 ] , 'EdgeColor' , 'none' ,'Curvature',0.2 ) ; 
        hold( h_ax , 'on' ) ;
        rectangle( h_ax2, 'Position',[ PRE_STIM_TIME+1 ,  (n_reps+1)*(kk-1) ,  stim_times( kk ) ,  (n_reps+1) ],'FaceColor',[ my_colormap(color_ind , : ) , 0.3 ] , 'EdgeColor' , 'none' ,'Curvature',0.2) ;
        hold( h_ax2 , 'on' ) ;
        
    end
    
    hold( h_ax, 'off' ) 
    hold( h_ax2, 'off' ) 



end

