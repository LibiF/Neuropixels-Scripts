function [lin_pop_sparseness , log_pop_sparseness] = calculate_FM_popsparseness(data_table ,...
         brain_regions , protocol_struct)


    MSEC_IN_SEC = 1000 ;
    n_stims_per_dir = length( protocol_struct.oct_speeds) ;
    stim_times = [ protocol_struct.time_of_stim , fliplr( protocol_struct.time_of_stim ) ] .* MSEC_IN_SEC' ; 
    lin_slopes = ( protocol_struct.upper_freq - protocol_struct.bottom_freq ) ./ stim_times ;
    lin_slopes( 1 : n_stims_per_dir ) = - lin_slopes( 1: n_stims_per_dir ) ;
    lin_slopes = round( lin_slopes .* 100 ) ./ 100 ; 

    log_slopes =[ - protocol_struct.oct_speeds , fliplr( protocol_struct.oct_speeds ) ] ; 

    lin_pop_spaseness = cell( size( brain_regions , 2 ) , 1 ) ;
    log_pop_spaseness = cell( size( brain_regions , 2 ) , 1 ) ;

    for mm = 1 : size( brain_regions , 2 ) 

        reg_units = data_table( data_table.acronym == brain_regions{ 1 , mm } , : ) ;
        lin_sig_mat = cell2mat( reg_units.lin_sig_mat ) ;
        log_sig_mat = cell2mat( reg_units.log_sig_mat ) ;

        lin_sig_fract = sum( lin_sig_mat , 1 ) ./ size ( lin_sig_mat , 1 ) ;
        log_sig_fract = sum( log_sig_mat , 1 ) ./ size ( log_sig_mat , 1 ) ;

        lin_pop_spaseness{ mm ,1 } = lin_sig_fract ;
        log_pop_spaseness{ mm ,1 } = log_sig_fract ;

    end

    lin_pop_sparseness = cell2mat( FM_rez.lin_pop_spaseness ) ;
    h_fig = figure() ;
    h_ax = axes('Parent' , h_fig ) ;
    h_bar = bar( h_ax , categorical( lin_slopes ), lin_pop_sparseness' , 'FaceColor' , 'flat' ) ;
    for mm = 1 : size( brain_regions , 2 )

        h_bar(mm).CData = my_colormap( mm , : ) ;

    end

    ylabel( h_ax , 'Fraction of population responsive' ) ;
    xlabel( h_ax , 'Slopes [kHz/sec]' ) ;
    title( h_ax , 'Population sparseness - linear FMs' ) ;

    log_pop_sparseness = cell2mat( FM_rez.log_pop_spaseness ) ;
    h_fig2 = figure() ;
    h_ax2 = axes('Parent' , h_fig2 ) ;
    h_bar2 = bar( h_ax2 , categorical( log_slopes ), log_pop_sparseness' , 'FaceColor' , 'flat' ) ;
    for mm = 1 : size( brain_regions , 2 )

        h_bar2(mm).CData = my_colormap( mm , : ) ;

    end
    ylabel( h_ax2 , 'Fraction of population responsive' ) ;
    xlabel( h_ax2 , 'Slopes [Oct/sec]' ) ;
    title( h_ax2 , 'Population sparseness - logarithmic FMs' ) ;

end

