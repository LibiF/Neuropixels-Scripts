function [ significant_response_mat , significant_response_nums ] = response_distribution( data_table ,...
           brain_regions , IF_PRINT , title_text )

    % The function RESPONSE_DISTRIBUTION calculates statistics of the pure
    % tone response in the different brain regions and plots the
    % corresponding bar graph

    legend_text = { 'Excited' , 'Suppressed' , 'Both' } ;
    significant_response_mat = zeros( size( brain_regions , 2 ) , 3 ) ;
    significant_response_nums = zeros( size( brain_regions , 2 ) , 4 ) ;
    for kk = 1 : length( brain_regions)
        
        reg_units = data_table( data_table.acronym == brain_regions{ 1 , kk } , : ) ;
        excited_units = reg_units.is_excited ;
        inhibited_units = reg_units.is_inhibited ;
        excited_inhibited = nansum( excited_units + inhibited_units , 2 ) ;
        excited_fract = nansum( excited_units ) ./ size( excited_units ,1 ) ;
        inhibited_fract = nansum( inhibited_units ) ./ size( inhibited_units ,1 ) ;
        both_fract = sum( excited_inhibited == 2 ) ./ size( excited_inhibited , 1 ) ;
        
        significant_response_mat( kk , 1 ) = excited_fract ;
        significant_response_mat( kk , 2 ) = inhibited_fract ;
        significant_response_mat( kk , 3 ) = both_fract ;
        
        significant_response_nums( kk , 1 ) = nansum( excited_units ) ;
        significant_response_nums( kk , 2 ) = nansum( inhibited_units ) ;
        significant_response_nums( kk ,3 ) = sum( excited_inhibited == 2 ) ;
        significant_response_nums( kk ,4 ) = size( reg_units ,1 ) ;
        
        if IF_PRINT
            
            disp( ['Region - ' brain_regions{1,kk}  ] )
            disp( ['Fraction of excited units ' num2str(nansum( excited_units )) ' , '  num2str(excited_fract) ] ) ;
            disp( ['Fraction of inhibited units ' num2str(nansum(inhibited_units )) ' , '   num2str(inhibited_fract) ] ) ;
            disp( ['Fraction of excited and inhibited units '  num2str(nansum( excited_inhibited == 2 ) ) ' , '  num2str( both_fract) ] ) ;
            
        end    
        
    end   

    h_fig = figure() ;
    h_ax = axes( 'Parent' , h_fig ) ;
    X_vals = categorical( brain_regions ) ;
    bar( h_ax , X_vals , significant_response_mat ) ; 
    h_ax.FontSize = 16 ;
    h_ax.FontWeight = 'bold' ; 
    legend( h_ax, legend_text , 'Orientation' , 'horizontal' ,...
            'Location' , 'north' , 'FontSize' , 8 ) ;
    ylim( h_ax, [ 0 , 1 ] ) ;    
    title( h_ax,  title_text , 'FontSize' , 18) ;
  
    
end

