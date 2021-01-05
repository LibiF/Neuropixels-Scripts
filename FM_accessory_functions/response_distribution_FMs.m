function [ significant_response_mat_log , significant_response_nums_log ,...
           significant_response_mat_lin , significant_response_nums_lin ] = response_distribution_FMs( data_table ,...
           brain_regions , if_print , title_text )

    % The function RESPONSE_DISTRIBUTION_FMS summarizes the FM response by
    % calculating the fraction of units excited, suppressed or both for
    % linear and logarithmic FMs. The function outputs two types of
    % variables (for each of the FM stimulus types):
    % significant_response_mat - a matrix size r x 3 where r is the number
    % of brain regions and the three columns correspond to the three
    % response types (excitation, suppression, both). The values in the
    % matrix are the fraction of units in the region showing this type of
    % response.
    % significant_response_nums - a r x 4 matrix where r is the number of
    % brain regions. Values in each row x column correspond to the number
    % of units showing this type of response (excitatory, suppressive,
    % both) and the 4th column corresponds to the total number of units.
    % The function also generates a bar graphs summarizing the statistics.
    % If the input parameter "if_print" is set to true the summary
    % statistics will be printed to the MATLAB cmd window
       
    legend_text = { 'Excited' , 'Suppressed' , 'Both' } ;
    significant_response_mat_log = zeros( size( brain_regions , 2 ) , 3 ) ;
    significant_response_nums_log = zeros( size( brain_regions , 2 ) , 4 ) ;
    
    significant_response_mat_lin = zeros( size( brain_regions , 2 ) , 3 ) ;
    significant_response_nums_lin = zeros( size( brain_regions , 2 ) , 4 ) ;
    
    for kk = 1 : length( brain_regions)
        
        reg_units = data_table( data_table.acronym == brain_regions{ 1 , kk } , : ) ;
        
        % log
        excited_units_log = reg_units.is_excited_log ;
        inhibited_units_log = reg_units.is_inhibited_log ;
        excited_inhibited_log = nansum( excited_units_log + inhibited_units_log , 2 ) ;
        excited_fract_log = nansum( excited_units_log ) ./ size( excited_units_log ,1 ) ;
        inhibited_fract_log = nansum( inhibited_units_log ) ./ size( inhibited_units_log ,1 ) ;
        both_fract_log = sum( excited_inhibited_log == 2 ) ./ size( excited_inhibited_log , 1 ) ;
        
        significant_response_mat_log( kk , 1 ) = excited_fract_log ;
        significant_response_mat_log( kk , 2 ) = inhibited_fract_log ;
        significant_response_mat_log( kk , 3 ) = both_fract_log ;
        
        significant_response_nums_log( kk , 1 ) = nansum( excited_units_log ) ;
        significant_response_nums_log( kk , 2 ) = nansum( inhibited_units_log ) ;
        significant_response_nums_log( kk ,3 ) = sum( excited_inhibited_log == 2 ) ;
        significant_response_nums_log( kk ,4 ) = size( reg_units ,1 ) ;
        
        %lin
        excited_units_lin = reg_units.is_excited_lin ;
        inhibited_units_lin = reg_units.is_inhibited_lin ;
        excited_inhibited_lin = nansum( excited_units_lin + inhibited_units_lin , 2 ) ;
        excited_fract_lin = nansum( excited_units_lin ) ./ size( excited_units_lin ,1 ) ;
        inhibited_fract_lin = nansum( inhibited_units_lin ) ./ size( inhibited_units_lin ,1 ) ;
        both_fract_lin = sum( excited_inhibited_lin == 2 ) ./ size( excited_inhibited_lin , 1 ) ;
        
        significant_response_mat_lin( kk , 1 ) = excited_fract_lin ;
        significant_response_mat_lin( kk , 2 ) = inhibited_fract_lin ;
        significant_response_mat_lin( kk , 3 ) = both_fract_lin ;
        
        significant_response_nums_lin( kk , 1 ) = nansum( excited_units_lin ) ;
        significant_response_nums_lin( kk , 2 ) = nansum( inhibited_units_lin ) ;
        significant_response_nums_lin( kk ,3 ) = sum( excited_inhibited_lin == 2 ) ;
        significant_response_nums_lin( kk ,4 ) = size( reg_units ,1 ) ;
        
        if if_print
            
            disp( ['Region - ' brain_regions{1,kk}  ] )
            disp( ['Fraction of excited units (log) ' num2str(nansum( excited_units_log )) ' , '  num2str(excited_fract_log) ] ) ;
            disp( ['Fraction of inhibited units (log)' num2str(nansum(inhibited_units_log )) ' , '   num2str(inhibited_fract_log) ] ) ;
            disp( ['Fraction of excited and inhibited units (log)'  num2str(nansum( excited_inhibited_log == 2 ) ) ' , '  num2str( both_fract_log) ] ) ;
            
            disp( ['Region - ' brain_regions{1,kk}  ] )
            disp( ['Fraction of excited units (lin) ' num2str(nansum( excited_units_lin )) ' , '  num2str(excited_fract_lin) ] ) ;
            disp( ['Fraction of inhibited units (lin)' num2str(nansum(inhibited_units_lin )) ' , '   num2str(inhibited_fract_lin) ] ) ;
            disp( ['Fraction of excited and inhibited units (lin)'  num2str(nansum( excited_inhibited_lin == 2 ) ) ' , '  num2str( both_fract_lin) ] ) ;
            
        end    
        
    end   

    h_fig = figure() ;
    h_ax = axes( 'Parent' , h_fig ) ;
    X_vals = categorical( brain_regions ) ;
    bar( h_ax , X_vals , significant_response_mat_log ) ; 
    h_ax.FontSize = 16 ;
    h_ax.FontWeight = 'bold' ; 
    legend( h_ax, legend_text , 'Orientation' , 'horizontal' ,...
            'Location' , 'north' , 'FontSize' , 8 ) ;
    ylim( h_ax, [ 0 , 1 ] ) ;    
    title( h_ax,  [ title_text ' - log' ]  , 'FontSize' , 18) ;
    
    h_fig2 = figure() ;
    h_ax2 = axes( 'Parent' , h_fig2 ) ;
    X_vals = categorical( brain_regions ) ;
    bar( h_ax2 , X_vals , significant_response_mat_lin ) ; 
    h_ax2.FontSize = 16 ;
    h_ax2.FontWeight = 'bold' ; 
    legend( h_ax2, legend_text , 'Orientation' , 'horizontal' ,...
            'Location' , 'north' , 'FontSize' , 8 ) ;
    ylim( h_ax2, [ 0 , 1 ] ) ;    
    title( h_ax2,  [ title_text ' - lin' ] , 'FontSize' , 18) ;
  
end

