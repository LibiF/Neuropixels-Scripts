function fract_sum_mat = fraction_exc_inh_allFMs( data_table , brain_regions )

    % The function FRACTION_EXC_INH_ALLFMS receives the data table and
    % a cell of brain regions to analyze and outputs a n x 5 matrix with n-
    % the number of brain regions and 5 columns correspond to number of
    % units:
    % 1 - Excited
    % 2 - Suppressed
    % 3 - Both
    % 4 - Non responsive
    % 5 - Total num.
    % The summary is over all FMs and FM types meaning that a unit showing
    % an excitatory to at least one FM stimulus (linear or logarithmic)
    % will be considered excited and likewise for suppressive/both.
    % Units which are non-responsive did not show any significant response
    % to any of the total FM stimuli (20 in my protocol)
    
    RED = [	174, 96, 88 ] ./ 256 ; 
    BLUE = [ 1,85,151] ./ 256 ; 
    GRAY =  [211,211,211]./256 ;
    
    fract_sum_mat = zeros( size( brain_regions , 2 ) , 5 ) ;

    for kk = 1 : size( brain_regions , 2 )

        reg_table = data_table( data_table.acronym == brain_regions{ 1 , kk } , : ) ;
        exc_signif = [cell2mat( reg_table.lin_sig_mat) , cell2mat( reg_table.log_sig_mat ) ] ;
        exc_signif = any( exc_signif , 2 ) ; 
        
        inh_signif = [cell2mat( reg_table.lin_sig_inh_mat) , cell2mat( reg_table.log_sig_inh_mat ) ] ;
        inh_signif = any( inh_signif , 2 ) ;
        
        sum_exc = sum( exc_signif ) ;
        all_inh = sum( inh_signif ) ;
        only_inh = sum( not(exc_signif) & inh_signif == 1 ) ;
        exc_inh = all_inh - only_inh ;
        only_exc = sum_exc-exc_inh ;

        not_responsive = size( reg_table , 1 ) - sum_exc - only_inh ;
        fract_sum_mat( kk , 1 ) = only_exc ;
        fract_sum_mat( kk , 2) = exc_inh ; 
        fract_sum_mat( kk , 3 ) = only_inh ;
        fract_sum_mat( kk , 4 ) = not_responsive ; 
        fract_sum_mat( kk , 5 ) = sum( fract_sum_mat(  kk , 1:4 ) ) ;
        
    end

end

