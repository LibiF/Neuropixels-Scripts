function [all_dates_corr_cell , total_corr_cell , all_dates_shuffled_dist , total_corr_shuffled_cell ] =  calc_noise_correlations2( data_table , brain_regions ,...
               fluct_mat_column , bins_to_consider , iters_to_shuffle , is_for_excited )

    % The function calculates noise correlations for simultaneously
    % recorded units within the same brain region and across regions. The
    % noise correlations are calculated based on fluctuation-vecotrs
    % showing the trial-to-trial variability in stimulus responses
    % The function also calculates a shuffled distribution based on
    % shuffling the time bins of the response for "iters_to_shuffle"
    % iterations.
           
    record_dates = unique( data_table.rec_date ) ;
    inter_reg_combinations = nchoosek( 1: size( brain_regions , 2 ) , 2 ) ;
    
    all_dates_corr_cell = cell( size( brain_regions , 2 ) , size( brain_regions , 2 ) , size( record_dates , 1 ) ) ; 
    all_dates_shuffled_dist = cell( size( brain_regions , 2 ) , size( brain_regions , 2 ) , size( record_dates , 1 ) ) ;
    
    total_corr_cell = cell( size( brain_regions , 2 ) , size( brain_regions , 2 ) ) ;
    total_corr_shuffled_cell = cell( size( brain_regions , 2 ) , size( brain_regions , 2 ) ) ;
    
    for kk = 1 : size( record_dates , 1 )


        date_table = data_table( data_table.rec_date == record_dates(kk,:) , : ) ; 
        if is_for_excited
            
            date_table = date_table( not( cellfun( @(x) isempty( x ) , date_table.FRA) ) , : ) ;
        
        end
        for mm = 1 : size( brain_regions, 2 )
            
            reg1_table = date_table( date_table.acronym == brain_regions{1 , mm } , : ) ;
            reg1_fluct_mat = [] ;
            for nn = 1: size( reg1_table , 1 ) 
            
                unit_fluct = reg1_table{ nn , fluct_mat_column }{1,1}( : , bins_to_consider )  ;
                unit_fluct = unit_fluct' ;
                reg1_fluct_mat = [ reg1_fluct_mat , unit_fluct( : ) ] ; 
            
            end
            
            if ~isempty( reg1_fluct_mat )
            
                all_dates_corr_cell{ mm , mm , kk } = corrcoef( reg1_fluct_mat ) ; 
            
                shuffled_dist = [] ;
                for iter = 1 : iters_to_shuffle 
                
                    [~, reg1_fluct_randperm] = sort(rand(size( reg1_fluct_mat ,1 ) , size( reg1_fluct_mat , 2 ) )  , 1 )  ;
                    shuf_reg1_fluct = [] ;
                    for bb = 1 : size( reg1_fluct_mat , 2 ) 

                        shuf_reg1_fluct( : , bb ) = reg1_fluct_mat( reg1_fluct_randperm( : , bb ) , bb ) ; 

                    end
                    
                    shuff_noise_corrs = corrcoef( shuf_reg1_fluct ) ; 
                    shuff_noise_corrs = shuff_noise_corrs( logical( triu( ones( size(shuff_noise_corrs ) ) , 1 ) ) ) ;
                    shuffled_dist = [ shuffled_dist ; shuff_noise_corrs] ;
                    
                end
                 
                all_dates_shuffled_dist{ mm , mm , kk } = shuffled_dist ; 
                
            end
            
            comparisons = find( inter_reg_combinations(:,1) == mm )  ;
            
            if ~isempty( comparisons ) 
            
                for ll = 1: length( comparisons ) 
                
                    compared_reg = comparisons( ll ) ;
                    compared_reg = inter_reg_combinations( compared_reg , 2 ) ; 
                    reg2_table = date_table( date_table.acronym == brain_regions{1 , compared_reg } , : ) ;
                    reg2_fluct_mat = [] ;

                    for ii = 1: size( reg2_table , 1 ) 
            
                        unit_fluct2 = reg2_table{ ii , fluct_mat_column }{1,1}( : , bins_to_consider )  ;
                        unit_fluct2 = unit_fluct2' ;
                        reg2_fluct_mat = [ reg2_fluct_mat , unit_fluct2( : ) ] ; 
            
                    end
                    
                    if ~isempty( reg2_fluct_mat)
                   
                       all_dates_corr_cell{ mm , compared_reg , kk } = corr( reg1_fluct_mat , reg2_fluct_mat ) ; 
                       
                       shuffled_crossreg_dist = [] ;
                        for iter = 1 : iters_to_shuffle 

                            [~, reg2_fluct_randperm] = sort(rand(size( reg2_fluct_mat ,1 ) , size( reg2_fluct_mat , 2 ) )  , 1 )  ;
                            shuf_reg2_fluct = [] ;
                            for bb = 1 : size( reg2_fluct_mat , 2 ) 

                                shuf_reg2_fluct( : , bb ) = reg2_fluct_mat( reg2_fluct_randperm( : , bb ) , bb ) ; 

                            end

                            shuff_cross_reg_noise_corrs = corr( shuf_reg1_fluct, shuf_reg2_fluct ) ; 
                            shuff_cross_reg_noise_corrs = shuff_cross_reg_noise_corrs( : ) ;
                            shuffled_crossreg_dist = [  shuffled_crossreg_dist ; shuff_cross_reg_noise_corrs] ;

                        end

                        all_dates_shuffled_dist{ mm , compared_reg , kk } =shuffled_crossreg_dist ; 
               
                   end
                   
                end
                
            end    
                
        end

    end
    
    for cc = 1 : size( brain_regions , 2 )
        
        reg_corr_mats = all_dates_corr_cell( cc ,cc , : ) ;
        reg_corr_mats = squeeze( reg_corr_mats ) ; 
        reg_corr_indices = cellfun( @(x) triu(true(size(x)), 1 ) , reg_corr_mats , 'UniformOutput' , false ) ;
        all_corrs = [] ;
        for ff = 1 : size( reg_corr_mats , 1 )
            
            mat = reg_corr_mats{ff,1} ;
            indices = reg_corr_indices{ff,1} ;
            mat_indices = mat( indices ) ;
            all_corrs = [ all_corrs ; mat_indices ] ;  
        end
        total_corr_cell{ cc , cc } = all_corrs ;

%         reg_corr_mats = cellfun( @(x) x( cellfun(@(x) cell2mat(x) ,reg_corr_indices , 'UniformOutput' , false ) ) , reg_corr_mats , 'UniformOutput' , false ) ;
%         reg_corr_mats = cellfun( @(x) x(:) , reg_corr_mats  , 'UniformOutput' , false ) ;
        
        
%         reg_corr_mats = cellfun( @(x) triu( x, 1 )  , reg_corr_mats , 'UniformOutput' , false ) ;
%         reg_corr_mats = cellfun( @(x) x( x ~= 0) , reg_corr_mats , 'UniformOutput' , false ) ;
%         reg_corr_mats = cellfun( @(x) x(:) , reg_corr_mats  , 'UniformOutput' , false ) ;
%         total_corr_cell{ cc , cc } = cell2mat( reg_corr_mats ) ;
        
        reg_shuff_corr_mats = all_dates_shuffled_dist( cc ,cc , : ) ;
        reg_shuff_corr_mats = squeeze( reg_shuff_corr_mats ) ; 
        reg_shuff_corr_mats = cell2mat( reg_shuff_corr_mats ) ;
        total_corr_shuffled_cell{ cc , cc } = reg_shuff_corr_mats ;
         
         for dd = cc+1 : size( brain_regions , 2 )

            cross_reg_mat = all_dates_corr_cell( cc , dd , : ) ;
            cross_reg_mat = squeeze( cross_reg_mat ) ;
            cross_reg_mat = cellfun( @(x) x(:) , cross_reg_mat , 'UniformOutput' , false ) ;  
            total_corr_cell{ cc , dd } = cell2mat( cross_reg_mat ) ; 
            
            cross_reg_shuffled = all_dates_shuffled_dist( cc , dd , : ) ;
            cross_reg_shuffled = squeeze( cross_reg_shuffled ) ;
            total_corr_shuffled_cell{ cc , dd } = cell2mat( cross_reg_shuffled ) ; 
            
         end
         
    end

end
