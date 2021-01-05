function [ within_reg_corrs , all_dates_corr_cell ] = calculate_sig_corrs3( data_table, brain_regions , iters_for_shuffle )

    % This function calculates signal correlations between units
    % recorded within each brain regions.
    % Signal correlations are defined as the correlations between FRAs.
    % In addition, shuffled correlation distributions are calculated based
    % on #"iters_to_shuffle" iterations.
    % Correlations are calculated for all unit pairs and only for
    % simultaneously recorded unit pairs
    % The function output is:
    % within_reg_corrs - a cell with n rows and 4 columns:
    % n - corresponds to the number of cross region combinations
    % Columns:
    % 1 - correlations between all units
    % 2 - correlations between simultaneously recorded units
    % 3 - shuffled distribution for all units
    % 4 - shuffled distribution for simultaneously recorded pairs
    % The output variable all_dates_corr_cell is a 3D cells containing the
    % cross region pairwise correlations for all recording sessions
    % (meaning, only for simultaneously recorded units)
    
    record_dates = unique( data_table.rec_date ) ;
    within_reg_corrs = cell( size( brain_regions , 2 ) , 4 ) ;
    all_dates_corr_cell = cell( size( brain_regions , 2 ) , size( brain_regions , 2 ) , size( record_dates , 1 ) ) ; 
    
    for kk = 1 : size( brain_regions , 2 )
        
        region_table = data_table( data_table.acronym == brain_regions{ 1 , kk } , : ) ;
%         rec_dates = unique( region_table.rec_date ) ;
        rec_dates = record_dates ;
        reg_FRAs = [ region_table.FRA , region_table.suppressive_FRA ] ;
        suppressed_units = cellfun( @(x) isempty(x), reg_FRAs(:,1) ) ;
        reg_FRAs(suppressed_units, 1 ) = reg_FRAs(suppressed_units , 2 ) ;
        reg_FRAs = reg_FRAs( 1:size( reg_FRAs , 1 ) , 1 ) ; 
        reg_FRAs = cellfun( @(x) x(:)' , reg_FRAs , 'UniformOutput' , false ) ;
        reg_FRAs = cell2mat( reg_FRAs )' ; 
        reg_sig_corr = corr( reg_FRAs ) ;
        h_fig_mat = figure();
        h_ax_mat = axes('Parent' , h_fig_mat ) ;
        imagesc( h_ax_mat , reg_sig_corr ) ;
        title( [ brain_regions{ 1 , kk } ' - Signal Correlations ' ] )
         
        reg_corrs_vec = reg_sig_corr( logical( triu( ones( size( reg_sig_corr ) ) , 1 ) ) ) ;
        within_reg_corrs{ kk , 1 } = reg_corrs_vec ;
       
        shuffled_dist = [] ;
        for aa = 1 : iters_for_shuffle
        
            [~, FRA_randperm] = sort(rand(size( reg_FRAs ,1 ) , size( reg_FRAs , 2 ) )  , 1 )  ;
            
            for bb = 1 : size( reg_FRAs , 2 ) 
                
                shuf_FRA( : , bb ) = reg_FRAs( FRA_randperm( : , bb ) , bb ) ; 
            
            end
            
            shuf_corr = corr( shuf_FRA ) ;
            reg_shuf_corrs_vec = shuf_corr( logical( triu( ones( size( shuf_corr ) ) , 1 ) ) ) ;
            shuffled_dist = [ shuffled_dist ; reg_shuf_corrs_vec ] ;
        
        end
        within_reg_corrs{ kk , 3 } = shuffled_dist ; 
        
        by_date_corrs = [] ;    
        by_date_shuffled_corrs = [] ;        
        h_fig = figure() ;
       
        for mm = 1 : size( record_dates , 1 )

            date_table = region_table( region_table.rec_date == rec_dates( mm , : ) , : ) ;
            if ~isempty( date_table  )
                
                reg_date_FRAs = [ date_table.FRA , date_table.suppressive_FRA ] ;
                date_suppressed_units = cellfun( @(x) isempty(x), reg_date_FRAs(:,1) ) ;
                reg_date_FRAs( date_suppressed_units , 1 ) = reg_date_FRAs( date_suppressed_units , 2 ) ;
                reg_date_FRAs = reg_date_FRAs( 1 : size( reg_date_FRAs ) , 1 ) ;
                
                reg_date_FRAs = cellfun( @(x) x(:)' , reg_date_FRAs , 'UniformOutput' , false ) ;
                reg_date_FRAs = cell2mat( reg_date_FRAs )' ; 
                reg_date_sig_corr = corr( reg_date_FRAs ) ;
                h_ax = subplot( 2 , ceil( size( rec_dates , 1 )/ 2 ) , mm ) ;
                imagesc( h_ax , reg_date_sig_corr ) ;
                title( h_ax, char( rec_dates( mm , : ) ) ) ;

                all_dates_corr_cell{ kk , kk , mm } = reg_date_sig_corr ; 

                reg_date_sig_corrs_vec = reg_date_sig_corr( logical( triu( ones( size( reg_date_sig_corr ) ) , 1 ) ) ) ;
                by_date_corrs = [ by_date_corrs ; reg_date_sig_corrs_vec ] ;

                for aa = 1 : iters_for_shuffle

                    [~, FRA_date_randperm] = sort(rand(size( reg_date_FRAs ,1 ) , size( reg_date_FRAs , 2 ) )  , 1 )  ;

                    for bb = 1 : size( reg_date_FRAs , 2 ) 

                        shuf_date_FRA( : , bb ) = reg_date_FRAs( FRA_date_randperm( : , bb ) , bb ) ; 

                    end
                    shuf_date_corr = corr( shuf_date_FRA ) ;
                    reg_shuf_date_corrs_vec = shuf_date_corr( logical( triu( ones( size( shuf_date_corr ) ) , 1 ) ) ) ;
                    by_date_shuffled_corrs = [ by_date_shuffled_corrs ; reg_shuf_date_corrs_vec ] ;

                end


            end   

            within_reg_corrs{ kk , 2 } = by_date_corrs ; 
            within_reg_corrs{ kk , 4 } = by_date_shuffled_corrs ; 
            
        end
    
    end
    
end

