function [ across_reg_corrs , all_dates_corr_cell ] = calculate_crossreg_sig_corrs3( data_table , brain_regions , iters_for_shuffle )

    % This function calculates signal correlations between units
    % recorded from different brain regions.
    % Signal correlations are defined as the correlations between FRAs.
    % In addition, shuffled correlation distributions are calculated based
    % on #"iters_to_shuffle" iterations.
    % Correlations are calculated for all unit pairs and only for
    % simultaneously recorded unit pairs
    % The function output is:
    % across_reg_corrs - a cell with n rows and 4 columns:
    % n - corresponds to the number of cross region combinations
    % Columns:
    % 1 - correlations between all units
    % 2 - correlations between simultaneously recorded units
    % 3 - shuffled distribution for all units
    % 4 - shuffled distribution for simultaneously recorded pairs
    % The output variable all_dates_corr_cell is a 3D cells containing the
    % cross region pairwise correlations for all recording sessions
    % (meaning, only for simultaneously recorded units)
    
    reg_combos = nchoosek( 1 : size( brain_regions , 2 ) , 2 ) ;
    across_reg_corrs = cell( size( reg_combos , 1 ) , 4 ) ;
    rec_dates = unique( data_table.rec_date ) ;
    all_dates_corr_cell = cell( size( brain_regions , 2 ) , size( brain_regions , 2 ) , size( rec_dates , 1 ) ) ; 

    
    for kk = 1 : size( reg_combos , 1 )
        
        combo = reg_combos( kk , : ) ; 
        reg1_table = data_table( data_table.acronym == brain_regions{ 1 , combo(1) } , : ) ;
%         reg1_table = reg1_table( not( cellfun( @(x) isempty( x ) , reg1_table.FRA) ) , : ) ;
        reg1_FRAs = [ reg1_table.FRA , reg1_table.suppressive_FRA ] ;
        reg1_suppressed_units = cellfun( @(x) isempty(x),  reg1_FRAs(:,1) ) ;
        reg1_FRAs(  reg1_suppressed_units , 1 ) = reg1_FRAs(  reg1_suppressed_units , 2 ) ;
        reg1_FRAs = reg1_FRAs( 1:size( reg1_FRAs , 1 ) , 1 ) ;
        reg1_FRAs = cellfun( @(x) x(:)' , reg1_FRAs , 'UniformOutput' , false ) ;
        reg1_FRAs = cell2mat( reg1_FRAs )' ; 
        
        reg2_table = data_table( data_table.acronym == brain_regions{ 1 , combo(2) } , : ) ;
%         reg2_table = reg2_table( not( cellfun( @(x) isempty( x ) , reg2_table.FRA) ) , : ) ;
%         rec_dates = unique( reg2_table.rec_date ) ;
        
        reg2_FRAs = [ reg2_table.FRA , reg2_table.suppressive_FRA ] ;
        reg2_suppressed_units = cellfun( @(x) isempty(x),  reg2_FRAs(:,1) ) ;
        reg2_FRAs(  reg2_suppressed_units , 1 ) = reg2_FRAs(  reg2_suppressed_units , 2 ) ;
        reg2_FRAs = reg2_FRAs( 1:size( reg2_FRAs , 1 ) , 1 ) ;
        reg2_FRAs = cellfun( @(x) x(:)' , reg2_FRAs , 'UniformOutput' , false ) ;
        reg2_FRAs = cell2mat( reg2_FRAs )' ; 
        
        
        cross_reg_sig_corr = corr( reg1_FRAs , reg2_FRAs ) ;
        h_fig_mat = figure();
        h_ax_mat = axes( 'Parent' , h_fig_mat ) ;
        imagesc( h_ax_mat , cross_reg_sig_corr ) ;
        caxis( h_ax_mat , [-1,1] )
        title( [ brain_regions{ 1 ,  combo(1) } 'x' brain_regions{ 1 ,  combo(2) } ' - Signal Correlations ' ] )
         
        cross_reg_corrs_vec = cross_reg_sig_corr( : ) ;
        across_reg_corrs{ kk , 1 } = cross_reg_corrs_vec ;
       
        shuffled_dist = [] ;
        for aa = 1 : iters_for_shuffle
        
            [~, FRA1_randperm] = sort(rand(size( reg1_FRAs ,1 ) , size( reg1_FRAs , 2 ) )  , 1 )  ;
            [~, FRA2_randperm] = sort(rand(size( reg2_FRAs ,1 ) , size( reg2_FRAs , 2 ) )  , 1 )  ;
            
            for bb = 1 : size( reg1_FRAs , 2 ) 
                
                shuf_FRA1( : , bb ) = reg1_FRAs( FRA1_randperm( : , bb ) , bb ) ; 
            
            end
            
            for cc = 1 : size( reg2_FRAs , 2 ) 
                
                shuf_FRA2( : , cc ) = reg2_FRAs( FRA2_randperm( : , cc ) , cc ) ; 
            
            end
            
            shuf_corr = corr( shuf_FRA1 , shuf_FRA2 ) ;
            shuffled_dist = [ shuffled_dist ; shuf_corr(:) ] ;
        
        end
        
        across_reg_corrs{ kk , 3 } = shuffled_dist ; 
        
        
        by_date_corrs = [] ;        
        h_fig = figure() ;
       
        for mm = 1 : size( rec_dates , 1 )

            date_table = data_table( data_table.rec_date == rec_dates( mm , : ) , : ) ;
            
            reg1_table = date_table( date_table.acronym == brain_regions{ 1 , combo(1) } , : ) ;
            reg1_FRAs = [ reg1_table.FRA , reg1_table.suppressive_FRA ] ;
            reg1_suppressed_units = cellfun( @(x) isempty(x),  reg1_FRAs(:,1) ) ;
            reg1_FRAs(  reg1_suppressed_units , 1 ) = reg1_FRAs(  reg1_suppressed_units , 2 ) ;
            reg1_FRAs = reg1_FRAs( 1:size( reg1_FRAs , 1 ) , 1 ) ;
            reg1_FRAs = cellfun( @(x) x(:)' , reg1_FRAs , 'UniformOutput' , false ) ;
            reg1_FRAs = cell2mat( reg1_FRAs )' ;  

            reg2_table = date_table( date_table.acronym == brain_regions{ 1 , combo(2) } , : ) ;
            reg2_FRAs = [ reg2_table.FRA , reg2_table.suppressive_FRA ] ;
            reg2_suppressed_units = cellfun( @(x) isempty(x),  reg2_FRAs(:,1) ) ;
            reg2_FRAs(  reg2_suppressed_units , 1 ) = reg2_FRAs(  reg2_suppressed_units , 2 ) ;
            reg2_FRAs = reg2_FRAs( 1:size( reg2_FRAs , 1 ) , 1 ) ;
            reg2_FRAs = cellfun( @(x) x(:)' , reg2_FRAs , 'UniformOutput' , false ) ;
            reg2_FRAs = cell2mat( reg2_FRAs )' ; 
            
            if ~isempty( reg1_FRAs ) && ~isempty( reg2_FRAs ) 
                
                cross_reg_date_sig_corr = corr( reg1_FRAs , reg2_FRAs ) ; 
                h_ax = subplot( 2 , ceil( size( rec_dates , 1 )/ 2 ) , mm ) ;
                imagesc( h_ax , cross_reg_date_sig_corr ) ;
                
                all_dates_corr_cell{ combo(1) , combo(2) , mm } = cross_reg_date_sig_corr ; 
                caxis( h_ax, [-1,1] )
                title( h_ax, char( rec_dates( mm , : ) ) ) ;
                cross_reg_date_sig_corr = cross_reg_date_sig_corr( : ) ;
                by_date_corrs = [ by_date_corrs ; cross_reg_date_sig_corr ] ;
              
                
                shuffled_bydate_dist = [] ;
                for aa = 1 : iters_for_shuffle

                    [~, FRA1_bydate_randperm] = sort(rand(size( reg1_FRAs ,1 ) , size( reg1_FRAs , 2 ) )  , 1 )  ;
                    [~, FRA2_bydate_randperm] = sort(rand(size( reg2_FRAs ,1 ) , size( reg2_FRAs , 2 ) )  , 1 )  ;

                    for bb = 1 : size( reg1_FRAs , 2 ) 

                        shuf_bydate_FRA1( : , bb ) = reg1_FRAs( FRA1_bydate_randperm( : , bb ) , bb ) ; 

                    end

                    for cc = 1 : size( reg2_FRAs , 2 ) 

                        shuf_bydate_FRA2( : , cc ) = reg2_FRAs( FRA2_bydate_randperm( : , cc ) , cc ) ; 

                    end

                    shuf_bydate_corr = corr( shuf_bydate_FRA1 , shuf_bydate_FRA2 ) ;
                   shuffled_bydate_dist = [ shuffled_bydate_dist ; shuf_bydate_corr(:) ] ;

                end

                across_reg_corrs{ kk , 4 } = shuffled_dist ;      

            end
            
        end   
           
        across_reg_corrs{ kk , 2 } = by_date_corrs ; 
        
    end


end

