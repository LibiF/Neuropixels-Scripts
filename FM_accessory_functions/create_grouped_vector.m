function [ data_vec , group_vec ] = create_grouped_vector( data_cell ) 

    % This function is used to create a long data vector (data_vec) along
    % with a matching group assignment vector (group_vec).
    % The data vector and group assignment vectors can be used to run
    % statistical analysis such as KW

    n_groups = size( data_cell , 1 ) ;
    data_vec = cell2mat( cellfun( @(x) x(:) , data_cell  , 'UniformOutput' , false ) ) ;
    group_sizes = cellfun( @(x) size( x, 1 ) , data_cell ) ;
    group_vec = [] ;
    for kk = 1 : n_groups
    
        group_vec = [ group_vec ; kk.* ones( group_sizes(kk,1 ) ,1 ) ] ;
        
    end
    
end

