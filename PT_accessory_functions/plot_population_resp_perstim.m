function stimuli_PSTHs = plot_population_resp_perstim( data_table , brain_regions , PRE_STIM_MSEC, protocol_struct, SMOOTH_WIN , my_colormap , title_text )

    % The function PLOT_POPULATION_RESP_PERSTIM plots the mean population
    % response within each brain region averaged over stimuli.
    % The function plots the mean response+SEM using the function
    % "plot_areaerrorbar" by Victor Martinez-Cagigal   
    % https://www.mathworks.com/matlabcentral/fileexchange/58262-shaded-area-error-bar-plot

    MSEC_IN_SEC = 1000 ;
    
    % Preparatory step to plotting
    options.handle = figure() ;
    h_ax= axes( 'Parent' , options.handle ) ; 
    new_legend = cell( size( brain_regions , 2 )*2 , 1 ) ;

    % Extract parameters from protocol
    n_stimuli = length( protocol_struct.atten ) * length( protocol_struct.freqs ) ;
    n_reps = protocol_struct.n_reps ; 
    raster_size = n_stimuli * n_reps ; 
    
    % Define empty data structures
    stimuli_PSTHs = cell( n_stimuli , size( brain_regions , 2 ) ) ; 
    
    for kk = 1 : size( brain_regions , 2 )

        new_legend{ 2*(kk-1) + 1 , 1 } = brain_regions{ 1,  kk } ;
        new_legend{ 2*(kk-1) + 2 , 1 } = 'SEM' ;
        
        wanted_table = data_table( data_table.acronym == brain_regions{ 1 , kk } , : ) ;
        n_units = size( wanted_table , 1 ) ;
      
        rasters = cell2mat( wanted_table.PT_raster ) ;
        
        trial_vecs = repmat( [ 1 : n_reps ]' , 1 , n_units ) ;
        units_vec = 0 : 1 : n_units - 1 ; 
        all_stim_trials = trial_vecs + ( units_vec .* raster_size ) ;
        
        for tt = 1 : n_stimuli 
        
            stim_trials = all_stim_trials + (tt-1) * n_reps ;
            stim_trials = stim_trials(:) ; 
            
            raster_trials = rasters( stim_trials , : ) ; 
            mean_curves = nanmean( raster_trials ) ;
            mean_curves = mean_curves - nanmean( mean_curves( 1 : PRE_STIM_MSEC ) ) ;
            mean_curves = smooth( mean_curves , SMOOTH_WIN )' ;
            stimuli_PSTHs{ tt, kk } = mean_curves ; 
            
        end    
            
        gca() ;
        options.color_line = my_colormap(kk,:) ;
        options.error = 'sem' ;
        options.color_area = my_colormap(kk,:) ; 
        options.alpha = 0.5 ; 
        options.line_width  = 2.5 ;
        options.x_axis = 1 : length( mean_curves ) ; 
        plot_areaerrorbar( cell2mat( stimuli_PSTHs(1:end,kk) ), options ) ; 
        hold( h_ax , 'on' ) ; 
        
    end
              
     title( h_ax, title_text ) ;
end

