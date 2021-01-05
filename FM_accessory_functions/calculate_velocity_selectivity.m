function [ mean_FR_lin_cell ,  Fano_FR_lin_cell , mean_FR_log_cell , Fano_FR_log_cell ,...
           mean_FR_lin_cell_no_dir , Fano_FR_lin_cell_no_dir , mean_FR_log_cell_no_dir ,...
           Fano_FR_log_cell_no_dir ] = calculate_velocity_selectivity( data_table ,...
           protocol_struct , PRE_STIM_MSEC  , WINDOW_SIZE , calculation_mode )
      
    % The function CALCULATE_VELOCITY_SELECTIVITY receives the data table,
    % the FM protocol struct, the pre stimulus baseline time
    % (PRE_STIM_MSEC) and the response window size (WINDOW_SIZE) and
    % calculates the degree of FM velocity selectivy by calculating the
    % summarized response to FMs of different velocities.
    % The function also receives the parameter "calculation_mode" which 
    % refers to the way we want to sum the spikes to calculate velocity 
    % selectivity:
    % calculation_mode = 1 - calculate over all the stimulus time and then
    % average.
    % calculation_mode = 2 - calculate over the optimal window 
    % The function calculates for every FM velocity (separately for linear
    % and logarithmic FMs) the mean response and the fano factor of the
    % response.
    % In addition, the function calculates the mean response and
    % fano-factor when taking into consideration FM speeds with no
    % direction (rising/falling).
    
    % Define parameters and constants
    MSEC_IN_SEC = 1000 ;
    Hz_IN_kHz = 1000 ;
    RED = [	174, 96, 88 ] ./ 256 ; 
    BLUE = [ 1,85,151] ./ 256 ; 
    
    n_reps = protocol_struct.n_reps ;
    n_stims_per_prot = 2 * length( protocol_struct.oct_speeds ) ;
    ind_vec = [ 1 : n_reps : n_reps * n_stims_per_prot ] ;
    stim_times = round( [ protocol_struct.time_of_stim , fliplr( protocol_struct.time_of_stim ) ] .* MSEC_IN_SEC )' ; 
    
    % Prepare empty data structures
    mean_FR_lin_cell = cell( size( data_table ,1 ) , 1 ) ;
    mean_FR_log_cell = cell( size( data_table ,1 ) , 1 ) ;
    Fano_FR_lin_cell = cell( size( data_table ,1 ) , 1 ) ;
    Fano_FR_log_cell = cell( size( data_table ,1 ) , 1 ) ;
    
    mean_FR_lin_cell_no_dir = cell( size( data_table ,1 ) , 1 ) ;
    mean_FR_log_cell_no_dir = cell( size( data_table ,1 ) , 1 ) ;   
    Fano_FR_lin_cell_no_dir = cell( size( data_table ,1 ) , 1 ) ;
    Fano_FR_log_cell_no_dir = cell( size( data_table ,1 ) , 1 ) ;   
      
  for kk = 1 : size( data_table , 1 ) 

        mean_FR_lin = zeros( 1 , n_stims_per_prot ) ;
        mean_FR_log = zeros( 1 , n_stims_per_prot ) ;
        Fano_FR_lin = zeros( 1 , n_stims_per_prot ) ;      
        Fano_FR_log = zeros( 1 , n_stims_per_prot ) ;      

        log_resp_mat = data_table.log_responses{ kk , 1 } ;        
        lin_resp_mat = data_table.lin_responses{ kk , 1 } ;

        lin_windows = data_table.best_lin_win{ kk , 1 } ; 
        log_windows = data_table.best_log_win{ kk , 1 } ; 
   
            for mm = 1 : n_stims_per_prot
            
                lin_trials = lin_resp_mat( [ ind_vec(mm) : ind_vec(mm) + n_reps - 1 ] , : ) ;
                log_trials = log_resp_mat( [ ind_vec(mm) : ind_vec(mm) + n_reps - 1 ] , : ) ;

                if calculation_mode == 1 
                
                    total_spikes_lin = sum( lin_trials( : , PRE_STIM_MSEC +1 : PRE_STIM_MSEC + stim_times( mm )  ) , 2 )./ stim_times( mm ) .* Hz_IN_kHz ;
                    total_spikes_log = sum( log_trials( : , PRE_STIM_MSEC + 1 : PRE_STIM_MSEC + stim_times( mm ) ) , 2 )./stim_times( mm ) .* Hz_IN_kHz ;

                elseif calculation_mode == 2
                    
                    total_spikes_lin = sum( lin_trials( : , PRE_STIM_MSEC + lin_windows(mm) : PRE_STIM_MSEC + lin_windows(mm) + WINDOW_SIZE - 1 ) , 2 )./ WINDOW_SIZE .* Hz_IN_kHz;
                    total_spikes_log = sum( log_trials( : , PRE_STIM_MSEC + log_windows(mm) : PRE_STIM_MSEC + log_windows(mm) + WINDOW_SIZE - 1 ) , 2 )./ WINDOW_SIZE .* Hz_IN_kHz ;
                    
                end
                 
                mean_spikes_lin = mean( total_spikes_lin , 1 ) ;
                FF_spikes_lin = var( total_spikes_lin , 1 ) ./ mean_spikes_lin ;

                mean_spikes_log = mean( total_spikes_log , 1 ) ;
                FF_spikes_log = var( total_spikes_log , 1 ) ./ mean_spikes_log ;

                mean_FR_lin( 1 , mm ) = mean_spikes_lin ;
                mean_FR_log( 1 , mm ) = mean_spikes_log ;
                Fano_FR_lin( 1 , mm ) = FF_spikes_lin ;
                Fano_FR_log( 1 , mm ) = FF_spikes_log ;
                
            end
        
            mean_FR_lin_cell{ kk , 1 } = mean_FR_lin ;
            mean_FR_log_cell{ kk , 1 } = mean_FR_log ;
            Fano_FR_lin_cell{ kk , 1 } = Fano_FR_lin ;
            Fano_FR_log_cell{ kk , 1 } = Fano_FR_log ; 
            
            mean_FR_no_dir_lin = reshape( mean_FR_lin , length( protocol_struct.oct_speeds ) , 2 ) ;
            mean_FR_no_dir_lin( : , 2 ) = flipud( mean_FR_no_dir_lin ( : , 2 ) ) ;
            mean_FR_no_dir_lin = sum( mean_FR_no_dir_lin , 2 )./2 ;
            
            mean_FR_no_dir_log = reshape( mean_FR_log , length( protocol_struct.oct_speeds ) , 2 ) ;
            mean_FR_no_dir_log( : , 2 ) = flipud( mean_FR_no_dir_log ( : , 2 ) ) ;
            mean_FR_no_dir_log = sum( mean_FR_no_dir_log , 2 )./2 ;
            
            Fano_FR_no_dir_lin = reshape( Fano_FR_lin , length( protocol_struct.oct_speeds ) , 2 ) ;
            Fano_FR_no_dir_lin( : , 2 ) = flipud( Fano_FR_no_dir_lin ( : , 2 ) ) ;
            Fano_FR_no_dir_lin = sum( Fano_FR_no_dir_lin , 2 )./2 ;
            
            Fano_FR_no_dir_log = reshape( Fano_FR_log , length( protocol_struct.oct_speeds ) , 2 ) ;
            Fano_FR_no_dir_log( : , 2 ) = flipud( Fano_FR_no_dir_log ( : , 2 ) ) ;
            Fano_FR_no_dir_log = sum( Fano_FR_no_dir_log , 2 )./2 ;
            
            mean_FR_lin_cell_no_dir{ kk ,1 } = mean_FR_no_dir_lin' ; 
            mean_FR_log_cell_no_dir{ kk ,1 } = mean_FR_no_dir_log' ; 
            Fano_FR_lin_cell_no_dir{ kk , 1 } = Fano_FR_no_dir_lin' ;
            Fano_FR_log_cell_no_dir{ kk , 1 } = Fano_FR_no_dir_log' ; 
            
    end    

end

