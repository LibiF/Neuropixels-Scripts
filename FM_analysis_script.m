%% Script for analysis of FM responses - Feigin et al. 2020

% Load data table
addpath( 'E:\Neuropixels_scripts\New_scripts\may2020' ) ;
NV_table = 'E:\Final_datasets_forthesis\Awake\Naives\NV_pop_table.mat' ;
NV_table = load( NV_table ) ;
NV_table = NV_table.sorted_population_table ;

% Load FM stimulus protocol
protocol_file_path = 'E:\Neuropixels_scripts\Awake_analysis_scipts\' ;
protocol_file_name = 'FMprotocol.mat' ;
protocol_struct = load( fullfile( protocol_file_path , protocol_file_name ) ) ; 

% Define colormap and useful figure colors
my_colormap = [ 54, 144, 192 ; 140 , 107, 177 ; 102, 194, 164 ]./256 ;
RED = [	174, 96, 88 ] ./ 256 ; 
BLUE = [ 1,85,151] ./ 256 ; 
%% Pre-processing

% Remove units with less than "SPIKE_THRESH" spikes in FM protocol/missing
% protocol
SPIKE_THRESH = 10 ;
FM_columns = [ 15 , 16 ] ;
responsive_units = remove_silent_neurons_fromFMs( NV_table , SPIKE_THRESH , FM_columns ) ; 

% Change the Acronym of TeA just for purposes of visibility and graph
% legends
responsive_units.acronym( responsive_units.acronym == 'TEa' ) = 'TeA' ;

%% Split to SUs and MUs 

responsive_SUs = responsive_units( responsive_units.clus_classification == 2 , : ) ;
responsive_MUs = responsive_units( responsive_units.clus_classification == 1 , : ) ;

%% Keep only desired regions (AUDp,AUDv,TeA)

responsive_SUs = responsive_SUs( responsive_SUs.acronym == 'TeA' | responsive_SUs.acronym == 'AUDp' | responsive_SUs.acronym == 'AUDv' , : ) ; 
responsive_MUs = responsive_MUs( responsive_MUs.acronym == 'TeA' | responsive_MUs.acronym == 'AUDp' | responsive_MUs.acronym == 'AUDv' , : ) ; 

%% Find optimal (max and min) response window for each FM

WINDOW_SIZE = 50 ;                              % Set response window size 
POST_STIM_WINDOW = 100 ;                        % Max time post stimulus end
PRE_STIM_MSEC = 100 ;                           % Time of pre-stim baseline activity
brain_regions =  {'AUDp' , 'AUDv' , 'TeA' }  ;  % Desired brain regions

% Find maximal response window to define excited units
 [ responsive_SUs.best_lin_win, responsive_SUs.best_log_win , responsive_SUs.lin_sig_mat ,...
   responsive_SUs.log_sig_mat , responsive_SUs.lin_lat_topeak,...
   responsive_SUs.log_lat_topeak] = find_best_window_in_FMresponses( responsive_SUs ,...
                protocol_struct , WINDOW_SIZE , POST_STIM_WINDOW , PRE_STIM_MSEC ) ; 
            
% Find minimal response window to define suppressed units
 [ responsive_SUs.min_lin_win, responsive_SUs.min_log_win , responsive_SUs.lin_sig_inh_mat ,...
   responsive_SUs.log_sig_inh_mat ] = find_min_window_in_FMresponses( responsive_SUs ,...
                protocol_struct , WINDOW_SIZE , POST_STIM_WINDOW , PRE_STIM_MSEC ) ; 
%% Summary statistics for excited and suppressed units

responsive_SUs.is_excited_log = sum( cell2mat( responsive_SUs.log_sig_mat ) , 2 ) > 0 ;
responsive_SUs.is_excited_lin = sum( cell2mat( responsive_SUs.lin_sig_mat ) , 2 ) > 0 ;

responsive_SUs.is_inhibited_log = sum( cell2mat( responsive_SUs.log_sig_inh_mat ) , 2 ) > 0 ;
responsive_SUs.is_inhibited_lin = sum( cell2mat( responsive_SUs.lin_sig_inh_mat ) , 2 ) > 0 ;

FM_rez.log_inhib_only = find( responsive_SUs.is_excited_log == 0 & responsive_SUs.is_inhibited_log > 0 ) ;
FM_rez.lin_inhib_only = find( responsive_SUs.is_excited_lin == 0 & responsive_SUs.is_inhibited_lin > 0 ) ;

IF_PRINT = 1 ;
title_text = 'Fraction of responsive SUs' ;
[ FM_rez.significant_response_mat_log , FM_rez.significant_response_nums_log ,...
  FM_rez.significant_response_mat_lin , FM_rez.significant_response_nums_lin ] = response_distribution_FMs( responsive_SUs ,...
  brain_regions , IF_PRINT , title_text ) ; 

%% Summary statistics over all FMs

% Summarize response types over all FMs
FM_rez.response_sum_mat = fraction_exc_inh_allFMs( responsive_SUs , brain_regions ) ;
% Define columns to run stats on:
% 1 - Excited
% 2 - Suppressed
% 3 - E & S
% 4 - Non responsive
TOTAL_CELLS_COL = 5 ;
STATS_OVER_COLS = 1:2 ;

% Summarize responses and run KW statistics
response_sum_vec = [] ; group_vec = [] ;
for kk =1 : size( brain_regions ,2 )
    
    response_sum_vec =  [ response_sum_vec ; ones( sum( FM_rez.response_sum_mat( kk , STATS_OVER_COLS ) , 2) , 1 ) ] ;
    response_sum_vec =  [ response_sum_vec ; zeros( ( FM_rez.response_sum_mat( kk ,TOTAL_CELLS_COL ) - sum( FM_rez.response_sum_mat( kk , STATS_OVER_COLS ) , 2 ) ) , 1 ) ] ;
    group_vec = [ group_vec ; kk.* ones(FM_rez.response_sum_mat( kk ,TOTAL_CELLS_COL ) , 1 ) ] ;
    
end    

[p,tbl,stats] = kruskalwallis( response_sum_vec , group_vec )
[c,~,~,gnames] = multcompare(stats)

%% Clear variables

clear response_sum_vec group_vec p tbl stats c gnames
clc

%% Quantify FM response bandwidth

% Calculate FM bandwidth and plot CDFs
RESP_MODE = 1 ;
[ FM_rez.log_signif_exc , FM_rez.lin_signif_exc ] = plot_significance_CDF_per_region( responsive_SUs ,...
  brain_regions , my_colormap , RESP_MODE ) ; 
[ FM_rez.signif_FMs_per_region ] = plot_significance_CDF_per_region_allstimuli( responsive_SUs ,...
  brain_regions , my_colormap ) ;

% KW for fraction signiifcant FMs per region
[ signif_FMs_per_reg_data ,signif_FMs_per_reg_groups ] = create_grouped_vector( FM_rez.signif_FMs_per_region ) ;
[p,tbl,stats] = kruskalwallis( signif_FMs_per_reg_data , signif_FMs_per_reg_groups )
[c,~,~,gnames] = multcompare(stats)

%% Clear variables
clear signif_FMs_per_reg_data signif_FMs_per_reg_groups p tbl stats c gnames
clc

%% Compare excitation for linear vs. logarithmic FMs

% Calculate the difference in responses to linear and logarithmic FMs
FM_rez.linlog_signifresp_diff = paired_responses_linlog( responsive_SUs ,...
 brain_regions , my_colormap ) ;
title_text = 'lin-log FMs number of significant stimuli' ;
create_bp_from_cells( FM_rez.linlog_signifresp_diff( 1:length(brain_regions) , 1 ) ,...
 brain_regions , my_colormap , title_text )

% Validate that distributions have 0 medians
for kk = 1 : size( brain_regions , 2 )
   
    disp ( brain_regions{ 1 , kk } ) ;
    [ p , h ] = signrank(FM_rez.linlog_signifresp_diff{kk,1})

end

% KW for fraction signifcant FMs per region
[ linlog_signifresp_diff_data ,linlog_signifresp_diff_groups ] = create_grouped_vector( FM_rez.linlog_signifresp_diff ) ;
[p,tbl,stats] = kruskalwallis( linlog_signifresp_diff_data , linlog_signifresp_diff_groups )
[c,~,~,gnames] = multcompare(stats)

%% Clear variables 

clear linlog_signifresp_diff_data linlog_signifresp_diff_groups p tbl stats c gnames
clc

%% Calculate lifetime sparseness for FM stimuli

% Calculate lifetime sparseness for FMs
[ responsive_SUs.S_stat_lin , responsive_SUs.S_stat_log , responsive_SUs.S_stat_tot ] = calculate_FM_lifetimesparsness( responsive_SUs ,...
  protocol_struct , PRE_STIM_MSEC , WINDOW_SIZE ) ;

% Prepare empty data structures to summarize results across regions
FM_rez.S_stat_lin = cell( size( brain_regions , 2 ) ,1 ) ;
FM_rez.S_stat_log = cell( size( brain_regions , 2 ) ,1 ) ;
FM_rez.S_stat_tot = cell( size( brain_regions , 2 ) ,1 ) ;

for mm = 1 : size( brain_regions , 2 )
    
    reg_units = responsive_SUs( responsive_SUs.acronym == brain_regions{ 1 , mm } , : ) ;
    reg_units_all = reg_units( reg_units.is_excited_lin == 1 | reg_units.is_excited_log == 1 , : ) ;
    reg_units_lin = reg_units( reg_units.is_excited_lin == 1 , : ) ;
    reg_units_log = reg_units( reg_units.is_excited_log == 1 , : ) ;

    FM_rez.S_stat_lin{ mm , 1 } = reg_units_lin.S_stat_lin ;
    FM_rez.S_stat_log{ mm , 1 } = reg_units_log.S_stat_lin ;
    FM_rez.S_stat_tot{ mm , 1 } = reg_units_all.S_stat_tot ; 
    
end

% Plot boxplot
title_text = 'Lifetime Sparseness - all stimuli' ;
DATA_COL = 1 ;
create_bp_from_cells( FM_rez.S_stat_tot , brain_regions , my_colormap , title_text , DATA_COL  ) ;

% KW for total FMs lifetime sparseness
[ tot_lifetimesparse_data , tot_lifetimesparse_groups ] = create_grouped_vector( FM_rez.S_stat_tot ) ;
[p,tbl,stats] = kruskalwallis( tot_lifetimesparse_data , tot_lifetimesparse_groups )
[c,~,~,gnames] = multcompare(stats)

%% Clear variables

clear tot_lifetimesparse_data  tot_lifetimesparse_groups p tbl stats c gnames
clc

%% Calculate Population sparseness

% Preparatory phase - extract parameters from the FM protocol
MSEC_IN_SEC = 1000 ;
n_stims_per_dir = length( protocol_struct.oct_speeds) ;
stim_times = [ protocol_struct.time_of_stim , fliplr( protocol_struct.time_of_stim ) ] .* MSEC_IN_SEC' ; 
lin_slopes = ( protocol_struct.upper_freq - protocol_struct.bottom_freq ) ./ stim_times ;
lin_slopes( 1 : n_stims_per_dir ) = - lin_slopes( 1: n_stims_per_dir ) ;
lin_slopes = round( lin_slopes .* 100 ) ./ 100 ; 
log_slopes =[ - protocol_struct.oct_speeds , fliplr( protocol_struct.oct_speeds ) ] ; 

% Create empty data structures
FM_rez.lin_pop_spaseness = cell( size( brain_regions , 2 ) , 1 ) ;
FM_rez.log_pop_spaseness = cell( size( brain_regions , 2 ) , 1 ) ;

% Calculate the population sparseness across regions
for mm = 1 : size( brain_regions , 2 ) 

    reg_units = responsive_SUs( responsive_SUs.acronym == brain_regions{ 1 , mm } , : ) ;
    lin_sig_mat = cell2mat( reg_units.lin_sig_mat ) ;
    log_sig_mat = cell2mat( reg_units.log_sig_mat ) ;
    
    lin_sig_fract = sum( lin_sig_mat , 1 ) ./ size ( lin_sig_mat , 1 ) ;
    log_sig_fract = sum( log_sig_mat , 1 ) ./ size ( log_sig_mat , 1 ) ;
    
    FM_rez.lin_pop_spaseness{ mm ,1 } = lin_sig_fract ;
    FM_rez.log_pop_spaseness{ mm ,1 } = log_sig_fract ;
    
end

% Plot linear FM population sparseness results
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

% Plot logarithmic FM population sparseness results
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

%% Calculate direction selectivity index

% Calculate the direction selectivity index for all units
PRE_STIM_MSEC = 100 ;
[ responsive_SUs.total_selectivity_index_log ,  responsive_SUs.total_selectivity_index_lin , ...
  responsive_SUs.total_DSI ] = calculate_direction_selectivity( responsive_SUs ,...
  protocol_struct , PRE_STIM_MSEC ,WINDOW_SIZE) ; 

% Plot DSI distributions
 [FM_rez.total_selec_ind_log , FM_rez.total_selec_ind_lin ,...
  FM_rez.total_DSI_allstims ] = plot_selectivity_ind_distribution( responsive_SUs ,...
  brain_regions , my_colormap );

% Plot boxplots for each DSI parameter calculated
title_text = 'Total lin S.I'
create_bp_from_cells( FM_rez.total_selec_ind_lin( 1:size( brain_regions , 2 ) , 1 ) ,...
                      brain_regions , my_colormap , title_text )
title_text = 'Total log S.I'
create_bp_from_cells( FM_rez.total_selec_ind_log( 1:size( brain_regions , 2 ) , 1 ) ,...
                      brain_regions , my_colormap , title_text )
title_text = 'Total DSI - all stimuli'
create_bp_from_cells( FM_rez.total_DSI_allstims( 1:size( brain_regions , 2 ) , 1 ) ,...
                      brain_regions , my_colormap , title_text )
create_scatterbp_forcells( FM_rez.total_DSI_allstims( 1:size( brain_regions , 2 ) , 1 ) ,...
                           brain_regions , my_colormap , title_text )

% Calculate signrank stat to see whether distributions are centered around 
for mm = 1 : size( brain_regions , 2 )
    
    disp( brain_regions{ 1 , mm } ) ;
    [ p , h , stats ] = signrank( FM_rez.total_DSI_allstims{ mm , 1 } )
    
end

% Calculate distribution of absolute DSI values
FM_rez.total_DSI_absval = cellfun( @(x) abs(x) , FM_rez.total_DSI_allstims ,'UniformOutput' , false) ;
title_text = 'DSI absolute value' ;
create_bp_from_cells( FM_rez.total_DSI_absval( 1:3 , 1 ) , brain_regions , my_colormap , title_text )

% KW for direction selectivity all stims
[ DSI_allstims_data , DSI_allstims_groups ] = create_grouped_vector(  FM_rez.total_DSI_allstims ) ;
[p,tbl,stats] = kruskalwallis( DSI_allstims_data , DSI_allstims_groups )
[c,~,~,gnames] = multcompare(stats)

%% KW for absolute DSI all stims
clear DSI_allstims_data  DSI_allstims_groups p tbl stats c gnames
clc

[ DSI_abs_allstims_data , DSI_abs_allstims_groups ] = create_grouped_vector(  FM_rez.total_DSI_absval ) ;
[p,tbl,stats] = kruskalwallis( DSI_abs_allstims_data , DSI_abs_allstims_groups )
[c,~,~,gnames] = multcompare(stats)

%% Clean KW
clear DSI_abs_allstims_data DSI_abs_allstims_groups p tbl stats c gnames
clc

%% Direction selectivity vs. sparseness

% Screen for FM selective units
FM_selec_SUs = responsive_SUs( sum( [responsive_SUs.is_excited_lin ,...
    responsive_SUs.is_excited_log , responsive_SUs.is_inhibited_lin , responsive_SUs.is_inhibited_log ] , 2 ) > 0 , : ) ;

% Plot absolute value of DSI vs. the FM lifetime sparseness
h_fig = figure() ;
h_ax = axes('Parent' ,h_fig ) ;
for kk = 1 : size( brain_regions , 2 ) 
    
    reg_units = FM_selec_SUs( FM_selec_SUs.acronym == brain_regions{ 1 , kk } , : ) ;
    scatter( h_ax , reg_units.S_stat_tot , abs( reg_units.total_DSI ) , 20 , my_colormap(kk,:) )  ;
    hold( h_ax ,'on' ) ;
    
end

% Fit a linear model for the data and plot it over the scatter plot
mdl = fitlm(FM_selec_SUs.S_stat_tot,abs( FM_selec_SUs.total_DSI ))
plot( h_ax , [0:0.1:1] ,[0:0.1:1].*mdl.Coefficients.Estimate(2) + mdl.Coefficients.Estimate(1), '--k' ) ;
% plot( h_ax , [0:0.1:1] ,[0:0.1:1].*mdl.Coefficients.Estimate(1) , '--k' ) ;
xlabel( h_ax , 'Lifetime Sparseness' ) ;
ylabel( h_ax , '|DSI|' ) ; 
[R,p] = corrcoef(FM_selec_SUs.S_stat_tot,abs( FM_selec_SUs.total_DSI ))

%% Calculate velocity selectivity

% Define calculation mode for "calculate_velocity_selectivity"
CALC_MODE = 2 ;
% Calculate velocity and speed selectivity
[ responsive_SUs.mean_FR_lin , responsive_SUs.Fano_FR_lin ,  responsive_SUs.mean_FR_log ,...
  responsive_SUs.Fano_FR_log , responsive_SUs.FR_lin_nodir , responsive_SUs.FF_lin_nodir ,...
  responsive_SUs.FR_log_nodir , responsive_SUs.FF_log_nodir ] = calculate_velocity_selectivity( responsive_SUs ,...
  protocol_struct , PRE_STIM_MSEC  , WINDOW_SIZE , CALC_MODE ) ; 

% Reorganize data vectors according to "reordered_FMs"
reordered_FMs = [ n_stims_per_dir : -1 : 1 , 2 * n_stims_per_dir : -1 : n_stims_per_dir + 1 ] ; 
mean_FR_log = cell2mat( responsive_SUs.mean_FR_log ) ;           
responsive_SUs.reordered_log_FR = num2cell( mean_FR_log(: , reordered_FMs ) , 2 ) ; 
mean_FR_lin = cell2mat( responsive_SUs.mean_FR_lin ) ;           
responsive_SUs.reordered_lin_FR = num2cell( mean_FR_lin(: , reordered_FMs ) , 2 ) ; 

%% Plot the speed selectivity results and statistical comparisons

% For logarithmic FM sweeps
DATA_COL = 45 ;         
DATA_SUB_COL = 1 ; 
SMOOTH_WIN = 1 ;
IS_NORMALIZED = 1 ;
title_text = 'Mean population FR per sweep speed - log sweeps ' ; 
FM_rez.log_FRperspeed_cell = plot_population_resp_FMspeed( responsive_SUs( responsive_SUs.is_excited_log == 1, : ) ,...
  brain_regions , IS_NORMALIZED , DATA_COL , DATA_SUB_COL , my_colormap , title_text ) ;
title_text = 'Log sweeps signif mat' ; 
[ FM_rez.logFMs_speed_ttest_scores , FM_rez.logFM_speed_ttest_perregion ]= compare_speedselectivity_via_ttest( FM_rez.log_FRperspeed_cell(1:size( brain_regions , 2 ),1) ,...
  brain_regions , title_text ) ;

% For linear FM sweeps
DATA_COL = 43 ;
SMOOTH_WIN = 1 ;
title_text = 'Mean population FR per sweep speed - lin sweeps ' ; 
FM_rez.lin_FRperspeed_cell = plot_population_resp_FMspeed( responsive_SUs(  responsive_SUs.is_excited_lin == 1 , : ) ,...
  brain_regions , IS_NORMALIZED, DATA_COL , DATA_SUB_COL , my_colormap , title_text ) ;
title_text = 'Lin sweeps signif mat' ; 
[ FM_rez.linFMs_speed_ttest_score , FM_rez.linFM_speed_ttest_perregion ] = compare_speedselectivity_via_ttest( FM_rez.lin_FRperspeed_cell(1:size( brain_regions , 2 ),1) ,...
  brain_regions , title_text ) ;

%% Calculate FM d prime

% Create suitable response matrices
[ logFM_region_matrices , linFM_region_matrices ] = create_FM_responsemat_fordprime( responsive_SUs ,...
        protocol_struct , brain_regions , WINDOW_SIZE , PRE_STIM_MSEC ) ;

% Calculate d primes and plot corresponding matrices
for mm = 1 : size( brain_regions , 2 )
    
    logFM_d_prime_matrices{mm,1} = findDprime( logFM_region_matrices{ mm,1 } ) ; 
    h_fig = figure() ;
    h_ax = axes('Parent' , h_fig ) ;
    imagesc( h_ax , logFM_d_prime_matrices{mm,1} ) ;
    title( h_ax , [  brain_regions{1,mm}  ' - log FM d prime' ] ) ;
    caxis( h_ax , [0 , 1.7 ] ) ;
    
    linFM_d_prime_matrices{mm,1} = findDprime( linFM_region_matrices{ mm,1 } ) ; 
    h_fig2 = figure() ;
    h_ax2 = axes('Parent' , h_fig2 ) ;
    imagesc( h_ax2 , linFM_d_prime_matrices{mm,1} ) ;
    title( h_ax2 , [  brain_regions{1,mm}  ' - lin FM d prime' ] ) ; 
    caxis( h_ax2 , [0 , 1.7 ] ) ;
    
end    

% Extract all unique d' values from symmetric matrix and plot box plot
% of d' distributions - log and linear FMs
rez_struct.logFMd_prime_vecs = cellfun( @(x) x(triu(boolean( ones( length(protocol_struct.oct_speeds )*2 ) ), 1 ) ) , logFM_d_prime_matrices , 'UniformOutput' , false ) ;
title_text = 'log FM d prime' ; 
DATA_COL = 1 ;
create_bp_from_cells( rez_struct.logFMd_prime_vecs , brain_regions , my_colormap , title_text , DATA_COL ) ;    

rez_struct.linFMd_prime_vecs = cellfun( @(x) x(triu(boolean( ones( length(protocol_struct.oct_speeds )*2 ) ), 1 ) ) , linFM_d_prime_matrices , 'UniformOutput' , false ) ;
title_text = 'lin FM d prime' ; 
DATA_COL = 1 ;
create_bp_from_cells( rez_struct.linFMd_prime_vecs , brain_regions , my_colormap , title_text , DATA_COL ) ;    

% Statistically compare pairwise d primes between AUDp and TeA
[h,p,CIs,tstat] = ttest( rez_struct.logFMd_prime_vecs{3,1} , rez_struct.logFMd_prime_vecs{1,1} ) 

%% Scatter plot for FM d prime

h_fig = figure() ;
h_ax = axes( 'Parent' , h_fig ) ;
scatter( h_ax, rez_struct.logFMd_prime_vecs{1,1} , rez_struct.logFMd_prime_vecs{3,1} , 20) ;
xlabel( h_ax , 'AUDp d prime' ) ;
ylabel( h_ax , 'TeA d prime' ) ;
title( h_ax , 'Logarithmic FMs d prime' ) ;
hold( h_ax , 'on' ) ;
plot( h_ax , 0:0.1:1.7  , 0:0.1:1.7 , 'k' ) ;
xlim( h_ax, [0,1.7] );
ylim( h_ax, [0,1.7] ) ; 

%% Plot FM figures - This is an accessory section to help generate FM rasters

oct_speeds = protocol_struct.oct_speeds ; 
n_reps = protocol_struct.n_reps ; 
MSEC_IN_SEC = 1000 ;                                                                                                          % Constant
POST_STIM_TIME = 300 ;                                                                                                       % Time post stimulus that we want to show
PSTH_BIN_SIZE = 10 ; 
lin_speeds = ( protocol_struct.upper_freq - protocol_struct.bottom_freq )./ protocol_struct.time_of_stim ; 
lin_speeds = lin_speeds./MSEC_IN_SEC ; 

% index_list = [ 79 , 68, 71 , 30 , 194, 179, 268, 206 , 272 , 275, 421 , 393, 413, 416, 395 ] ;
index_list = [ 30, 342, 378,329,380 ] ;


for mm = 1 : length( index_list )     
    kk = index_list( mm ) ;
    
    log_FM = responsive_SUs.log_responses{ kk , 1 } ;
    log_PSTH_unit = responsive_SUs.log_PSTH{ kk , 1 } ;
    lin_FM = responsive_SUs.lin_responses{ kk , 1 } ;
    lin_PSTH_unit = responsive_SUs.lin_PSTH{ kk , 1 } ;
    cluster_num = responsive_SUs.clus_num( kk , 1 ) ;
    brain_region = responsive_SUs.acronym( kk , : ) ;
    color_ind = find( contains(brain_regions,  char(brain_region ) ) ) ;
    
    
    h_fig = create_FM_fig2( log_FM, lin_FM, oct_speeds, lin_speeds, n_reps , PRE_STIM_MSEC , cluster_num , log_PSTH_unit , lin_PSTH_unit , protocol_struct , my_colormap , color_ind ) ;

end