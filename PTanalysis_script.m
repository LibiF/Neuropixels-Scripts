%% Script for Pure tone analysis - Feigin et. al. 2020

% Load data table
NV_table = 'E:\Final_datasets_forthesis\Awake\Naives\NV_pop_table.mat' ;
NV_table = load( NV_table ) ;
NV_table= NV_table.sorted_population_table ;

% Load protocol struct
protocol_file_path = 'E:\Neuropixels_scripts\Awake_analysis_scipts\' ;
protocol_file_name = 'Puretone_protocol.mat' ;
protocol_struct = load( fullfile( protocol_file_path , protocol_file_name ) ) ; 

% Defining a colormap for figures
my_colormap = [ 54, 144, 192 ; 140 , 107, 177 ; 102, 194, 164 ]./256 ;

%% Remove cells with less than some spike threshold

% Remove units with <SPIKE_THRESH spikes during entire PT protocol
SPIKE_THRESH = 10 ;
responsive_units = remove_silent_neurons( NV_table , SPIKE_THRESH ) ;

% Change the Acronym of TeA just for purposes of visibility and graph
% legends
responsive_units.acronym( responsive_units.acronym == 'TEa' ) = 'TeA' ;

%% Split to SUs and MUs & Isolate units from AUDp, AUDv and TeA only

responsive_SUs = responsive_units( responsive_units.clus_classification == 2 , : ) ;
responsive_MUs = responsive_units( responsive_units.clus_classification == 1 , : ) ;

responsive_SUs = responsive_SUs( responsive_SUs.acronym == 'TeA' | responsive_SUs.acronym == 'AUDp' | responsive_SUs.acronym == 'AUDv' , : ) ; 
responsive_MUs = responsive_MUs( responsive_MUs.acronym == 'TeA' | responsive_MUs.acronym == 'AUDp' | responsive_MUs.acronym == 'AUDv' , : ) ; 

%% Spontaneous firing rate distribution

PRE_STIM_MSEC = 100 ;                                                       % Time of baseline activity pre-trial
brain_regions =  {'AUDp' , 'AUDv' , 'TeA' }  ;                              % Define brain regions for analysis

% Calculate spontaneous firing rates and plot distributions and summary
% boxplots
IF_PRINT = 1 ;
[  responsive_SUs.spont_FR ,  responsive_SUs.spontFR_std ] = calc_spont_FR( responsive_SUs , PRE_STIM_MSEC  ) ; 
rez_struct.spont_FR_cell = plot_spont_FR(  responsive_SUs , brain_regions , IF_PRINT , my_colormap ) ;

% Create summary box plot
title_text =  'Spontaneous FR in auditory cortices' ;
create_bp_from_cells_logscale( rez_struct.spont_FR_cell , brain_regions , my_colormap , title_text ) ;

% Calculate summary KW statistics
[ spont_FR_data , spont_FR_groups ] = create_grouped_vector(rez_struct.spont_FR_cell) ;
[p,tbl,stats] = kruskalwallis(spont_FR_data, spont_FR_groups ) ;
[c,~,~,gnames] = multcompare(stats) ;

%% Clear variables
clear spont_FR_data spont_FR_groups p tbl stats c gnames
clc

%% Test fit of spontaneous FR to log-normal distributions

for kk = 1 : size( brain_regions , 2 )
    
    disp( brain_regions{1,kk} ) ;
    log_dist = log(rez_struct.spont_FR_cell{kk, 1}) ;
    log_dist( isnan( log_dist ) | isinf( log_dist ) ) = [] ;
    [h,p,k,crit]=lillietest(log_dist, 'Distribution' , 'normal')
    
end

clear h p k crit log_dist

%% Find maximal and minimal response window

WINDOW_SIZE = 50 ;                                                         % Set size of response window
RESPONSE_WINDOW = 200 ;                                                    % Max time post stim of response window
% Find excitatory and suppressive response windows
[ responsive_SUs.max_resp_wind , responsive_SUs.is_excited ] = find_max_resp_window( responsive_SUs ,...
  PRE_STIM_MSEC , WINDOW_SIZE , RESPONSE_WINDOW ) ;
[ responsive_SUs.min_resp_wind , responsive_SUs.is_inhibited ] = find_min_resp_window( responsive_SUs ,...
  PRE_STIM_MSEC , WINDOW_SIZE , RESPONSE_WINDOW ) ;

% Create a copy containing all SUs before eliminating non-responsive ones
% (will be used later in code)
all_region_SUs = responsive_SUs ;

% Leave only units showing a significt response to pure tones (either
% excitatory or suppressive
responsive_SUs = responsive_SUs( responsive_SUs.is_excited == 1 | responsive_SUs.is_inhibited == 1 , : ) ; 

% Plot bar graph summarizing response statistics in all brain regions
IF_PRINT = 1 ;
title_text = 'Fraction of responsive SUs' ;
[rez_struct.significant_response_mat_SU , rez_struct.significant_response_nums_SU ] = response_distribution( responsive_SUs , brain_regions , IF_PRINT , title_text ) ;

%% Fraction E and I unitstotal

% Calculate stats for Excited (E), Suppressed (S) or both (E&S) populations 
% by setting the variable COL_NUM to:
% 1 - Excitated
% 2 - Suppressed
% 3 - Both
COL_NUM = 3 ;
TOTAL_COL = 4 ;     % The column representing the total amount of units in each region
response_sum_vec = [] ; group_vec = [] ;
for kk =1 : size( brain_regions ,2 )
    
    response_sum_vec =  [ response_sum_vec ; ones( rez_struct.significant_response_nums_SU(kk,COL_NUM) , 1 ) ] ;
    response_sum_vec =  [ response_sum_vec ; zeros( rez_struct.significant_response_nums_SU(kk,TOTAL_COL)-rez_struct.significant_response_nums_SU(kk,COL_NUM) , 1 ) ] ;
    group_vec = [ group_vec ; kk.* ones(rez_struct.significant_response_nums_SU( kk ,TOTAL_COL ) , 1 ) ] ;
    
end    

[p,tbl,stats] = kruskalwallis( response_sum_vec , group_vec )
[c,~,~,gnames] = multcompare(stats)

%% Clear variables

clc
clear p tbl stats c gnames response_sum_vec group_vec COL_NUM TOTAL_COL 

%% Find BF for excited units

CHECK_EXCITATION = 1 ;                          
% Setting this parameter to 1 means excitatory FRAs are calculated only for excited units

% Calculate excitatory and suppressive FRAs
[ responsive_SUs.FRA , responsive_SUs.FRA_signif ,responsive_SUs.BF ,~ ] = find_BF_in_best_resp_wind_3( responsive_SUs , PRE_STIM_MSEC , protocol_struct , WINDOW_SIZE , CHECK_EXCITATION ) ;
[ responsive_SUs.suppressive_FRA , responsive_SUs.suppressive_FRA_signif ] = calculate_suppressive_FRAs(  responsive_SUs , PRE_STIM_MSEC , protocol_struct , WINDOW_SIZE ) ;

% Summarize BF data
for kk = 1 : size( brain_regions , 2 ) 

    reg_units = responsive_SUs( responsive_SUs.acronym == brain_regions{ 1 , kk } , : ) ;
    reg_units = reg_units( reg_units.is_excited ==1  , : ) ;
    rez_struct.BF{ kk , 1 } = reg_units.BF ;    
    
end

% Create Box plots summarizing BF distributions
title_text = 'BF' ;
create_BFbp_from_cells_log_scale( rez_struct.BF , brain_regions , my_colormap , title_text , protocol_struct )

% Calculate KW stats for BF
[ BF_data , BF_groups ] = create_grouped_vector(rez_struct.BF) ;
[p,tbl,stats] = kruskalwallis( BF_data, BF_groups )
[c,~,~,gnames] = multcompare(stats)


%% Clear variables

clear BF_data BF_groups p tbl stats c gnames 
clc

%% Plot units' BF by individual mice

excited_units = responsive_SUs( responsive_SUs.is_excited == 1 , : ) ; 
excited_units.rec_date = categorical( excited_units.rec_date ) ;
rec_dates = unique( excited_units.rec_date ) ;

% Note - colormap created using the brewermap tool 
by_date_colormap = brewermap(  size( rec_dates , 1 ) , 'Spectral' ) ;
rez_struct.BF_by_date = plot_BF_bymouse( excited_units , rec_dates , brain_regions , by_date_colormap , protocol_struct ) ;

clear excited_units rec_dates

%% Quantify differenced in BF distributions by date - 2Way anova

BF = [] ;
rec_ind = [] ;
brain_reg_group = [] ;
n_recs = length( unique( responsive_SUs.rec_date ) ) ;
brain_reg_short = {'p','v','T'} ;
for kk = 1 : size(brain_regions , 2 )
    
    for mm = 1: n_recs 

        BF = [ BF ; rez_struct.BF_by_date{kk,1}{1,mm} ] ;
        brain_reg_group = [ brain_reg_group ; repmat( char(brain_reg_short{kk}) , size( rez_struct.BF_by_date{kk,1}{1,mm} , 1 ) , 1 ) ] ;
        rec_ind = [ rec_ind ; mm* ones( size( rez_struct.BF_by_date{kk,1}{1,mm} , 1 ) , 1 ) ] ;
       
    end
    
end

% [p,tbl,stats] = anovan(lat_to_peak,{brain_reg_group,FM_type , FM_dir},'model','interaction','varnames',{'br-reg','lin-log','FM-dir'})
[p,tbl,stats] = anovan(BF,{brain_reg_group,rec_ind})
results = multcompare(stats,'Dimension',[1 2])

%% Compare BFs on a per recording basis (no meaning for region)

response_sum_vec = [] ; group_vec = [] ;
for kk =1 : size( brain_regions ,2 )

    reg_BFs = rez_struct.BF_by_date{kk,1} ;
    for mm = 1:n_recs
    
        response_sum_vec =  [ response_sum_vec ; reg_BFs{1,mm} ] ;
        group_vec = [ group_vec ; mm.* ones( size( reg_BFs{1,mm} , 1 ) , 1 ) ] ;
    
    end
    
end

[p,tbl,stats] = kruskalwallis( response_sum_vec , group_vec )
[c,~,~,gnames] = multcompare(stats)

%% Compare BFs on a per region basis

response_sum_vec = [] ; group_vec = [] ;
% The variable kk willl correspond to the brain region index in the cell
% "brain_regions"
for kk = 2:2

    reg_BFs = rez_struct.BF_by_date{kk,1} ;
    for mm = 1:n_recs
    
        response_sum_vec =  [ response_sum_vec ; reg_BFs{1,mm} ] ;
        group_vec = [ group_vec ; mm.* ones( size( reg_BFs{1,mm} , 1 ) , 1 ) ] ;
    
    end
    
end

[p,tbl,stats] = kruskalwallis( response_sum_vec , group_vec )
[c,~,~,gnames] = multcompare(stats)

%% Clear variables

clear p tbl stats c gnames response_sum_vec group_vec rec_dates reg_BFs n_recs brain_reg_short
clc

%% Difference in BF of adjacent units

RED = [	174, 96, 88 ] ./ 256 ; 
BLUE = [ 1,85,151] ./ 256 ; 
    
CHANNEL_RADIUS = 8 ;
rez_struct.BF_diffs_nearby_channels = calculate_octaveBFdiff_nearunits(  responsive_SUs , CHANNEL_RADIUS , [ RED ; BLUE ]   ) ;
[p,h,U_stat] = ranksum( rez_struct.BF_diffs_nearby_channels{1,1} , rez_struct.BF_diffs_nearby_channels{ 2,1} ) 

%% Clear variables

clear p h U_stat RED BLUE CHANNEL_RADIUS
clc

%% Evoked FRs for BF

% Extract evoked FRs 
responsive_SUs.evoked_BF_FR = calc_BFevoked_FRs( responsive_SUs ) ;

% Summarize evoked FRs by brain region
rez_struct.evoked_BF_FR_cell = cell( size(brain_regions , 2 ) , 1 ) ;
for kk = 1 : size( brain_regions , 2 )
    
    reg_units = responsive_SUs( responsive_SUs.acronym == brain_regions{ 1, kk } , : ) ;
    reg_units = reg_units( reg_units.is_excited == 1 , : ) ;
    rez_struct.evoked_BF_FR_cell{kk,1} = reg_units.evoked_BF_FR ;
    rez_struct.evokedvsspont_FR_cell{kk,1} = reg_units.evoked_BF_FR./ reg_units.spont_FR ;
    rez_struct.spont_FR_ofexcited_units{kk,1} = reg_units.spont_FR ;
    
end    

% Plot evoked_FR box plot and distributions
title_text = [ 'Evoked BF FR in auditory cortices ' ] ;
create_bp_from_cells_logscale(  rez_struct.evoked_BF_FR_cell , brain_regions , my_colormap , title_text ) ;
plot_evoked_FRs_distribution(  rez_struct.evoked_BF_FR_cell , brain_regions , title_text  , my_colormap ) ;

% Kruskal-Wallis for evoked FRs
[ evoked_FR_data , evoked_FR_groups ] = create_grouped_vector(rez_struct.evoked_BF_FR_cell) ;
[p,tbl,stats] = kruskalwallis( evoked_FR_data,evoked_FR_groups )
[c,~,~,gnames] = multcompare(stats)

%% Clean KW

clear evoked_FR_data evoked_FR_groups p tbl stats c gnames
clc

%% Test fit of evoked FR to log-normal distributions


for kk = 1 : size( brain_regions , 2 )
    
    disp( brain_regions{1,kk} ) ;
    log_dist = log(rez_struct.evoked_BF_FR_cell{kk, 1}) ;
    log_dist( isnan( log_dist ) | isinf( log_dist ) ) = [] ;
    [h,p,k,crit]=lillietest(log_dist, 'Distribution' , 'normal')
    
end

%% Clear variable

clear h p k crit 
clc

%% Calculate population discriminability (d') between pairs of stimuli

% Prepare data and empty data structures
PT_trial_matrices = create_PT_responsemat_fordprime( all_region_SUs , protocol_struct , brain_regions , WINDOW_SIZE , PRE_STIM_MSEC ) ;
d_prime_matrices = cell( size( brain_regions , 2 ) , 1 ) ;

% Calculate d primes and plot matrices
for mm = 1 : size( brain_regions , 2 )
    
    d_prime_matrices{mm,1} = findDprime( PT_trial_matrices{ mm,1 } ) ; 
    h_fig = figure() ;
    h_ax = axes('Parent' , h_fig ) ;
    imagesc( h_ax , d_prime_matrices{mm,1} ) ;
    title( h_ax , [  brain_regions{1,mm}  ' - PT d prime' ] ) ;
    caxis( h_ax , [0 , 1.7 ] ) ;
    
end    

% Calculate stats for d'
rez_struct.d_prime_vecs = cellfun( @(x) x(triu(boolean( ones( length(protocol_struct.freqs ) ) ), 1 ) ) , d_prime_matrices , 'UniformOutput' , false ) ;
[ h , p  , CI, stat] = ttest( rez_struct.d_prime_vecs{1,1} , rez_struct.d_prime_vecs{3,1} )

%% Clear variables
clear h p CI stat
clc
%% Latency to 1st spike of early responsive units (find significant responses during first 110 msec)

% Since maximal response window can be very late, restrict the excitation
% to the early response window (by controling the variable RESPONSE_WINDOW)
PRE_STIM_MSEC = 100 ;
WINDOW_SIZE = 50 ;
RESPONSE_WINDOW = 110 ; 
CHECK_EXCITATION = 1 ;
[ responsive_SUs.max_early_resp_wind , responsive_SUs.is_early_excited ] = find_max_resp_window( responsive_SUs , PRE_STIM_MSEC , WINDOW_SIZE , RESPONSE_WINDOW ) ;
[ responsive_SUs.early_FRA_cell , responsive_SUs.early_FRA_signif_cell ,...
  responsive_SUs.early_BF_mat, responsive_SUs.early_lat_to_BFfirstspike ] = find_BF_in_early_resp_wind( responsive_SUs ,...
  PRE_STIM_MSEC , protocol_struct , WINDOW_SIZE , CHECK_EXCITATION ) ;

% Summarize and plot boxplot
for kk = 1 : size( brain_regions , 2 ) 

    reg_units = responsive_SUs( responsive_SUs.acronym == brain_regions{ 1 , kk } , : ) ;
    rez_struct.num_early_excited(kk , 2 ) = sum( reg_units.is_excited == 1 ) ;
    reg_units = reg_units( reg_units.is_early_excited  ==1  , : ) ;
    rez_struct.num_early_excited(kk , 1 ) = size( reg_units , 1 ) ;
    rez_struct.early_BF{ kk , 1 } = reg_units.early_BF_mat ;    
    reg_units = reg_units( reg_units.early_lat_to_BFfirstspike(: , 4 ) > 0 , : ) ;
    rez_struct.lat_to_earlyBF1st_spike{ kk , 1 } = reg_units.early_lat_to_BFfirstspike ;

    
end

title_text = 'Early window latency to BF 1st spike' ;
DATA_COL = 1 ;
create_bp_from_cells( rez_struct.lat_to_earlyBF1st_spike , brain_regions , my_colormap , title_text , DATA_COL )

% KW for latency to 1st spike
latency_to_1st_spike = cellfun( @(x) x( : , 1) , rez_struct.lat_to_earlyBF1st_spike , 'UniformOutput' , false ) ;
[ latency_to_1st_spike_data ,latency_to_1st_spike_groups ] = create_grouped_vector( latency_to_1st_spike ) ;

[p,tbl,stats] = kruskalwallis( latency_to_1st_spike_data, latency_to_1st_spike_groups )
[c,~,~,gnames] = multcompare(stats)

%% Clear variables

clear latency_to_1st_spike_data latency_to_1st_spike_groups p tbl stats c gnames
clc

%% Tuning bandwidth

signif_level = 0.05 ;
% Calculate num of expected false-positives
FP_of_responses = signif_level * length( protocol_struct.atten ) * length( protocol_struct.freqs ) ;

% Calculate bandwidth
DISPLAY_MODE = 0 ;
rez_struct.bandwidth_cell = calculate_plot_bandwidth( responsive_SUs , brain_regions , FP_of_responses , protocol_struct , my_colormap , IF_PRINT , DISPLAY_MODE ) ;

% Plot summary box plot
title_text = 'Number of significant ( freq x atten ) combinations' ;
DATA_COL = 1 ;
create_bp_from_cells( rez_struct.bandwidth_cell , brain_regions , my_colormap , title_text , DATA_COL  ) ;

% KW for tuning width

[ bandwidth_data ,bandwidth_groups ] = create_grouped_vector( rez_struct.bandwidth_cell ) ;
[p,tbl,stats] = kruskalwallis( bandwidth_data , bandwidth_groups )
[c,~,~,gnames] = multcompare(stats)

%% Clear variables
clear bandwidth_data bandwidth_groups p tbl stats c gnames
clc

%% Check fraction of untuned units

% Untuned units defined as units which have an auditory response but no
% pair of frequency x attenuation by itself elicited a significant response
rez_struct.untuned_SUs_mat = cellfun( @(x) sum( x==0 ) , rez_struct.bandwidth_cell ) ;
rez_struct.untuned_SUs_mat( : , 2 ) = cellfun( @(x) size( x , 1) , rez_struct.bandwidth_cell ) ;

%% Calculate broadness of tuning using lifetime sparseness

% Lifetime sparseness
CHECK_EXCITATION = 1 ;
WINDOW_SIZE = 50 ;
responsive_SUs.lifetime_sparse = calculate_lifetimesparsness( responsive_SUs , PRE_STIM_MSEC , WINDOW_SIZE , protocol_struct , CHECK_EXCITATION) ;

% Summary lifetime sparseness across brain region and plot boxplot
rez_struct.lifetime_sparseness = cell(size( brain_regions , 2 ) , 1 ) ;
for kk = 1 : size( brain_regions , 2 )
    
    reg_units = responsive_SUs( responsive_SUs.acronym == brain_regions{ 1, kk } , : ) ;
    reg_units = reg_units( reg_units.is_excited == 1 , : ) ;
    rez_struct.lifetime_sparseness{kk,1} = reg_units.lifetime_sparse ;
    
end  

title_text = 'Lifetime Sparseness' ;
DATA_COL = 1 ;
create_bp_from_cells( rez_struct.lifetime_sparseness , brain_regions , my_colormap , title_text , DATA_COL  )

% KW for lifetime sparseness
[ lifetime_sparseness_data ,lifetime_sparseness_groups ] = create_grouped_vector( rez_struct.lifetime_sparseness ) ;
[p,tbl,stats] = kruskalwallis( lifetime_sparseness_data , lifetime_sparseness_groups )
[c,~,~,gnames] = multcompare(stats)


%% Clear variables

clear lifetime_sparseness_data lifetime_sparseness_groups p tbl stats c gnames
clc

%% Calculate population sparseness

% Create empty data structures
rez_struct.reg_pop_sparseness = cell( size( brain_regions , 2 ) , 1 ) ;

% Set tick values for population sparseness matrix presentation
X_TICK_VECTOR_FRA = 1 : 5 : length( protocol_struct.freqs ) ;
X_TICK_LABEL_FRA = { num2str( round( protocol_struct.freqs( 1 :5 :  end )' * 10 ) / 10 ) } ; 
Y_TICK_VECTOR_FRA = 1 : length ( protocol_struct.atten ) ;
Y_TICK_LABEL_FRA = num2str(  protocol_struct.atten( 1 : end )' ) ;
X_TICK_ANGLE_FRA = 90 ;

% Calculate and plot population sparseness
h_fig = figure() ;
for mm =1 : size( brain_regions , 2 )
    
    region_units = responsive_SUs( responsive_SUs.acronym == brain_regions{ 1 , mm } , : ) ;
    signif_units = region_units( region_units.is_excited == 1 , : ) ;
    reg_signif_FRA = cell2mat( reshape( signif_units.FRA_signif , [ 1 , 1 , size( signif_units , 1 ) ] ) ) ;
    reg_signif_FRA = sum( reg_signif_FRA , 3 ) ;
    reg_pop_sparseness = reg_signif_FRA ./ size( region_units , 1 ) ;
    rez_struct.reg_pop_sparseness{ mm , 1 } = reg_pop_sparseness ;
    
    h_ax = subplot( size( brain_regions ,2  ) , 1 , mm ) ;
    imagesc( h_ax , reg_pop_sparseness ) ;
    caxis( h_ax , [ 0 ,0.3 ] ) ; 

    set( h_ax , 'XTick' , X_TICK_VECTOR_FRA , 'XTickLabel', X_TICK_LABEL_FRA ) ;
    xtickangle( h_ax , X_TICK_ANGLE_FRA ) ;
    set( h_ax , 'YTick' , Y_TICK_VECTOR_FRA , 'YTickLabel' , Y_TICK_LABEL_FRA ) ;
    h_ax.FontSize = 14 ;
    
end

% Plot population sparseness box plot
pop_sparseness_bp = cellfun( @(x) x(:) , rez_struct.reg_pop_sparseness , 'UniformOutput' , false ) ;
title_text = 'Population Sparseness' ;
DATA_COL = 1 ;
create_bp_from_cells( pop_sparseness_bp , brain_regions , my_colormap , title_text , DATA_COL  )
 
% KW for tuning width
[ popsparse_data ,popsparse_groups ] = create_grouped_vector( pop_sparseness_bp ) ;
[p,tbl,stats] = kruskalwallis( popsparse_data , popsparse_groups )
[c,~,~,gnames] = multcompare(stats)

%% Clear variables

clear popsparse_data popsparse_groups p tbl stats c gnames X_TICK_ANGLE_FRA X_TICK_LABEL_FRA X_TICK_VECTOR_FRA
clear Y_TICK_LABEL_FRA Y_TICK_VECTOR_FRA
clc 

%% Latency to peak

% Calculate latency to peak
SMOOTHING_WINDOW = 9 ;
responsive_SUs.latency_to_peak = find_latency_to_peak( responsive_SUs , PRE_STIM_MSEC , SMOOTHING_WINDOW ) ;

% Summarize latency across regions and plot box plot
for kk = 1 : size( brain_regions , 2 ) 

    reg_units = responsive_SUs( responsive_SUs.acronym == brain_regions{ 1 , kk } , : ) ;
    reg_units = reg_units( reg_units.is_excited  == 1  , : ) ;
    rez_struct.lat_to_peak{ kk , 1 } = reg_units.latency_to_peak ;

end

title_text = 'Latency to PSTH peak' ;
DATA_COL = 1 ;
create_bp_from_cells( rez_struct.lat_to_peak , brain_regions , my_colormap , title_text , DATA_COL )

% KW for latency to lat to peak
[ lat_to_peak_data ,lat_to_peak_groups ] = create_grouped_vector( rez_struct.lat_to_peak ) ;
[p,tbl,stats] = kruskalwallis( lat_to_peak_data , lat_to_peak_groups )
[c,~,~,gnames] = multcompare(stats)

%% Clear variables

clear lat_to_peak_data lat_to_peak_groups p tbl stats c gnames
clc

%% Latency to min

% Calculate latency to min
SMOOTHING_WINDOW = 9 ;
responsive_SUs.latency_to_min = find_latency_to_min( responsive_SUs , PRE_STIM_MSEC , SMOOTHING_WINDOW ) ;

% Summary across brain regions and plot boxplot
for kk = 1 : size( brain_regions , 2 ) 

    reg_units = responsive_SUs( responsive_SUs.acronym == brain_regions{ 1 , kk } , : ) ;
    reg_units = reg_units( reg_units.is_inhibited  == 1  & reg_units.is_excited == 0, : ) ;
    rez_struct.lat_to_min{ kk , 1 } = reg_units.latency_to_min ;

end

title_text = 'Latency to PSTH min' ;
DATA_COL = 1 ;
create_bp_from_cells( rez_struct.lat_to_min , brain_regions , my_colormap , title_text , DATA_COL )

% KW for latency to lat to min
[ lat_to_min_data ,lat_to_min_groups ] = create_grouped_vector( rez_struct.lat_to_min ) ;
[p,tbl,stats] = kruskalwallis( lat_to_min_data , lat_to_min_groups )
[c,~,~,gnames] = multcompare(stats)

%% Clear variables

clear lat_to_min_data lat_to_min_groups p tbl stats c gnames
clc

%% Latency vs. sparseness

% Inspect the relationship between latency and lifetime sparseness by
% plotting a scatter, fitting a linear model and calculating their
% corresponding correlation coefficient 

excited_SUs = responsive_SUs( responsive_SUs.is_excited == 1 , : ) ;
h_fig = figure() ;
h_ax = axes('Parent' ,h_fig ) ;
for kk = 1 : size( brain_regions , 2 ) 
    
    reg_units = excited_SUs( excited_SUs.acronym == brain_regions{ 1 , kk } , : ) ;
    scatter( h_ax , reg_units.lifetime_sparse , reg_units.latency_to_peak  , 20 , my_colormap(kk,:) )  ;
    hold( h_ax ,'on' ) ;
    
end

mdl = fitlm(excited_SUs.lifetime_sparse , excited_SUs.latency_to_peak ) 
plot( h_ax , [0:0.1:1] ,[0:0.1:1].*mdl.Coefficients.Estimate(2) + mdl.Coefficients.Estimate(1), '--k' ) ;
% plot( h_ax , [0:0.1:1] ,[0:0.1:1].*mdl.Coefficients.Estimate(1) , '--k' ) ;
ylabel( h_ax , 'Latency to peak [ms]' ) ;
xlabel( h_ax , 'Lifetime Sparseness' ) ; 
h_ax.FontName = 'Arial' ;

[ r , p ] = corrcoef( excited_SUs.lifetime_sparse , excited_SUs.latency_to_peak ) ;

%% Clear variables

clear r p mdl excited_SUs
clc

%% Plot population PSTH averaged over sitmuli

SMOOTH_BIN = 9 ;
title_text = 'Population PSTH per stim - excited units' ;
rez_struct.exc_pop_PSTHs_perstim = plot_population_resp_perstim( responsive_SUs( responsive_SUs.is_excited == 1 , :  )  , brain_regions , PRE_STIM_MSEC, protocol_struct, SMOOTH_BIN , my_colormap , title_text ) ;
rez_struct.exc_pop_PSTHcompare_viattest = compare_PSTHs_via_ttest( rez_struct.exc_pop_PSTHs_perstim , brain_regions) ;

title_text = 'Population PSTH per stim - inhibited units' ;
rez_struct.inh_pop_PSTHs_perstim = plot_population_resp_perstim( responsive_SUs( (  responsive_SUs.is_excited == 0 & responsive_SUs.is_inhibited == 1  ) , :  )  , brain_regions , PRE_STIM_MSEC, protocol_struct, SMOOTH_BIN , my_colormap , title_text ) ;
rez_struct.inh_pop_PSTHcompare_viattest = compare_PSTHs_via_ttest( rez_struct.inh_pop_PSTHs_perstim , brain_regions ) ;

pval_thresh = 0.05 ./ ( size( rez_struct.exc_pop_PSTHcompare_viattest , 1 ) * size( rez_struct.exc_pop_PSTHcompare_viattest , 2 ) ) ;
[ ~ , time_to_differntiate_exc ]= max( rez_struct.exc_pop_PSTHcompare_viattest < pval_thresh , [] , 2 ) ;
[ ~ , time_to_differntiate_inh ]= max( rez_struct.inh_pop_PSTHcompare_viattest < pval_thresh , [] , 2 ) ;

%% Compare to baseline

TIME_DELAY = 150 ;

[ rez_struct.exc_pop_PSTHcomparetobaseline_viattest , rez_struct.exc_signif_times ] = find_time_to_baseline_ttest( rez_struct.exc_pop_PSTHs_perstim ,...
  brain_regions , PRE_STIM_MSEC , pval_thresh , TIME_DELAY ) ; 
[ rez_struct.inh_pop_PSTHcomparetobaseline_viattest  , rez_struct.inh_signif_times ] = find_time_to_baseline_ttest( rez_struct.inh_pop_PSTHs_perstim ,...
  brain_regions , PRE_STIM_MSEC , pval_thresh , TIME_DELAY ) ; 
%% PSTH_dynamics Matrices

% Re-bin PSTHs to larger temporal bins
OLD_BIN = 1 ;                                                                                                                       % PSTH bins are 1 msec long
NEW_BIN = 10 ;                                                                                                                     % New bin size for rasters       
SMOOTH_BIN = 3 ;
responsive_SUs.PSTH_binned = re_bin_PSTH( responsive_SUs , OLD_BIN , NEW_BIN , SMOOTH_BIN ) ;

% Plot PSTH heatmaps
WITH_INHIBITED = 0 ; 
PRE_STIM_BINS = 10 ;
TIME_FRAME = [ 0 , 600 ] ;
make_PSTH_heatmaps2( responsive_SUs , brain_regions , PRE_STIM_BINS , NEW_BIN , TIME_FRAME )

%% Create bins and calculate noise correlations


% Bin raster plots pre calculation of noise correlations
OLD_BIN = 1 ;                                                                                                                       % Raster bins are 1 msec long
NEW_BIN = 20 ;                                                                                                                     % New bin size for rasters                                                
responsive_SUs = re_bin_pure_tones( responsive_SUs , OLD_BIN , NEW_BIN ) ;
BINNED_RAST_COLUMN = size(responsive_SUs , 2 ) ;
% Calculate trial-to-trial variability/fluctuations
responsive_SUs.binned_fluct_mat =  calculate_fluct_mat( responsive_SUs , protocol_struct , BINNED_RAST_COLUMN ) ;

% Calculate noise correlations
brain_regions = { 'AUDp' , 'AUDv' , 'TeA' } ;
PRE_STIM_MSEC = 100 ;
FLUCT_MAT_COL = BINNED_RAST_COLUMN + 1 ; 
BINS_TO_CONSIDER = [ 6 : 15 ] ; 
IS_FOR_EXCITED = 0 ;
% Set number of iterations to calculate shuffled distribution
N_ITERS = 500 ;


begin_time = (BINS_TO_CONSIDER(1) - 1 ) .* NEW_BIN - PRE_STIM_MSEC ;
end_time = (BINS_TO_CONSIDER(end) - 1 ) .* NEW_BIN - PRE_STIM_MSEC ;
[rez_struct.all_dates_noise_corr_cell , rez_struct.total_noise_corr_cell ,...
 rez_struct.all_dates_shuff_noise_corr_cell , rez_struct.total_shuff_noise_corr_cell ] =  calc_noise_correlations2( responsive_SUs  ,...
 brain_regions , FLUCT_MAT_COL , BINS_TO_CONSIDER , N_ITERS , IS_FOR_EXCITED ) ;
%% Plot noise corr results within and across regions

% Within region noise correlations
title_text = [ 'Within region noise correlations , Bin Size ' num2str( NEW_BIN ) ' msec , Time window = ' num2str( begin_time ) '-' num2str( end_time ) ] ;
within_region = rez_struct.total_noise_corr_cell( 1 : ( size(rez_struct.total_noise_corr_cell , 1 ) +1 ) : size(rez_struct.total_noise_corr_cell , 1 )^2 ) ; 
within_region = within_region' ;
within_region_shuffle = rez_struct.total_shuff_noise_corr_cell( 1 : ( size(rez_struct.total_noise_corr_cell , 1 ) +1 ) : size(rez_struct.total_noise_corr_cell , 1 )^2 ) ; 
within_region_shuffle = within_region_shuffle' ;

create_bp_from_cells ( within_region , brain_regions , my_colormap , title_text )
YLIM = [0 , 0.5 ] ; 
rez_struct.within_reg_sig_corr_probs = plot_noisecorrelation_distribution2( within_region ,...
    within_region_shuffle , brain_regions , title_text , my_colormap ,YLIM ) ;

% Across region noise correlations
cross_reg_combos = {'AUDp x AUDv' , 'AUDp x TeA' , 'AUDv x TeA' } ;

title_text = [ 'Across region noise correlations , Bin Size ' num2str( NEW_BIN ) ' msec , Time window = ' num2str( begin_time ) '-' num2str( end_time ) ] ;
across_region_indices = [ 4 ,7, 8 ] ;
across_regions = rez_struct.total_noise_corr_cell( across_region_indices ) ; 
across_regions = across_regions' ;
across_regions_shuffle = rez_struct.total_shuff_noise_corr_cell( across_region_indices ) ; 
across_regions_shuffle =across_regions_shuffle' ;
create_bp_from_cells (across_regions , cross_reg_combos , my_colormap , title_text )
YLIM = [0 , 0.45 ] ; 
rez_struct.across_reg_sig_corr_probs = plot_noisecorrelation_distribution2( across_regions ,...
    across_regions_shuffle , cross_reg_combos , title_text , my_colormap ,YLIM ) ;

%% Calculate within region Signal Correlations

N_ITERS = 500 ;
[ rez_struct.within_reg_sig_corrs , rez_struct.within_reg_sig_corr_bydate ] = calculate_sig_corrs3(responsive_SUs, brain_regions , N_ITERS ) ;
title_text = [ 'Within region Signal Correlations - by date' ] ;
create_bp_from_cells( rez_struct.within_reg_sig_corrs(:,1) , brain_regions , my_colormap , title_text ) ;

% Stats A1 vs. TeA 
[ p , h ] = ranksum( rez_struct.within_reg_sig_corrs{1,1} , rez_struct.within_reg_sig_corrs{3,1} ) 
% Stats A2 vs. TeA 
[ p , h ] = ranksum( rez_struct.within_reg_sig_corrs{2,1} , rez_struct.within_reg_sig_corrs{3,1} ) 
% Stats A1 vs. A2 
[ p , h ] = ranksum( rez_struct.within_reg_sig_corrs{1,1} , rez_struct.within_reg_sig_corrs{2,1} ) 

title_text = [ 'Within region Signal Correlations' ] ; 
rez_struct.within_reg_sig_corr_probs = plot_correlation_distribution( rez_struct.within_reg_sig_corrs , brain_regions , title_text , my_colormap) ; 

%% Clear variables

clear p h 
clc 

%% Calculate across region Signal Correlations

reg_combos = {'AUDpxAUDv' , 'AUDpxTeA' , 'AUDvxTeA' } ;
[ rez_struct.across_reg_sig_corrs , rez_struct.across_reg_sig_corrs_bydate ] = calculate_crossreg_sig_corrs3( responsive_SUs, brain_regions , N_ITERS ) ;
title_text = [ 'Cross regions region Signal Correlations' ] ;
create_bp_from_cells( rez_struct.across_reg_sig_corrs(:,1) , reg_combos , my_colormap , title_text ) ;
% Stats A1xA2 vs. A1xTeA 
[ p , h ] = ranksum( rez_struct.across_reg_sig_corrs{1,2} , rez_struct.across_reg_sig_corrs{2,2} ) 
% Stats A1xTeA vs. A2xTeA 
[ p , h ] = ranksum( rez_struct.across_reg_sig_corrs{2,2} , rez_struct.across_reg_sig_corrs{3,2} ) 
% Stats A1xA2 vs. A2xTeA 
[ p , h ] = ranksum( rez_struct.across_reg_sig_corrs{1,2} , rez_struct.across_reg_sig_corrs{3,2} ) 

title_text = [ 'Across region Signal Correlations' ] ; 
rez_struct.across_reg_sig_corr_probs = plot_correlation_distribution( rez_struct.across_reg_sig_corrs , reg_combos , title_text , my_colormap ) ; 

%% Clear variables

clear p h 
clc

%% Save rez_struct

save( 'E:\Final_datasets_forthesis\Awake\Naives\NV_awake_rez_struct.mat' , '-v7.3' , 'rez_struct' ) ;
