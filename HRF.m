%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% BRAIN RESPONSE FUNCTIONS AND NEUROVASCULAR COUPLING IN TYPE 2 DIABETES:
% INSIGHTS FROM FMRI
% 
%                       Catarina Guerra | 2015240209
%                               December 2020
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
close all
clc

load('ROI_betas.mat');
load('ROI_SEs.mat');

%% CONSTANTS:

n_points = 8;                             % data points
time = (0:n_points-1)* 2.5;               % volumes per subject * TR (2.5s)
n_subjects = 141;
n_T2DM = 64; 
n_CNT = 77;
n_rois = 22;

ROI_regions = {'L IPL BA40', 'L Insula BA13', 'L Precuneus BA7', 'R IFG BA9', 'R MFG BA8', 'R MFG BA46', 'R MT BA19', ... 
               'R SFG BA6', 'R SPL BA7', 'R V2 BA18', 'L AC BA32', 'L CG BA31', 'L PC BA30', 'L PrcL BA5', 'L PrhG BA36', ...
               'R CG BA24', 'R Insula BA13', 'R PC BA23', 'R PrecG BA4', 'R PrecG BA43', 'R PrimSens BA1', 'R STG BA39'};


%% RAW DATA EXTRACTION:

% In this section, we retrieve the raw data (betas and standard errors) 
% from both sides for each group and condition

T2DM_Thr_betas = [];
T2DM_Sub_betas = [];

T2DM_Thr_SE = [];
T2DM_Sub_SE = [];

CNT_Thr_betas = [];
CNT_Sub_betas = [];

CNT_Thr_SE = [];
CNT_Sub_SE = [];

% 2065 - where left hemisphere activations from the Threshold condition in
% CNT subjects begin
% 2073 - where left hemisphere activations from the Submaximum condition in
% CNT subjects begin
% 2049 - where right hemisphere activations from the Threshold condition in
% CNT subjects begin
% 2057 - where right hemisphere activations from the Submaximum condition 
% in CNT subjects begin
% However, these two conditions are interleaved, hence the (c-1)*32 
% expression.

for c = 1:n_CNT
    CNT_Thr_betas = [CNT_Thr_betas; ROIbetas(2065+((c-1)*32): 2072+((c-1)*32),:); ROIbetas(2049+((c-1)*32): 2056+((c-1)*32),:)];
    CNT_Sub_betas = [CNT_Sub_betas; ROIbetas(2073+((c-1)*32): 2080+((c-1)*32),:); ROIbetas(2057+((c-1)*32): 2064+((c-1)*32),:)];
    CNT_Thr_SE = [CNT_Thr_SE; ROISEs(2065+((c-1)*32): 2072+((c-1)*32),:); ROISEs(2049+((c-1)*32): 2056+((c-1)*32),:)];
    CNT_Sub_SE = [CNT_Sub_SE; ROISEs(2073+((c-1)*32): 2080+((c-1)*32),:); ROISEs(2057+((c-1)*32): 2064+((c-1)*32),:)];
        
end

% 17 - where left hemisphere activations from the Threshold condition in
% T2DM subjects begin
% 25 - where left hemisphere activations from the Submaximum condition in
% T2DM subjects begin
% 1 - where right hemisphere activations from the Threshold condition in
% T2DM subjects begin
% 9 - where right hemisphere activations from the Submaximum condition in 
% T2DM subjects begin
% Again, these two conditions are interleaved, hence the (c-1)*32 
% expression.

for t = 1:n_T2DM
    T2DM_Thr_betas = [T2DM_Thr_betas; ROIbetas(17+((t-1)*32): 24+((t-1)*32),:); ROIbetas(1+((t-1)*32): 8+((t-1)*32),:)];
    T2DM_Sub_betas = [T2DM_Sub_betas; ROIbetas(25+((t-1)*32): 32+((t-1)*32),:); ROIbetas(9+((t-1)*32): 16+((t-1)*32),:)];
    T2DM_Thr_SE = [T2DM_Thr_SE; ROISEs(17+((t-1)*32): 24+((t-1)*32),:); ROISEs(1+((t-1)*32): 8+((t-1)*32),:)];
    T2DM_Sub_SE = [T2DM_Sub_SE; ROISEs(25+((t-1)*32): 32+((t-1)*32),:); ROISEs(9+((t-1)*32): 16+((t-1)*32),:)];
end


%% RAW DATA TABLES:

% In this section, we create raw data tables regarding each condition and 
% group 

% Data information column for each group
CNT_info = cell(length(CNT_Thr_betas),1);
T2DM_info = cell(length(T2DM_Thr_betas),1);

for c = 1:n_CNT
    for side = 1:2
        for dp = 1:n_points
            if side == 1
                CNT_info{dp+(8*(side-1))+(16*(c-1))} = {sprintf('%s %d %s %d %s', 'Subject', c, 'datapoint', dp, 'Left')};
            else
                CNT_info{dp+(8*(side-1))+(16*(c-1))} = {sprintf('%s %d %s %d %s', 'Subject', c, 'datapoint', dp, 'Right')};
            end
        end
    end
end

for t = 1:n_T2DM
    for side = 1:2
        for dp = 1:n_points
            if side == 1
                T2DM_info{dp+(8*(side-1))+(16*(t-1))} = {sprintf('%s %d %s %d %s', 'Subject', t, 'datapoint', dp, 'Left')};
            else
                T2DM_info{dp+(8*(side-1))+(16*(t-1))} = {sprintf('%s %d %s %d %s', 'Subject', t, 'datapoint', dp, 'Right')};
            end
        end
    end
end

% Cell with the ROI's names which will be used as a header
headings_ROIs = cell(1,length(ROI_regions));

% Space removal
for r = 1:length(ROI_regions)
    headings_ROIs{r} = char(strrep(strcat(ROI_regions(r)),' ','_'));
end

% Every raw data in each ROI per condition and group
raw_data_CNT_Thr_betas = array2table(CNT_Thr_betas(:,2:end),'VariableNames',headings_ROIs);
raw_data_CNT_Thr_SE = array2table(CNT_Thr_SE(:,2:end),'VariableNames',headings_ROIs);
raw_data_CNT_Sub_betas = array2table(CNT_Sub_betas(:,2:end),'VariableNames',headings_ROIs);
raw_data_CNT_Sub_SE = array2table(CNT_Sub_SE(:,2:end),'VariableNames',headings_ROIs);

raw_data_T2DM_Thr_betas = array2table(T2DM_Thr_betas(:,2:end),'VariableNames',headings_ROIs);
raw_data_T2DM_Thr_SE = array2table(T2DM_Thr_SE(:,2:end),'VariableNames',headings_ROIs);
raw_data_T2DM_Sub_betas = array2table(T2DM_Sub_betas(:,2:end),'VariableNames',headings_ROIs);
raw_data_T2DM_Sub_SE = array2table(T2DM_Sub_SE(:,2:end),'VariableNames',headings_ROIs);

headings_CNT = array2table(CNT_info,'VariableNames', {'Subject_datapoint_side'});           
headings_T2DM = array2table(T2DM_info,'VariableNames', {'Subject_datapoint_side'});           

% Raw data tables per condition and group
table_raw_CNT_Thr_betas = [headings_CNT raw_data_CNT_Thr_betas];
table_raw_CNT_Thr_SE = [headings_CNT raw_data_CNT_Thr_SE];
table_raw_CNT_Sub_betas = [headings_CNT raw_data_CNT_Sub_betas];
table_raw_CNT_Sub_SE = [headings_CNT raw_data_CNT_Sub_SE];

table_raw_T2DM_Thr_betas = [headings_T2DM raw_data_T2DM_Thr_betas];
table_raw_T2DM_Thr_SE = [headings_T2DM raw_data_T2DM_Thr_SE];
table_raw_T2DM_Sub_betas = [headings_T2DM raw_data_T2DM_Sub_betas];
table_raw_T2DM_Sub_SE = [headings_T2DM raw_data_T2DM_Sub_SE];


%% AVERAGE AND MEDIAN HRF CURVES AND ITS PARAMETERS:

% In this section, we estimate the average and median HRF curves in each 
% condition and group per Regions of Interest (ROIs), as well as its
% parameters of interest.

thr_both_sides = [];
sub_both_sides = [];

mean_per_subject_thr = [];
mean_per_subject_sub = [];

mean_per_subject_thr_CNT = [];
mean_per_subject_sub_CNT = [];
mean_per_subject_thr_T2DM = [];
mean_per_subject_sub_T2DM = [];

total_avg_thr_CNT = [];
total_avg_thr_T2DM = [];
total_std_thr_CNT = [];
total_std_thr_T2DM = [];

total_avg_sub_CNT = [];
total_avg_sub_T2DM = [];
total_std_sub_CNT = [];
total_std_sub_T2DM = [];

total_median_thr_CNT = [];
total_median_thr_T2DM = [];
total_iqr_thr_CNT = [];
total_iqr_thr_T2DM = [];

total_median_sub_CNT = [];
total_median_sub_T2DM = [];
total_iqr_sub_CNT = [];
total_iqr_sub_T2DM = [];

total_average_HRF_parameters_thr_CNT = [];
total_average_HRF_parameters_thr_T2DM = [];
total_median_HRF_parameters_thr_CNT = [];
total_median_HRF_parameters_thr_T2DM = [];

total_average_HRF_parameters_sub_CNT = [];
total_average_HRF_parameters_sub_T2DM = [];
total_median_HRF_parameters_sub_CNT = [];
total_median_HRF_parameters_sub_T2DM = [];


for r=2:n_rois+1
    for i=1:n_subjects
        for j=1:2
            sided_thr = cell2mat(ROIbetas(1+(n_points*2)*(j-1)+(n_points*4)*(i-1):n_points+(n_points*2)*(j-1)+(n_points*4)*(i-1),r));    % each column is a side in each subject
            sided_sub = cell2mat(ROIbetas((n_points+1)+(n_points*2)*(j-1)+(n_points*4)*(i-1):(n_points*2)+(n_points*2)*(j-1)+(n_points*4)*(i-1),r));
            thr_both_sides = [thr_both_sides sided_thr]; 
            sub_both_sides = [sub_both_sides sided_sub]; 
        end

        mean_thr = mean(thr_both_sides,2);                                                      % mean over rows (2)
        mean_sub = mean(sub_both_sides,2);                                               
        mean_per_subject_thr = [mean_per_subject_thr mean_thr];                                 % vector with the mean of both sides for each Thr subject
        mean_per_subject_sub = [mean_per_subject_sub mean_sub];

        if i<65                                                                                 % until subject 65 --> T2DM
            mean_per_subject_thr_T2DM = [mean_per_subject_thr_T2DM mean_per_subject_thr];       % every mean value of beta per T2DM subject in the Thr condition for each ROI 
            mean_per_subject_sub_T2DM = [mean_per_subject_sub_T2DM mean_per_subject_sub];  
        else                                                                                    % after subject 65 --> CNT
            mean_per_subject_thr_CNT = [mean_per_subject_thr_CNT mean_per_subject_thr];         % every mean value of beta per CNT subject in the Thr condition for each ROI 
            mean_per_subject_sub_CNT = [mean_per_subject_sub_CNT mean_per_subject_sub];
        end
    
        % Delete to not duplicate or interfere with data
        thr_both_sides = [];
        sub_both_sides = [];
        mean_per_subject_thr = [];
        mean_per_subject_sub = [];
    end  
    
    
    % --------------------- THRESHOLD CONDITION -------------------------
    
    % Estimating the beta values' Average, Standard deviation, Median and
    % Interquartile range per datapoint in each group subject of the Thr
    % condition
    avg_thr_CNT = mean(mean_per_subject_thr_CNT,2);                        
    avg_thr_T2DM = mean(mean_per_subject_thr_T2DM,2);                     
    std_thr_CNT = std(mean_per_subject_thr_CNT,0,2);                      
    std_thr_T2DM = std(mean_per_subject_thr_T2DM,0,2);                    
    
    median_thr_CNT = median(mean_per_subject_thr_CNT,2);                  % median of all CNT individuals
    median_thr_T2DM = median(mean_per_subject_thr_T2DM,2);                % median of all T2DM individuals
    iqr_thr_CNT = iqr(mean_per_subject_thr_CNT,2);                        % inter-quartile range per rows (2) (for each CNT individual) in each timepoint
    iqr_thr_T2DM = iqr(mean_per_subject_thr_T2DM,2);                      % inter-quartile range per rows (2) (for each T2DM individual) in each timepoint

    
    total_avg_thr_CNT = [total_avg_thr_CNT avg_thr_CNT];                  % stores every average beta values per datapoint of CNT subjects in the Thr condition
    total_avg_thr_T2DM = [total_avg_thr_T2DM avg_thr_T2DM];               % stores every average beta values per datapoint of T2DM subjects in the Thr condition
    total_std_thr_CNT = [total_std_thr_CNT std_thr_CNT];                  % stores every standard deviation of the beta values over each datapoint of CNT subjects in the Thr condition
    total_std_thr_T2DM = [total_std_thr_T2DM std_thr_T2DM];               % stores every standard deviation of the beta values over each datapoint of T2DM subjects in the Thr condition
    
    total_median_thr_CNT = [total_median_thr_CNT median_thr_CNT];         % stores every median beta values per datapoint of CNT subjects in the Thr condition
    total_median_thr_T2DM = [total_median_thr_T2DM median_thr_T2DM];      % stores every median beta values per datapoint of T2DM subjects in the Thr condition
    total_iqr_thr_CNT = [total_iqr_thr_CNT iqr_thr_CNT];                  % stores every interquartile range of the beta values over each datapoint of CNT subjects in the Thr condition
    total_iqr_thr_T2DM = [total_iqr_thr_T2DM iqr_thr_T2DM];               % stores every interquartile range of the beta values over each datapoint of T2DM subjects in the Thr condition
        
    
        % «««««««««««««««««««««««« PARAMETERS »»»»»»»»»»»»»»»»»»»»»»»»»»»
        
        % Estimating the HRF parameters of interest per group and condition 
        % in each ROI based on the average / median HRF 
        
        % **************************** CNT ******************************

        [peak,volume] = max(avg_thr_CNT(2:5,1));                                        % peak amplitude: maximum within the sectioned time range
        peak_latency = volume * 2.5;                                                    % peak latency: time to peak (volume where the maximum takes place * TR) (we do not need to subtract 1 to the volume due to MATLAB indexing since the section we considered to estimate the peak started at the second datapoint, and thus, that datapoint acquires index 1 in this new array)
        slope_to_peak = minus(peak,avg_thr_CNT(1))/peak_latency;                        % slope to peak: (peak-1st point)/latency
        total_area = trapz((0:n_points-1)*2.5,avg_thr_CNT - min(avg_thr_CNT));          % area under the curve (trapezoidal method, from 0 to 17.5s)
        pos_area = positive_area((0:n_points-1)*2.5,avg_thr_CNT);                       % area of the positive sections of the HRF curve (from 0 to 15s)
        first_pos_area = positive_area((0:2)*2.5,avg_thr_CNT(1:3,1));                   % area of the first positive section of the HRF curve (from 0 to 5s)
        second_pos_area = positive_area((2:4)*2.5,avg_thr_CNT(3:5,1));                  % area of the second section of the HRF curve (from 5 to 10s)
        third_pos_area = positive_area((4:7)*2.5,avg_thr_CNT(5:8,1));                   % area of the third section of the HRF curve (from 10 to 17.5s)
        neg_area = negative_area((0:n_points-1)*2.5,avg_thr_CNT);                       % area of the negative sections of the HRF curve (from 0 to 17.5s)
        initial_dip_area = negative_area((0:volume)*2.5,avg_thr_CNT(1:(volume+1),1));   % area of the initial dip section of the HRF curve
        undershoot_area = negative_area((volume:7)*2.5,avg_thr_CNT((volume+1):8,1));    % area of the undershoot section of the HRF curve
        HRF_parameters = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];
        total_average_HRF_parameters_thr_CNT = [total_average_HRF_parameters_thr_CNT HRF_parameters];
        
        [peak,volume] = max(median_thr_CNT(2:5,1));
        peak_latency = volume * 2.5;
        slope_to_peak = minus(peak,median_thr_CNT(1))/peak_latency;
        total_area = trapz((0:n_points-1)*2.5,median_thr_CNT - min(median_thr_CNT));
        pos_area = positive_area((0:n_points-1)*2.5,median_thr_CNT);
        first_pos_area = positive_area((0:2)*2.5,median_thr_CNT(1:3,1));
        second_pos_area = positive_area((2:4)*2.5,median_thr_CNT(3:5,1));
        third_pos_area = positive_area((4:7)*2.5,median_thr_CNT(5:8,1));
        neg_area = negative_area((0:n_points-1)*2.5,median_thr_CNT);
        initial_dip_area = negative_area((0:volume)*2.5,median_thr_CNT(1:(volume+1),1));
        undershoot_area = negative_area((volume:7)*2.5,median_thr_CNT((volume+1):8,1));
        HRF_parameters = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];
        total_median_HRF_parameters_thr_CNT = [total_median_HRF_parameters_thr_CNT HRF_parameters];
        
        
        % *************************** T2DM ******************************
        
        [peak,volume] = max(avg_thr_T2DM(2:5,1));                                       
        peak_latency = volume * 2.5;                                                
        slope_to_peak = minus(peak,avg_thr_T2DM(1))/peak_latency;                       
        total_area = trapz((0:n_points-1)*2.5,avg_thr_T2DM - min(avg_thr_T2DM));        
        pos_area = positive_area((0:n_points-1)*2.5,avg_thr_T2DM);                      
        first_pos_area = positive_area((0:2)*2.5,avg_thr_T2DM(1:3,1));                  
        second_pos_area = positive_area((2:4)*2.5,avg_thr_T2DM(3:5,1));                 
        third_pos_area = positive_area((4:7)*2.5,avg_thr_T2DM(5:8,1));                  
        neg_area = negative_area((0:n_points-1)*2.5,avg_thr_T2DM);                      
        initial_dip_area = negative_area((0:volume)*2.5,avg_thr_T2DM(1:(volume+1),1));  
        undershoot_area = negative_area((volume:7)*2.5,avg_thr_T2DM((volume+1):8,1));   
        HRF_parameters = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];
        total_average_HRF_parameters_thr_T2DM = [total_average_HRF_parameters_thr_T2DM HRF_parameters];
        
        [peak,volume] = max(median_thr_T2DM(2:5,1));
        peak_latency = volume * 2.5;
        slope_to_peak = minus(peak,median_thr_T2DM(1))/peak_latency;
        total_area = trapz((0:n_points-1)*2.5,median_thr_T2DM - min(median_thr_T2DM));
        pos_area = positive_area((0:n_points-1)*2.5,median_thr_T2DM);
        first_pos_area = positive_area((0:2)*2.5,median_thr_T2DM(1:3,1));
        second_pos_area = positive_area((2:4)*2.5,median_thr_T2DM(3:5,1));
        third_pos_area = positive_area((4:7)*2.5,median_thr_T2DM(5:8,1));
        neg_area = negative_area((0:n_points-1)*2.5,median_thr_T2DM);
        initial_dip_area = negative_area((0:volume)*2.5,median_thr_T2DM(1:(volume+1),1));
        undershoot_area = negative_area((volume:7)*2.5,median_thr_T2DM((volume+1):8,1));
        HRF_parameters = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];
        total_median_HRF_parameters_thr_T2DM = [total_median_HRF_parameters_thr_T2DM HRF_parameters];
       
        
    % ---------------------- SUBMAXIMUM CONDITION ------------------------
   
    % Estimating the beta values' Average, Standard deviation, Median and
    % Interquartile range per datapoint in each group subject of the Sub
    % condition
    avg_sub_CNT = mean(mean_per_subject_sub_CNT,2);                      
    avg_sub_T2DM = mean(mean_per_subject_sub_T2DM,2);                    
    std_sub_CNT = std(mean_per_subject_sub_CNT,0,2);
    std_sub_T2DM = std(mean_per_subject_sub_T2DM,0,2);                    

    median_sub_CNT = median(mean_per_subject_sub_CNT,2);                  
    median_sub_T2DM = median(mean_per_subject_sub_T2DM,2);                
    iqr_sub_CNT = iqr(mean_per_subject_sub_CNT,2);                        
    iqr_sub_T2DM = iqr(mean_per_subject_sub_T2DM,2);                      

    
    total_avg_sub_CNT = [total_avg_sub_CNT avg_sub_CNT];         
    total_avg_sub_T2DM = [total_avg_sub_T2DM avg_sub_T2DM];      
    total_std_sub_CNT = [total_std_sub_CNT std_sub_CNT];            
    total_std_sub_T2DM = [total_std_sub_T2DM std_sub_T2DM];           
    
    total_median_sub_CNT = [total_median_sub_CNT median_sub_CNT];         
    total_median_sub_T2DM = [total_median_sub_T2DM median_sub_T2DM];      
    total_iqr_sub_CNT = [total_iqr_sub_CNT iqr_sub_CNT];                  
    total_iqr_sub_T2DM = [total_iqr_sub_T2DM iqr_sub_T2DM];               
       
    
        % «««««««««««««««««««««««« PARAMETERS »»»»»»»»»»»»»»»»»»»»»»»»»»»

        % Estimating the HRF parameters of interest per group and condition 
        % in each ROI based on the average / median HRF 
        
        % **************************** CNT ******************************
        
        [peak,volume] = max(avg_sub_CNT(2:5,1));
        peak_latency = (volume-1) * 2.5;
        slope_to_peak = minus(peak,avg_sub_CNT(1))/peak_latency;
        total_area = trapz((0:n_points-1)*2.5,avg_sub_CNT - min(avg_sub_CNT));
        pos_area = positive_area((0:n_points-1)*2.5,avg_sub_CNT);
        first_pos_area = positive_area((0:2)*2.5,avg_sub_CNT(1:3,1));
        second_pos_area = positive_area((2:4)*2.5,avg_sub_CNT(3:5,1));
        third_pos_area = positive_area((4:7)*2.5,avg_sub_CNT(5:8,1));
        neg_area = negative_area((0:n_points-1)*2.5,avg_sub_CNT);
        initial_dip_area = negative_area((0:volume)*2.5,avg_sub_CNT(1:(volume+1),1));
        undershoot_area = negative_area((volume:7)*2.5,avg_sub_CNT((volume+1):8,1));
        HRF_parameters = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];
        total_average_HRF_parameters_sub_CNT = [total_average_HRF_parameters_sub_CNT HRF_parameters];
        
        [peak,volume] = max(median_sub_CNT(2:5,1));
        peak_latency = (volume-1) * 2.5;
        slope_to_peak = minus(peak,median_sub_CNT(1))/peak_latency;
        total_area = trapz((0:n_points-1)*2.5,median_sub_CNT - min(median_sub_CNT));
        pos_area = positive_area((0:n_points-1)*2.5,median_sub_CNT);
        first_pos_area = positive_area((0:2)*2.5,median_sub_CNT(1:3,1));
        second_pos_area = positive_area((2:4)*2.5,median_sub_CNT(3:5,1));
        third_pos_area = positive_area((4:7)*2.5,median_sub_CNT(5:8,1));
        neg_area = negative_area((0:n_points-1)*2.5,median_sub_CNT);
        initial_dip_area = negative_area((0:volume)*2.5,median_sub_CNT(1:(volume+1),1));
        undershoot_area = negative_area((volume:7)*2.5,median_sub_CNT((volume+1):8,1));
        HRF_parameters = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];
        total_median_HRF_parameters_sub_CNT = [total_median_HRF_parameters_sub_CNT HRF_parameters];
    
        
        % *************************** T2DM ******************************

        [peak,volume] = max(avg_sub_T2DM(2:5,1));
        peak_latency = (volume-1) * 2.5;
        slope_to_peak = minus(peak,avg_sub_T2DM(1))/peak_latency;
        total_area = trapz((0:n_points-1)*2.5,avg_sub_T2DM - min(avg_sub_T2DM));
        pos_area = positive_area((0:n_points-1)*2.5,avg_sub_T2DM);
        first_pos_area = positive_area((0:2)*2.5,avg_sub_T2DM(1:3,1));
        second_pos_area = positive_area((2:4)*2.5,avg_sub_T2DM(3:5,1));
        third_pos_area = positive_area((4:7)*2.5,avg_sub_T2DM(5:8,1));
        neg_area = negative_area((0:n_points-1)*2.5,avg_sub_T2DM);
        initial_dip_area = negative_area((0:volume)*2.5,avg_sub_T2DM(1:(volume+1),1));
        undershoot_area = negative_area((volume:7)*2.5,avg_sub_T2DM((volume+1):8,1));
        HRF_parameters = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];
        total_average_HRF_parameters_sub_T2DM = [total_average_HRF_parameters_sub_T2DM HRF_parameters];
        
        [peak,volume] = max(median_sub_T2DM(2:5,1));
        peak_latency = (volume-1) * 2.5;
        slope_to_peak = minus(peak,median_sub_T2DM(1))/peak_latency;
        total_area = trapz((0:n_points-1)*2.5,median_sub_T2DM - min(median_sub_T2DM));
        pos_area = positive_area((0:n_points-1)*2.5,median_sub_T2DM);
        first_pos_area = positive_area((0:2)*2.5,median_sub_T2DM(1:3,1));
        second_pos_area = positive_area((2:4)*2.5,median_sub_T2DM(3:5,1));
        third_pos_area = positive_area((4:7)*2.5,median_sub_T2DM(5:8,1));
        neg_area = negative_area((0:n_points-1)*2.5,median_sub_T2DM);
        initial_dip_area = negative_area((0:volume)*2.5,median_sub_T2DM(1:(volume+1),1));
        undershoot_area = negative_area((volume:7)*2.5,median_sub_T2DM((volume+1):8,1));
        HRF_parameters = [peak; peak_latency; slope_to_peak; total_area; pos_area; first_pos_area; second_pos_area; third_pos_area; neg_area; initial_dip_area; undershoot_area];
        total_median_HRF_parameters_sub_T2DM = [total_median_HRF_parameters_sub_T2DM HRF_parameters];
    
        
    % Delete to not duplicate or interfere with data
    mean_per_subject_thr_CNT = [];
    mean_per_subject_thr_T2DM = [];
    mean_per_subject_sub_CNT = [];
    mean_per_subject_sub_T2DM = [];
end


%% PLOTTING THE AVERAGE AND MEDIAN HRF CURVES: 

% In this section, we plot each average and median HRF curve per group and
% condition of each ROI, separating positive from negative signal change
% regions (regions where the HRF curve increases after an activation, 
% and regions where the HRF curve decreases after an activation,
% respectively).

% Plotting the average HRF according to each set of ROIs
figure('Name', sprintf('Average HRF - Positive signal change ROIs'))
for a=1:10                                                                  % 10 positive signal change regions of interest
    subplot(2,5,a);
    shadedErrorBar(time, total_avg_thr_CNT(:,a), (total_std_thr_CNT(:,a))','lineprops',{'-b','LineWidth', 2},'patchSaturation',0.05)
    shadedErrorBar(time, total_avg_thr_T2DM(:,a), (total_std_thr_T2DM(:,a))','lineprops',{'-r','LineWidth', 2},'patchSaturation',0.05)
    shadedErrorBar(time, total_avg_sub_CNT(:,a), (total_std_sub_CNT(:,a))','lineprops',{'--b','LineWidth', 1},'patchSaturation',0.025)
    shadedErrorBar(time, total_avg_sub_T2DM(:,a), (total_std_sub_T2DM(:,a))','lineprops',{'--r','LineWidth', 1},'patchSaturation',0.025)
    title(ROI_regions{a});
    xlabel('Time (s)')
    ylabel('Beta values')    
end
legend('CNT Thr','T2DM Thr','CNT Sub','T2DM Sub');

figure('Name', sprintf('Average HRF - Negative signal change ROIs'))
for d=1:12                                                                  % 12 negative signal change regions of interest
    subplot(3,4,d);
    shadedErrorBar(time, total_avg_thr_CNT(:,10+d), (total_std_thr_CNT(:,10+d))','lineprops',{'-b','LineWidth', 2},'patchSaturation',0.05)
    shadedErrorBar(time, total_avg_thr_T2DM(:,10+d), (total_std_thr_T2DM(:,10+d))','lineprops',{'-r','LineWidth', 2},'patchSaturation',0.05)    
    shadedErrorBar(time, total_avg_sub_CNT(:,10+d), (total_std_sub_CNT(:,10+d))','lineprops',{'--b','LineWidth', 1},'patchSaturation',0.025)
    shadedErrorBar(time, total_avg_sub_T2DM(:,10+d), (total_std_sub_T2DM(:,10+d))','lineprops',{'--r','LineWidth', 1},'patchSaturation',0.025)    
    title(ROI_regions{10+d});
    xlabel('Time (s)')
    ylabel('Beta values')  
end
legend('CNT Thr','T2DM Thr','CNT Sub','T2DM Sub');


% Plotting the median HRF according to each set of ROIs
figure('Name', sprintf('Median HRF - Positive signal change ROIs'))
for a=1:10                                                       
    subplot(2,5,a);
    shadedErrorBar(time, total_median_thr_CNT(:,a), (total_iqr_thr_CNT(:,a))','lineprops',{'-b','LineWidth', 2},'patchSaturation',0.05)
    shadedErrorBar(time, total_median_thr_T2DM(:,a), (total_iqr_thr_T2DM(:,a))','lineprops',{'-r','LineWidth', 2},'patchSaturation',0.05)    
    shadedErrorBar(time, total_median_sub_CNT(:,a), (total_iqr_sub_CNT(:,a))','lineprops',{'--b','LineWidth', 1},'patchSaturation',0.025)
    shadedErrorBar(time, total_median_sub_T2DM(:,a), (total_iqr_sub_T2DM(:,a))','lineprops',{'--r','LineWidth', 1},'patchSaturation',0.025)    
    title(ROI_regions{a});
    xlabel('Time (s)')
    ylabel('Beta values')    
end
legend('CNT Thr','T2DM Thr','CNT Sub','T2DM Sub');

figure('Name', sprintf('Median HRF - Negative signal change ROIs'))
for d=1:12                                                       
    subplot(3,4,d);
    shadedErrorBar(time, total_median_thr_CNT(:,10+d), (total_iqr_thr_CNT(:,10+d))','lineprops',{'-b','LineWidth', 2},'patchSaturation',0.05)
    shadedErrorBar(time, total_median_thr_T2DM(:,10+d), (total_iqr_thr_T2DM(:,10+d))','lineprops',{'-r','LineWidth', 2},'patchSaturation',0.05)
    shadedErrorBar(time, total_median_sub_CNT(:,10+d), (total_iqr_sub_CNT(:,10+d))','lineprops',{'--b','LineWidth', 1},'patchSaturation',0.025)
    shadedErrorBar(time, total_median_sub_T2DM(:,10+d), (total_iqr_sub_T2DM(:,10+d))','lineprops',{'--r','LineWidth', 1},'patchSaturation',0.025)
    title(ROI_regions{10+d});
    xlabel('Time (s)')
    ylabel('Beta values')  
end
legend('CNT Thr','T2DM Thr','CNT Sub','T2DM Sub');


%% GRAND AVERAGE / MEDIAN ANALYSIS - PARAMETERS AND PLOTS:

% In this section, we do a overall average / median of the HRF curves seen
% in each ROI and its parameters of interest for each condition and group 


% Overall average and standard deviation / median and interquartile range
% of the parameters of the average / median HRF of the positive signal 
% change ROIs
avg_psc_thr_CNT_parameters = mean(total_average_HRF_parameters_thr_CNT(:,1:10), 2);
avg_psc_thr_T2DM_parameters = mean(total_average_HRF_parameters_thr_T2DM(:,1:10), 2);
std_psc_thr_CNT_parameters = std(total_average_HRF_parameters_thr_CNT(:,1:10), 0, 2);
std_psc_thr_T2DM_parameters = std(total_average_HRF_parameters_thr_T2DM(:,1:10), 0, 2);

avg_psc_sub_CNT_parameters = mean(total_average_HRF_parameters_sub_CNT(:,1:10), 2);
avg_psc_sub_T2DM_parameters = mean(total_average_HRF_parameters_sub_T2DM(:,1:10), 2);
std_psc_sub_CNT_parameters = std(total_average_HRF_parameters_sub_CNT(:,1:10), 0, 2);
std_psc_sub_T2DM_parameters = std(total_average_HRF_parameters_sub_T2DM(:,1:10), 0, 2);

median_psc_thr_CNT_parameters = median(total_median_HRF_parameters_thr_CNT(:,1:10), 2);
median_psc_thr_T2DM_parameters = median(total_median_HRF_parameters_thr_T2DM(:,1:10), 2);
iqr_psc_thr_CNT_parameters = iqr(total_median_HRF_parameters_thr_CNT(:,1:10), 2);
iqr_psc_thr_T2DM_parameters = iqr(total_median_HRF_parameters_thr_T2DM(:,1:10), 2);

median_psc_sub_CNT_parameters = median(total_median_HRF_parameters_sub_CNT(:,1:10), 2);
median_psc_sub_T2DM_parameters = median(total_median_HRF_parameters_sub_T2DM(:,1:10), 2);
iqr_psc_sub_CNT_parameters = iqr(total_median_HRF_parameters_sub_CNT(:,1:10), 2);
iqr_psc_sub_T2DM_parameters = iqr(total_median_HRF_parameters_sub_T2DM(:,1:10), 2);


% Overall average and standard deviation / median and interquartile range
% of the parameters of the average / median HRF of the negative signal 
% change ROIs
avg_nsc_thr_CNT_parameters = mean(total_average_HRF_parameters_thr_CNT(:,11:end), 2);
avg_nsc_thr_T2DM_parameters = mean(total_average_HRF_parameters_thr_T2DM(:,11:end), 2);
std_nsc_thr_CNT_parameters = std(total_average_HRF_parameters_thr_CNT(:,11:end), 0, 2);
std_nsc_thr_T2DM_parameters = std(total_average_HRF_parameters_thr_T2DM(:,11:end), 0, 2);

avg_nsc_sub_CNT_parameters = mean(total_average_HRF_parameters_sub_CNT(:,11:end), 2);
avg_nsc_sub_T2DM_parameters = mean(total_average_HRF_parameters_sub_T2DM(:,11:end), 2);
std_nsc_sub_CNT_parameters = std(total_average_HRF_parameters_sub_CNT(:,11:end), 0, 2);
std_nsc_sub_T2DM_parameters = std(total_average_HRF_parameters_sub_T2DM(:,11:end), 0, 2);

median_nsc_thr_CNT_parameters = median(total_median_HRF_parameters_thr_CNT(:,11:end), 2);
median_nsc_thr_T2DM_parameters = median(total_median_HRF_parameters_thr_T2DM(:,11:end), 2);
iqr_nsc_thr_CNT_parameters = iqr(total_median_HRF_parameters_thr_CNT(:,11:end), 2);
iqr_nsc_thr_T2DM_parameters = iqr(total_median_HRF_parameters_thr_T2DM(:,11:end), 2);

median_nsc_sub_CNT_parameters = median(total_median_HRF_parameters_sub_CNT(:,11:end), 2);
median_nsc_sub_T2DM_parameters = median(total_median_HRF_parameters_sub_T2DM(:,11:end), 2);
iqr_nsc_sub_CNT_parameters = iqr(total_median_HRF_parameters_sub_CNT(:,11:end), 2);
iqr_nsc_sub_T2DM_parameters = iqr(total_median_HRF_parameters_sub_T2DM(:,11:end), 2);


% Plots
figure('Name', sprintf('Grand Average HRF - Positive signal change ROIs'))
shadedErrorBar(time, mean(total_avg_thr_CNT(:,1:10),2), (std(total_avg_thr_CNT(:,1:10),0,2))','lineprops',{'-b','LineWidth', 2},'patchSaturation',0.05)
shadedErrorBar(time, mean(total_avg_thr_T2DM(:,1:10),2), (std(total_avg_thr_T2DM(:,1:10),0,2))','lineprops',{'-r','LineWidth', 2},'patchSaturation',0.05)
shadedErrorBar(time, mean(total_avg_sub_CNT(:,1:10),2), (std(total_avg_sub_CNT(:,1:10),0,2))','lineprops',{'--b','LineWidth', 1},'patchSaturation',0.025)
shadedErrorBar(time, mean(total_avg_sub_T2DM(:,1:10),2), (std(total_avg_sub_T2DM(:,1:10),0,2))','lineprops',{'--r','LineWidth', 1},'patchSaturation',0.025)
xlabel('Time (s)')
ylabel('Beta weights') 
title('Grand Average HRF - Positive signal change ROIs');
legend('CNT Thr','T2DM Thr','CNT Sub','T2DM Sub');

figure('Name', sprintf('Grand Average HRF - Negative signal change ROIs'))
shadedErrorBar(time, mean(total_avg_thr_CNT(:,11:end),2), (std(total_avg_thr_CNT(:,11:end),0,2))','lineprops',{'-b','LineWidth', 2},'patchSaturation',0.05)
shadedErrorBar(time, mean(total_avg_thr_T2DM(:,11:end),2), (std(total_avg_thr_T2DM(:,11:end),0,2))','lineprops',{'-r','LineWidth', 2},'patchSaturation',0.05)
shadedErrorBar(time, mean(total_avg_sub_CNT(:,11:end),2), (std(total_avg_sub_CNT(:,11:end),0,2))','lineprops',{'--b','LineWidth', 1},'patchSaturation',0.025)
shadedErrorBar(time, mean(total_avg_sub_T2DM(:,11:end),2), (std(total_avg_sub_T2DM(:,11:end),0,2))','lineprops',{'--r','LineWidth', 1},'patchSaturation',0.025)
xlabel('Time (s)')
ylabel('Beta weights')  
title('Grand Average HRF - Negative signal change ROIs');
legend('CNT Thr','T2DM Thr','CNT Sub','T2DM Sub');

figure('Name', sprintf('Grand Median HRF - Positive signal change ROIs'))
shadedErrorBar(time, median(total_median_thr_CNT(:,1:10),2), (std(total_median_thr_CNT(:,1:10),0,2))','lineprops',{'-b','LineWidth', 2},'patchSaturation',0.05)
shadedErrorBar(time, median(total_median_thr_T2DM(:,1:10),2), (std(total_median_thr_T2DM(:,1:10),0,2))','lineprops',{'-r','LineWidth', 2},'patchSaturation',0.05)
shadedErrorBar(time, median(total_median_sub_CNT(:,1:10),2), (std(total_median_sub_CNT(:,1:10),0,2))','lineprops',{'--b','LineWidth', 1},'patchSaturation',0.025)
shadedErrorBar(time, median(total_median_sub_T2DM(:,1:10),2), (std(total_median_sub_T2DM(:,1:10),0,2))','lineprops',{'--r','LineWidth', 1},'patchSaturation',0.025)
xlabel('Time (s)')
ylabel('Beta weights')  
title('Grand Median HRF - Positive signal change ROIs');
legend('CNT Thr','T2DM Thr','CNT Sub','T2DM Sub');

figure('Name', sprintf('Grand Median HRF - Negative signal change ROIs'))
shadedErrorBar(time, median(total_median_thr_CNT(:,11:end),2), (std(total_median_thr_CNT(:,11:end),0,2))','lineprops',{'-b','LineWidth', 2},'patchSaturation',0.05)
shadedErrorBar(time, median(total_median_thr_T2DM(:,11:end),2), (std(total_median_thr_T2DM(:,11:end),0,2))','lineprops',{'-r','LineWidth', 2},'patchSaturation',0.05)
shadedErrorBar(time, median(total_median_sub_CNT(:,11:end),2), (std(total_median_sub_CNT(:,11:end),0,2))','lineprops',{'--b','LineWidth', 1},'patchSaturation',0.025)
shadedErrorBar(time, median(total_median_sub_T2DM(:,11:end),2), (std(total_median_sub_T2DM(:,11:end),0,2))','lineprops',{'--r','LineWidth', 1},'patchSaturation',0.025)
xlabel('Time (s)')
ylabel('Beta weights')  
title('Grand Median HRF - Negative signal change ROIs');
legend('CNT Thr','T2DM Thr','CNT Sub','T2DM Sub');


%% COEFFICIENT OF VARIATION:

% In this section we estimate the coefficient of variation for the HRF peak
% amplitude and peak latency

thr_both_sides = [];
sub_both_sides = [];

peak_thr_CNT_vector = [];
peak_thr_T2DM_vector = [];
peak_latency_thr_CNT_vector = [];
peak_latency_thr_T2DM_vector = [];

peak_sub_CNT_vector = [];
peak_sub_T2DM_vector = [];
peak_latency_sub_CNT_vector = [];
peak_latency_sub_T2DM_vector = [];

total_peak_thr_CNT = [];
total_peak_thr_T2DM = [];
total_peak_latency_thr_CNT = [];
total_peak_latency_thr_T2DM = [];

total_peak_sub_CNT = [];
total_peak_sub_T2DM = [];
total_peak_latency_sub_CNT = [];
total_peak_latency_sub_T2DM = [];
    
CV_peak_thr_CNT = [];
CV_peak_thr_T2DM = [];
CV_peak_latency_thr_CNT = [];
CV_peak_latency_thr_T2DM = [];

CV_peak_sub_CNT = [];
CV_peak_sub_T2DM = [];
CV_peak_latency_sub_CNT = [];
CV_peak_latency_sub_T2DM = [];


for r=2:(n_rois+1)
    for i=1:n_subjects
        for j=1:2
            sided_thr = cell2mat(ROIbetas(1+(n_points*2)*(j-1)+(n_points*4)*(i-1):n_points+(n_points*2)*(j-1)+(n_points*4)*(i-1),r));    % each column is a side in each subject
            sided_sub = cell2mat(ROIbetas((n_points+1)+(n_points*2)*(j-1)+(n_points*4)*(i-1):(n_points*2)+(n_points*2)*(j-1)+(n_points*4)*(i-1),r));
            thr_both_sides = [thr_both_sides sided_thr]; 
            sub_both_sides = [sub_both_sides sided_sub]; 
        end

        mean_thr = mean(thr_both_sides,2);                                                                  % mean over rows (2)
        mean_sub = mean(sub_both_sides,2);                                               

        if i<65                                                                                             % until subject 65 --> T2DM
            [peak_thr_T2DM, peak_latency_thr_T2DM] = max(mean_thr(2:5,1));                                  % gets the peak amplitude and peak latency within the 2.5 - 10 second range (where the HRF peak is commonly noticed) in T2DM individuals
            [peak_sub_T2DM, peak_latency_sub_T2DM] = max(mean_sub(2:5,1));                           
            peak_thr_T2DM_vector = [peak_thr_T2DM_vector; peak_thr_T2DM];                                   % stores the HRF peak amplitude of T2DM individuals in the Thr condition
            peak_latency_thr_T2DM_vector = [peak_latency_thr_T2DM_vector; peak_latency_thr_T2DM*2.5];       % stores the HRF peak latency of T2DM individuals in the Thr condition (We multiply by 2.5 to get the time. We don't need to subtract -1 due to the aforementioned reasons)           
            peak_sub_T2DM_vector = [peak_sub_T2DM_vector; peak_sub_T2DM];
            peak_latency_sub_T2DM_vector = [peak_latency_sub_T2DM_vector; peak_latency_sub_T2DM*2.5];
        else                                                                                                % after subject 65 --> CNT
            [peak_thr_CNT, peak_latency_thr_CNT] = max(mean_thr(2:5,1));                                    % gets the peak amplitude and peak latency within the 2.5 - 10 second range (where the HRF peak is commonly noticed) in CNT individuals
            [peak_sub_CNT, peak_latency_sub_CNT] = max(mean_sub(2:5,1));
            peak_thr_CNT_vector = [peak_thr_CNT_vector; peak_thr_CNT];                                      % stores the HRF peak amplitude of CNT individuals in the Thr condition
            peak_latency_thr_CNT_vector = [peak_latency_thr_CNT_vector; peak_latency_thr_CNT*2.5];          % stores the HRF peak latency of CNT individuals in the Thr condition (We multiply by 2.5 to get the time. We don't need to subtract -1 due to the aforementioned reasons)
            peak_sub_CNT_vector = [peak_sub_CNT_vector; peak_sub_CNT];
            peak_latency_sub_CNT_vector = [peak_latency_sub_CNT_vector; peak_latency_sub_CNT*2.5];
        end

        % Delete to not duplicate or interfere with data
        thr_both_sides = [];
        sub_both_sides = [];
    end  

    total_peak_thr_CNT = [total_peak_thr_CNT peak_thr_CNT_vector];                                          % stores every HRF peak amplitude value of CNT subjects in the Thr condition per ROI 
    total_peak_thr_T2DM = [total_peak_thr_T2DM peak_thr_T2DM_vector];                                       % stores every HRF peak amplitude value of T2DM subjects in the Thr condition per ROI 
    total_peak_latency_thr_CNT = [total_peak_latency_thr_CNT peak_latency_thr_CNT_vector];                  % stores every HRF peak latency value of CNT subjects in the Thr condition per ROI 
    total_peak_latency_thr_T2DM = [total_peak_latency_thr_T2DM peak_latency_thr_T2DM_vector];               % stores every HRF peak latency value of T2DM subjects in the Thr condition per ROI 
    
    total_peak_sub_CNT = [total_peak_sub_CNT peak_sub_CNT_vector];                                          % stores every HRF peak amplitude value of CNT subjects in the Sub condition per ROI 
    total_peak_sub_T2DM = [total_peak_sub_T2DM peak_sub_T2DM_vector];                                       % stores every HRF peak amplitude value of T2DM subjects in the Sub condition per ROI 
    total_peak_latency_sub_CNT = [total_peak_latency_sub_CNT peak_latency_sub_CNT_vector];                  % stores every HRF peak latency value of CNT subjects in the Sub condition per ROI 
    total_peak_latency_sub_T2DM = [total_peak_latency_sub_T2DM peak_latency_sub_T2DM_vector];               % stores every HRF peak latency value of T2DM subjects in the Sub condition per ROI 
    
    
    % Delete to not duplicate or interfere with data
    peak_thr_CNT_vector = [];
    peak_thr_T2DM_vector = [];
    peak_latency_thr_CNT_vector = [];    
    peak_latency_thr_T2DM_vector = [];
    
    peak_sub_CNT_vector = [];
    peak_sub_T2DM_vector = [];
    peak_latency_sub_CNT_vector = [];    
    peak_latency_sub_T2DM_vector = [];
end


% Estimating the average and standard deviation of the HRF peak amplitude 
% and peak latency in each group and condition per ROI
avg_peak_thr_CNT = mean(total_peak_thr_CNT,1);
avg_peak_thr_T2DM = mean(total_peak_thr_T2DM,1);
std_peak_thr_CNT = std(total_peak_thr_CNT,0,1);
std_peak_thr_T2DM = std(total_peak_thr_T2DM,0,1);

avg_peak_latency_thr_CNT = mean(total_peak_latency_thr_CNT,1);
avg_peak_latency_thr_T2DM = mean(total_peak_latency_thr_T2DM,1);
std_peak_latency_thr_CNT = std(total_peak_latency_thr_CNT,0,1);
std_peak_latency_thr_T2DM = std(total_peak_latency_thr_T2DM,0,1);

avg_peak_sub_CNT = mean(total_peak_sub_CNT,1);
avg_peak_sub_T2DM = mean(total_peak_sub_T2DM,1);
std_peak_sub_CNT = std(total_peak_sub_CNT,0,1);
std_peak_sub_T2DM = std(total_peak_sub_T2DM,0,1);

avg_peak_latency_sub_CNT = mean(total_peak_latency_sub_CNT,1);
avg_peak_latency_sub_T2DM = mean(total_peak_latency_sub_T2DM,1);
std_peak_latency_sub_CNT = std(total_peak_latency_sub_CNT,0,1);
std_peak_latency_sub_T2DM = std(total_peak_latency_sub_T2DM,0,1);


for r=1:22
    
    % Estimating the Coefficient of Variation for the HRF peak amplitude 
    % and peak latency in each group and condition per ROI

    
    % ***************************** CNT *********************************
    
    % --------------------- THRESHOLD CONDITION -------------------------
    
    avg_peak = avg_peak_thr_CNT(:,r);
    avg_peak_latency = avg_peak_latency_thr_CNT(:,r);
    std_peak = std_peak_thr_CNT(:,r);
    std_peak_latency = std_peak_latency_thr_CNT(:,r);
    
    CV_peak = std_peak/avg_peak;
    CV_latency = std_peak_latency/avg_peak_latency;

    CV_peak_thr_CNT = [CV_peak_thr_CNT CV_peak];
    CV_peak_latency_thr_CNT = [CV_peak_latency_thr_CNT CV_latency];
    
    
    % -------------------- SUBMAXIMUM CONDITION ------------------------
    
    avg_peak = avg_peak_sub_CNT(:,r);
    avg_peak_latency = avg_peak_latency_sub_CNT(:,r);
    std_peak = std_peak_sub_CNT(:,r);
    std_peak_latency = std_peak_latency_sub_CNT(:,r);
    
    CV_peak = std_peak/avg_peak;
    CV_latency = std_peak_latency/avg_peak_latency;
    
    CV_peak_sub_CNT = [CV_peak_sub_CNT CV_peak];
    CV_peak_latency_sub_CNT = [CV_peak_latency_sub_CNT CV_latency];
    
    
    % **************************** T2DM *********************************
    
    % --------------------- THRESHOLD CONDITION -------------------------
    
    avg_peak = avg_peak_thr_T2DM(:,r);
    avg_peak_latency = avg_peak_latency_thr_T2DM(:,r);
    std_peak = std_peak_thr_T2DM(:,r);
    std_peak_latency = std_peak_latency_thr_T2DM(:,r);
    
    CV_peak = std_peak/avg_peak;
    CV_latency = std_peak_latency/avg_peak_latency;
    
    CV_peak_thr_T2DM = [CV_peak_thr_T2DM CV_peak];
    CV_peak_latency_thr_T2DM = [CV_peak_latency_thr_T2DM CV_latency];
    
    
    % -------------------- SUBMAXIMUM CONDITION ------------------------
    
    avg_peak = avg_peak_sub_T2DM(:,r);
    avg_peak_latency = avg_peak_latency_sub_T2DM(:,r);
    std_peak = std_peak_sub_T2DM(:,r);
    std_peak_latency = std_peak_latency_sub_T2DM(:,r);
    
    CV_peak = std_peak/avg_peak;
    CV_latency = std_peak_latency/avg_peak_latency;
    
    CV_peak_sub_T2DM = [CV_peak_sub_T2DM CV_peak];
    CV_peak_latency_sub_T2DM = [CV_peak_latency_sub_T2DM CV_latency];
end


% Plotting the Coefficient of Variation graphs with data concerning each
% condition and group per ROI

% ««««««««««««««««««««««««« Peak Amplitude »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

figure
scatter(1:22,CV_peak_thr_CNT,30,'MarkerEdgeColor',[0 .5 .5],...
                                'MarkerFaceColor',[0 .7 .7],...
                                'LineWidth',1)
hold on
scatter(1:22,CV_peak_thr_T2DM,30,'MarkerEdgeColor',[0.64 0.08 0.18],...
                                 'MarkerFaceColor',[0.87 0.19 0.19],...
                                 'LineWidth',1)
hold on
scatter(1:22,CV_peak_sub_CNT,30,'MarkerEdgeColor',[0 .5 .5],...
                                'LineWidth',0.75)
hold on
scatter(1:22,CV_peak_sub_T2DM,30,'MarkerEdgeColor',[0.64 0.08 0.18],...
                                 'LineWidth',0.75)
                            
set(gca, 'XTick', 1:22, 'XTickLabel', ROI_regions);     % labels the x axis as the ROI names
xtickangle(45)                                          % rotates the labels by 45º
title('Coefficient of Variation - Peak Amplitude');
xlabel('ROI')
ylabel('Coefficient of Variation')  
legend ('CNT Thr','T2DM Thr', 'CNT Sub', 'T2DM Sub')


% «««««««««««««««««««««««««« Peak Latency »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

figure
scatter(1:22,CV_peak_latency_thr_CNT,30,'MarkerEdgeColor',[0 .5 .5],...
                                        'MarkerFaceColor',[0 .7 .7],...
                                        'LineWidth',1)
hold on
scatter(1:22,CV_peak_latency_thr_T2DM,30,'MarkerEdgeColor',[0.64 0.08 0.18],...
                                         'MarkerFaceColor',[0.87 0.19 0.19],...
                                         'LineWidth',1)
hold on                            
scatter(1:22,CV_peak_latency_sub_CNT,30,'MarkerEdgeColor',[0 .5 .5],...
                                        'LineWidth',0.75)
hold on
scatter(1:22,CV_peak_latency_sub_T2DM,30,'MarkerEdgeColor',[0.64 0.08 0.18],...
                                         'LineWidth',0.75)
                            
set(gca, 'XTick', 1:22, 'XTickLabel', ROI_regions);
xtickangle(45)
title('Coefficient of Variation - Peak Latency');
xlabel('ROI')
ylabel('Coefficient of Variation')  
legend ('CNT Thr','T2DM Thr', 'CNT Sub', 'T2DM Sub')