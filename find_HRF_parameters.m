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


%% CONSTANTS:

n_subjects = 141;
n_T2DM = 64;
n_CNT = 77;
n_rois = 22;
n_points = 8;
time = (0:n_points-1)* 2.5;                              % volumes per subject * TR (2.5 s)

ROI_regions = {'L IPL BA40', 'L Insula BA13', 'L Precuneus BA7', 'R IFG BA9', 'R MFG BA8', 'R MFG BA46', 'R MT BA19', ... 
               'R SFG BA6', 'R SPL BA7', 'R V2 BA18', 'L AC BA32', 'L CG BA31', 'L PC BA30', 'L PrcL BA5', 'L PrhG BA36', ...
               'R CG BA24', 'R Insula BA13', 'R PC BA23', 'R PrecG BA4', 'R PrecG BA43', 'R PrimSens BA1', 'R STG BA39'};

          
%% INITIALIZATION:

% In this section, we get every mean value of beta per subject in each
% group and condition for each ROI

thr_both_sides = zeros(n_points,2);                         % 2 since each column represents a stimulation side - left and right
sub_both_sides = zeros(n_points,2);

mean_per_subject_thr_CNT = zeros(n_points,n_CNT*n_rois);
mean_per_subject_sub_CNT = zeros(n_points,n_CNT*n_rois);
mean_per_subject_thr_T2DM = zeros(n_points,n_T2DM*n_rois);
mean_per_subject_sub_T2DM = zeros(n_points,n_T2DM*n_rois);


for r=2:n_rois+1
    for i=1:n_subjects
        for j=1:2
            sided_thr = cell2mat(ROIbetas(1+(n_points*2)*(j-1)+(n_points*4)*(i-1):n_points+(n_points*2)*(j-1)+(n_points*4)*(i-1),r));    % each column is a side in each subject
            sided_sub = cell2mat(ROIbetas((n_points+1)+(n_points*2)*(j-1)+(n_points*4)*(i-1):(n_points*2)+(n_points*2)*(j-1)+(n_points*4)*(i-1),r));
            thr_both_sides(:,j) = sided_thr; 
            sub_both_sides(:,j) = sided_sub; 
        end

        mean_thr = mean(thr_both_sides,2);                                                      % mean over rows (2)
        mean_sub = mean(sub_both_sides,2);                                               

        if i<65                                                                                 % until subject 65 --> T2DM 
            mean_per_subject_thr_T2DM(:,i+n_T2DM*(r-2)) = mean_thr;                             % every mean value of beta per T2DM subject in the Thr condition for each ROI 
            mean_per_subject_sub_T2DM(:,i+n_T2DM*(r-2)) = mean_sub;
        else                                                                                    % after subject 65 --> CNT
            mean_per_subject_thr_CNT(:,i-n_T2DM+n_CNT*(r-2)) = mean_thr;                        % every mean value of beta per CNT subject in the Thr condition for each ROI 
            mean_per_subject_sub_CNT(:,i-n_T2DM+n_CNT*(r-2)) = mean_sub;
        end
    
        % Delete to not duplicate or interfere with data
        thr_both_sides = zeros(n_points,2);
        sub_both_sides = zeros(n_points,2);
    end  
end


%% PARAMETER ESTIMATION:

% In this section, we estimate and store the HRF parameters of interest 
% per subject in each group and condition in every ROI

total_parameters_thr_CNT = zeros(n_CNT,11*n_rois);
total_parameters_sub_CNT = zeros(n_CNT,11*n_rois);
total_parameters_thr_T2DM = zeros(n_T2DM,11*n_rois);
total_parameters_sub_T2DM = zeros(n_T2DM,11*n_rois);


% ********************************* CNT **********************************
        
for c = 1:n_CNT
    for r = 1:n_rois                                                                                                                    % 22 ROI
        
        % Stores HRF curves per CNT subjects in each condition and in 
        % every ROI 
        data.CNT.subject(c).ROI(r).Threshold.HRFCurve = mean_per_subject_thr_CNT(1:n_points,r+22*(c-1));
        data.CNT.subject(c).ROI(r).Submaximum.HRFCurve = mean_per_subject_sub_CNT(1:n_points,r+22*(c-1));

        
        % -------------------- THRESHOLD CONDITION ----------------------

        [peak,volume] = max(data.CNT.subject(c).ROI(r).Threshold.HRFCurve(2:6));                                                        % peak: maximum within the sectioned time range (2.5 to 12.5 s - 2nd to 6th index since indexing starts at 1)
        latency = volume * 2.5;                                                                                                         % peak latency: time to peak (volume where the maximum takes place * TR)
        slope_to_peak = minus(peak,data.CNT.subject(c).ROI(r).Threshold.HRFCurve(1))/latency;                                           % relative slope to peak: (peak-1st point)/peak latency
        total_area = trapz((0:n_points-1)*2.5,data.CNT.subject(c).ROI(r).Threshold.HRFCurve(1:n_points)- ... 
                            min(data.CNT.subject(c).ROI(r).Threshold.HRFCurve(1:n_points)));                                            % area under the curve (trapezoidal method)
        pos_area = positive_area((0:n_points-1)*2.5,data.CNT.subject(c).ROI(r).Threshold.HRFCurve(1:n_points));                         % area of the positive sections of the HRF curve
        first_pos_area = positive_area((0:2)*2.5,data.CNT.subject(c).ROI(r).Threshold.HRFCurve(1:3));                                   % area of the first positive section of the HRF curve (from 0 to 5 s)
        second_pos_area = positive_area((2:4)*2.5,data.CNT.subject(c).ROI(r).Threshold.HRFCurve(3:5));                                  % area of the second section of the HRF curve (from 5 to 10 s)
        third_pos_area = positive_area((4:7)*2.5,data.CNT.subject(c).ROI(r).Threshold.HRFCurve(5:n_points));                            % area of the third section of the HRF curve (from 10 to 17.5 s)
        neg_area = negative_area((0:n_points-1)*2.5,data.CNT.subject(c).ROI(r).Threshold.HRFCurve(1:n_points));                         % area of the negative sections of the HRF curve 
        initial_dip_area = negative_area((0:volume)*2.5,data.CNT.subject(c).ROI(r).Threshold.HRFCurve(1:(volume+1)));                   % area of the initial dip of the HRF curve
        undershoot_area = negative_area((volume:7)*2.5,data.CNT.subject(c).ROI(r).Threshold.HRFCurve((volume+1):n_points));             % area of the undershoot of the HRF curve
                
        % Gathers and stores every HRF parameter in the Thr condition per 
        % ROI
        data.CNT.subject(c).ROI(r).Threshold.HRFParameters = [peak,latency,slope_to_peak,total_area,pos_area,first_pos_area,second_pos_area,third_pos_area,neg_area,initial_dip_area,undershoot_area];
        total_parameters_thr_CNT(c,1+11*(r-1):11+11*(r-1)) = data.CNT.subject(c).ROI(r).Threshold.HRFParameters;
       
        
        % -------------------- SUBMAXIMUM CONDITION ---------------------

        [peak,volume] = max(data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(2:6));                             
        latency = volume * 2.5;                                                                             
        slope_to_peak = minus(peak,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(1))/latency;                                  
        total_area = trapz((0:n_points-1)*2.5,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(1:n_points)- ... 
                            min(data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(1:n_points)));    
        pos_area = positive_area((0:n_points-1)*2.5,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(1:n_points));      
        first_pos_area = positive_area((0:2)*2.5,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(1:3));                   
        second_pos_area = positive_area((2:4)*2.5,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(3:5));                  
        third_pos_area = positive_area((4:7)*2.5,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(5:n_points));            
        neg_area = negative_area((0:n_points-1)*2.5,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(1:n_points));      
        initial_dip_area = negative_area((0:volume)*2.5,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve(1:(volume+1)));                   
        undershoot_area = negative_area((volume:7)*2.5,data.CNT.subject(c).ROI(r).Submaximum.HRFCurve((volume+1):n_points));           
        
        % Gathers and stores every HRF parameter in the Sub condition per 
        % ROI
        data.CNT.subject(c).ROI(r).Submaximum.HRFParameters = [peak,latency,slope_to_peak,total_area,pos_area,first_pos_area,second_pos_area,third_pos_area,neg_area,initial_dip_area,undershoot_area];       
        total_parameters_sub_CNT(c,1+11*(r-1):11+11*(r-1)) = data.CNT.subject(c).ROI(r).Submaximum.HRFParameters;
    end
    
    % Average HRF parameters per condition
    avg_parameters_thr_CNT = mean(total_parameters_thr_CNT, 1);
    avg_parameters_sub_CNT = mean(total_parameters_sub_CNT, 1);   
end


% ******************************** T2DM **********************************


for t=1:n_T2DM
    for r = 1:n_rois
        
        % Storing HRF curves per T2DM subjects in each condition and in 
        % each ROI
        data.T2DM.subject(t).ROI(r).Threshold.HRFCurve = mean_per_subject_thr_T2DM(1:n_points,r+22*(t-1));
        data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve = mean_per_subject_sub_T2DM(1:n_points,r+22*(t-1));
        
        % -------------------- THRESHOLD CONDITION ----------------------

        [peak,volume] = max(data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(2:6));                              
        latency = volume * 2.5;                                                                             
        slope_to_peak = minus(peak,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(1))/latency;
        total_area = trapz((0:n_points-1)*2.5,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(1:n_points)- ... 
                            min(data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(1:n_points)));  
        pos_area = positive_area((0:n_points-1)*2.5,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(1:n_points));
        first_pos_area = positive_area((0:2)*2.5,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(1:3));                   
        second_pos_area = positive_area((2:4)*2.5,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(3:5));                  
        third_pos_area = positive_area((4:7)*2.5,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(5:n_points));            
        neg_area = negative_area((0:n_points-1)*2.5,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(1:n_points));         
        initial_dip_area = negative_area((0:volume)*2.5,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve(1:(volume+1)));                  
        undershoot_area = negative_area((volume:7)*2.5,data.T2DM.subject(t).ROI(r).Threshold.HRFCurve((volume+1):n_points));           
        
        % Gathers and stores every HRF parameter in the Thr condition per 
        % ROI
        data.T2DM.subject(t).ROI(r).Threshold.HRFParameters = [peak,latency,slope_to_peak,total_area,pos_area,first_pos_area,second_pos_area,third_pos_area,neg_area,initial_dip_area,undershoot_area];        
        total_parameters_thr_T2DM(t,1+11*(r-1):11+11*(r-1)) = data.T2DM.subject(t).ROI(r).Threshold.HRFParameters;

        
        % -------------------- SUBMAXIMUM CONDITION ---------------------

        [peak,volume] = max(data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(2:6));                             
        latency = volume * 2.5;                                                                             
        slope_to_peak = minus(peak,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(1))/latency;                                  
        total_area = trapz((0:n_points-1)*2.5,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(1:n_points)- ... 
                            min(data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(1:n_points)));                         
        pos_area = positive_area((0:n_points-1)*2.5,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(1:n_points));      
        first_pos_area = positive_area((0:2)*2.5,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(1:3));                   
        second_pos_area = positive_area((2:4)*2.5,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(3:5));                  
        third_pos_area = positive_area((4:7)*2.5,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(5:n_points)); 
        neg_area = negative_area((0:n_points-1)*2.5,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(1:n_points));      
        initial_dip_area = negative_area((0:volume)*2.5,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve(1:(volume+1)));                   
        undershoot_area = negative_area((volume:7)*2.5,data.T2DM.subject(t).ROI(r).Submaximum.HRFCurve((volume+1):n_points));       
        
        % Gathers and stores every HRF parameter in the Sub condition per 
        % ROI        
        data.T2DM.subject(t).ROI(r).Submaximum.HRFParameters = [peak,latency,slope_to_peak,total_area,pos_area,first_pos_area,second_pos_area,third_pos_area,neg_area,initial_dip_area,undershoot_area];
        total_parameters_sub_T2DM(t,1+11*(r-1):11+11*(r-1)) = data.CNT.subject(t).ROI(r).Submaximum.HRFParameters;

    end
    
    % Average HRF parameters per condition
    avg_parameters_thr_T2DM = mean(total_parameters_thr_T2DM, 1);
    avg_parameters_sub_T2DM = mean(total_parameters_sub_T2DM, 1);
end


%% FINAL STORAGING:

% In this section, we create tables regarding the HRF parameters of 
% interest and its information per group and condition in each subject

% Storaging order:        ROI1 p.1 || ROI1 p.2 || ... || ROI22 p.11
%                                       1. CNT Thr
%                                       2. CNT Sub
%                                       3. T2DM Thr
%                                       4. T2DM Sub


% ************************** Individual data ****************************

% «««««««««««««««««««««««««««« Table data »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

group = cell(n_subjects*2,1);                                                               

CNT_total_parameters = zeros(n_CNT*2, 11*n_rois);                                               % gets the HRF parameters of all CNT subjects in both conditions and in every ROI
T2DM_total_parameters = zeros(n_T2DM*2, 11*n_rois);


for c = 1:n_CNT
    
    % Gathers every parameter of interest of every CNT subject in each condition and in every ROI
    for r = 1:n_rois
        CNT_total_parameters((c*2)-1:c*2,1+11*(r-1):11+11*(r-1)) = [data.CNT.subject(c).ROI(r).Threshold.HRFParameters; data.CNT.subject(c).ROI(r).Submaximum.HRFParameters];
    end
    
    % Identifies each subject's group for each condition
    for l = 1:2
        group{l+2*(c-1)}='CNT';
    end
end


for t = 1:n_T2DM
    
    % Gathers every parameter of interest of every T2DM subject in each condition and in every ROI
    for r = 1:n_rois
        T2DM_total_parameters((t*2)-1:t*2,1+11*(r-1):11+11*(r-1)) = [data.T2DM.subject(t).ROI(r).Threshold.HRFParameters; data.T2DM.subject(t).ROI(r).Submaximum.HRFParameters];
    end
    
    % Identifies each subject's group for each condition
    for l = 1:2
        group{(n_CNT*2)+(l+2*(t-1))}='T2DM';
    end
end

% Gathers every parameter of interest of all subjects in each condition and in every ROI
HRF_data = [CNT_total_parameters; T2DM_total_parameters];


% ««««««««««««««««««««««««««« Table setting »»»»»»»»»»»»»»»»»»»»»»»»»»»»»

headers = cell(1,242);                                                                                                                                                                     % cell array for the table headings
parameters = {' peak amplitude', ' peak latency', ' relative slope to peak', ' AUC', ' APCS', ' fPCSA', ' sPCSA', ' tPCSA', ' ANCS', ' initial dip area', ' undershoot area'};              % parameters of interest

for r = 1:n_rois

    % Automatically sets the table headers
    for par = 1:length(parameters)
        headers{par+11*(r-1)} = char(strrep(strcat(ROI_regions(r),parameters(par)),' ','_'));       % concatenates strings, replaces spaces by _ and converts to char
    end
end

% Forms tables
covariates_HRFdata = array2table(HRF_data,'VariableNames', headers);                                % table with the HRF data from every subject in each condition and in every ROI as well as its headers 
covariates_group = array2table(group,'VariableNames', {'Group'});                                   % table with groups (CNT and T2DM)

% Merges tables
covariates_HRFdata = [covariates_group covariates_HRFdata];


% *************************** Average data *****************************

% «««««««««««««««««««««««««««« Table data »»»»»»»»»»»»»»»»»»»»»»»»»»»»»»»

avg_HRF_data = [avg_parameters_thr_CNT; avg_parameters_sub_CNT; avg_parameters_thr_T2DM; avg_parameters_sub_T2DM];

% ««««««««««««««««««««««««««« Table setting »»»»»»»»»»»»»»»»»»»»»»»»»»»»»

avg_headers = headers;
avg_group = {'CNT','CNT','T2DM','T2DM'}';

% Forms tables
avg_covariates_HRFdata = array2table(avg_HRF_data,'VariableNames', avg_headers);                   % table with the average HRF data in each condition and in every ROI as well as its headers       
avg_covariates_group = array2table(avg_group,'VariableNames', {'Group'});                

% Merges tables
avg_covariates_HRFdata = [avg_covariates_group avg_covariates_HRFdata];