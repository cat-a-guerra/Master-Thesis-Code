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

load('ROI_regions.mat');
load('HRF_parameters.mat');
load('headers.mat');
load('covariates_HRFdata.mat');


%% CONSTANTS:

n_subjects = 141;
n_T2DM = 64;
n_CNT = 77;
n_rois = 22;


%% STATISTICS:

% In this section, we estimate the statistics of the HRF parameters of
% interest from every subject in each group and condition and in every ROI. 
% First, we apply a Shapiro-Wilk test in order to assess if the data has a 
% normal distribution per group and condition. Then, we assess the 
% significant difference for the same HRF parameters between groups. To do 
% so, we apply a Two-sample T-test or a Wilcoxon ranksum test, depending on 
% whether the data regarding the same HRF parameter in both groups of the 
% same condition has a normal distribution or not, respectively.


% HRF parameters data from each group and condition in every ROI
CNT_thr = zeros(n_CNT,n_rois*length(HRF_parameters));
T2DM_thr = zeros(n_T2DM,n_rois*length(HRF_parameters));
CNT_sub = zeros(n_CNT,n_rois*length(HRF_parameters));
T2DM_sub = zeros(n_T2DM,n_rois*length(HRF_parameters));

% Shapiro-Wilk test results of the HRF parameters from each group and 
% condition in every ROI 
sw_CNT_thr = zeros(2,n_rois*length(HRF_parameters));
sw_T2DM_thr = zeros(2,n_rois*length(HRF_parameters));
sw_CNT_sub = zeros(2,n_rois*length(HRF_parameters));
sw_T2DM_sub = zeros(2,n_rois*length(HRF_parameters));


% Column numbers for the HRF parameters with normal distribution from each 
% group and condition in every ROI
avg_std_indexes_CNT_thr = [];
avg_std_indexes_T2DM_thr = [];
avg_std_indexes_CNT_sub = [];
avg_std_indexes_T2DM_sub = [];

% Column numbers for the HRF parameters with non-normal distribution from 
% each group and condition in every ROI
median_iqr_indexes_CNT_thr = [];
median_iqr_indexes_T2DM_thr = [];
median_iqr_indexes_CNT_sub = [];
median_iqr_indexes_T2DM_sub = [];


% Average + standard-deviation of the normal HRF parameters // Median + 
% interquartile range of the non-normal HRF parameters from each group and
% condition in every ROI
avg_parameters_CNT_thr = [];
avg_parameters_T2DM_thr = [];
std_parameters_CNT_thr = [];
std_parameters_T2DM_thr = [];

median_parameters_CNT_thr = [];
median_parameters_T2DM_thr = [];
iqr_parameters_CNT_thr = [];
iqr_parameters_T2DM_thr = [];

avg_parameters_CNT_sub = [];
avg_parameters_T2DM_sub = [];
std_parameters_CNT_sub = [];
std_parameters_T2DM_sub = [];

median_parameters_CNT_sub = [];
median_parameters_T2DM_sub = [];
iqr_parameters_CNT_sub = [];
iqr_parameters_T2DM_sub = [];


% Column numbers for the HRF parameters with normal // non-normal 
% distribution from different groups but same condition in every ROI
ttest_indexes_thr = [];
ttest_indexes_sub = [];

wilcoxon_indexes_thr = [];
wilcoxon_indexes_sub = [];


% T-test and Wilcoxon test results of the HRF parameters from each 
% condition in every ROI
ttest_thr = [];
ttest_sub = [];

wilcoxon_thr = [];
wilcoxon_sub = [];



for w=2:width(covariates_HRFdata)
    
    
    % --------------------- THRESHOLD CONDITION -------------------------
        
    CNT_thr(1:n_CNT,w-1) = table2array(covariates_HRFdata(1:2:154,w));                                    
    
    [h,p] = swtest(table2array(covariates_HRFdata(1:2:154,w)), 0.05, true);                                             % estimates h and p of the Shapiro-Wilk test for each HRF parameter of the CNT subjects in every ROI
    sw_CNT_thr(1:2,w-1) = [h;p];
    
    if h==0                                                                                                             % in case of a normal distribution --> group average of the HRF parameters
        avg_std_indexes_CNT_thr = [avg_std_indexes_CNT_thr (w-1)];
        avg_parameters_CNT_thr = [avg_parameters_CNT_thr mean(table2array(covariates_HRFdata(1:2:154,w)),1)];
        std_parameters_CNT_thr = [std_parameters_CNT_thr std(table2array(covariates_HRFdata(1:2:154,w)),0,1)];
    else                                                                                                                % in case of a non-normal distribution --> group median of the HRF parameters
        median_iqr_indexes_CNT_thr = [median_iqr_indexes_CNT_thr (w-1)];
        median_parameters_CNT_thr = [median_parameters_CNT_thr median(table2array(covariates_HRFdata(1:2:154,w)),1)];
        iqr_parameters_CNT_thr = [iqr_parameters_CNT_thr iqr(table2array(covariates_HRFdata(1:2:154,w)),1)];
    end
            

    T2DM_thr(1:n_T2DM,w-1) = table2array(covariates_HRFdata(155:2:end,w));
    
    [h,p] = swtest(table2array(covariates_HRFdata(155:2:end,w)), 0.05, true);                                           % estimates h and p of the Shapiro-Wilk test for each HRF parameter of the T2DM subjects in every ROI                               
    sw_T2DM_thr(1:2,w-1) = [h;p];
    
    if h==0
        avg_std_indexes_T2DM_thr = [avg_std_indexes_T2DM_thr (w-1)];
        avg_parameters_T2DM_thr = [avg_parameters_T2DM_thr mean(table2array(covariates_HRFdata(155:2:end,w)),1)];
        std_parameters_T2DM_thr = [std_parameters_T2DM_thr std(table2array(covariates_HRFdata(155:2:end,w)),0,1)];
    else
        median_iqr_indexes_T2DM_thr = [median_iqr_indexes_T2DM_thr (w-1)];
        median_parameters_T2DM_thr = [median_parameters_T2DM_thr median(table2array(covariates_HRFdata(155:2:end,w)),1)];
        iqr_parameters_T2DM_thr = [iqr_parameters_T2DM_thr iqr(table2array(covariates_HRFdata(155:2:end,w)),1)];
    end
    
    
    % Testing significant differences for the same HRF parameters between 
    % groups
    if sw_CNT_thr(1,w-1)==0 && sw_T2DM_thr(1,w-1)==0                                                                    % in case of a normal distribution in the same HRF parameter from different groups --> Two sample T-test
        ttest_indexes_thr = [ttest_indexes_thr (w-1)];
        [h,p] = ttest2(CNT_thr(:,w-1),T2DM_thr(:,w-1));
        ttest_thr = [ttest_thr [h;p]];
    else                                                                                                                % in case of a non-normal distribution in the same HRF parameter from different groups --> Wilcoxon ranksum test
        wilcoxon_indexes_thr = [wilcoxon_indexes_thr (w-1)];
        [p,h] = ranksum(CNT_thr(:,w-1),T2DM_thr(:,w-1));
        wilcoxon_thr = [wilcoxon_thr [h;p]];
    end
        
    
    % --------------------- SUBMAXIMUM CONDITION ------------------------
    
    CNT_sub(1:n_CNT,w-1) = table2array(covariates_HRFdata(2:2:154,w));                                    
    
    [h,p] = swtest(table2array(covariates_HRFdata(2:2:154,w)), 0.05, true);                                             
    sw_CNT_sub(1:2,w-1) = [h;p];
    
    if h==0                                                                                                             
        avg_std_indexes_CNT_sub = [avg_std_indexes_CNT_sub (w-1)];
        avg_parameters_CNT_sub = [avg_parameters_CNT_sub mean(table2array(covariates_HRFdata(2:2:154,w)),1)];
        std_parameters_CNT_sub = [std_parameters_CNT_sub std(table2array(covariates_HRFdata(2:2:154,w)),0,1)];
    else                                                                                                                
        median_iqr_indexes_CNT_sub = [median_iqr_indexes_CNT_sub (w-1)];
        median_parameters_CNT_sub = [median_parameters_CNT_sub median(table2array(covariates_HRFdata(2:2:154,w)),1)];
        iqr_parameters_CNT_sub = [iqr_parameters_CNT_sub iqr(table2array(covariates_HRFdata(2:2:154,w)),1)];
    end  
    
    
    T2DM_sub(1:n_T2DM,w-1) = table2array(covariates_HRFdata(156:2:end,w));

    [h,p] = swtest(table2array(covariates_HRFdata(156:2:end,w)), 0.05, true);                                           
    sw_T2DM_sub(1:2,w-1) = [h;p];
    
    if h==0
        avg_std_indexes_T2DM_sub = [avg_std_indexes_T2DM_sub (w-1)];
        avg_parameters_T2DM_sub = [avg_parameters_T2DM_sub mean(table2array(covariates_HRFdata(156:2:end,w)),1)];
        std_parameters_T2DM_sub = [std_parameters_T2DM_sub std(table2array(covariates_HRFdata(156:2:end,w)),0,1)];
    else
        median_iqr_indexes_T2DM_sub = [median_iqr_indexes_T2DM_sub (w-1)];
        median_parameters_T2DM_sub = [median_parameters_T2DM_sub median(table2array(covariates_HRFdata(156:2:end,w)),1)];
        iqr_parameters_T2DM_sub = [iqr_parameters_T2DM_sub iqr(table2array(covariates_HRFdata(156:2:end,w)),1)];
    end
    
    
    % Testing significant differences for the same HRF parameters between 
    % groups
    if sw_CNT_sub(1,w-1)==0 && sw_T2DM_sub(1,w-1)==0                        
        ttest_indexes_sub = [ttest_indexes_sub (w-1)];
        [h,p] = ttest2(CNT_sub(:,w-1),T2DM_sub(:,w-1));
        ttest_sub = [ttest_sub [h;p]];
    else                                                                    
        wilcoxon_indexes_sub = [wilcoxon_indexes_sub (w-1)];
        [p,h] = ranksum(CNT_sub(:,w-1),T2DM_sub(:,w-1));
        wilcoxon_sub = [wilcoxon_sub [h;p]];
    end
end



%% FINAL STORAGING:

% In this section, we create tables regarding the HRF parameters of
% interest, as well as the Shapiro-Wilk test, Two-sample T-test and 
% Wilcoxon ranksum test results for them.


% Column with the statistical test's parameter names to identify its 
% corresponding results    
test_values = array2table({'Hypothesis test result', 'p-value'}','VariableNames', {'Value'});  


% ************* HRF parameters data and Shapiro-Wilk tables **************
    
% ««««««««««««««««««««« HRF parameters data tables »»»»»»»»»»»»»»»»»»»»»»» 
    
% Forms tables with the HRF parameters data from every subject in each 
% group and condition and in every ROI as well as its headers 
CNT_thr_table = array2table(CNT_thr,'VariableNames', headers);              
CNT_sub_table = array2table(CNT_sub,'VariableNames', headers);          
T2DM_thr_table = array2table(T2DM_thr,'VariableNames', headers);         
T2DM_sub_table = array2table(T2DM_sub,'VariableNames', headers);         
    
    
% «««««««««««««««««««««««« Shapiro-Wilk tables »»»»»»»»»»»»»»»»»»»»»»»»»»» 

% Forms tables with the Shapiro-Wilk test results for the HRF parameters 
% data from every subject in each group and condition and in every ROI 
% as well as its headers 
sw_CNT_thr_table = array2table(sw_CNT_thr,'VariableNames', headers);               
sw_CNT_sub_table = array2table(sw_CNT_sub,'VariableNames', headers);       
sw_T2DM_thr_table = array2table(sw_T2DM_thr,'VariableNames', headers);      
sw_T2DM_sub_table = array2table(sw_T2DM_sub,'VariableNames', headers);       
                
% Merges tables
sw_CNT_thr_table = [test_values sw_CNT_thr_table];
sw_CNT_sub_table = [test_values sw_CNT_sub_table];
sw_T2DM_thr_table = [test_values sw_T2DM_thr_table];
sw_T2DM_sub_table = [test_values sw_T2DM_sub_table];
    
    

% ********* Two sample T-test and Wilcoxon ranksum test tables ***********
    
% ««««««««««««««««««««««««« Two sample T-test »»»»»»»»»»»»»»»»»»»»»»»»»»» 
    
if isempty(ttest_thr)==0                                                                    % when a HRF parameter from different groups in the Thr condition has normal distribution  
        
    ttest_thr_headers = cell(1,size(ttest_thr,2));                                          % headers for the table

    for l = 1:length(ttest_indexes_thr)
        ttest_thr_headers(l) = headers(ttest_indexes_thr(l));                               % gets the headers of the normal Thr parameters
    end

    % Forms tables with the Two sample T-test results for the HRF
    % parameters with normal distribution from different groups in the Thr 
    % condition as well as its headers 
    ttest_thr_table = array2table(ttest_thr,'VariableNames', ttest_thr_headers);            

    % Merges the previous table with the statistical test's parameter names 
    ttest_thr_table = [test_values ttest_thr_table];
end
    
    
if isempty(ttest_sub)==0                                                                    % when a HRF parameter from different groups in the Sub condition has normal distribution       

    ttest_sub_headers = cell(1,size(ttest_sub,2));                                          

    for l = 1:length(ttest_indexes_sub)
        ttest_sub_headers(l) = headers(ttest_indexes_sub(l));                               % gets the headers of the normal Sub parameters
    end
    
    % Forms tables with the Two sample T-test results for the HRF
    % parameters with normal distribution from different groups in the Sub 
    % condition as well as its headers
    ttest_sub_table = array2table(ttest_sub,'VariableNames', ttest_sub_headers);            

    % Merges the previous table with the statistical test's parameter names 
    ttest_sub_table = [test_values ttest_sub_table];
end
    
    
% ««««««««««««««««««««««« Wilcoxon ranksum test »»»»»»»»»»»»»»»»»»»»»»»»»» 

if isempty(wilcoxon_thr)==0                                                                     % when a HRF parameter from different groups in the Thr condition has non-normal distribution 

    wilcoxon_thr_headers = cell(1,size(wilcoxon_thr,2));                                       

    for l = 1:length(wilcoxon_indexes_thr)
        wilcoxon_thr_headers(l) = headers(wilcoxon_indexes_thr(l));                             % gets the headers of non-normal Thr HRF parameters
    end

    % Forms tables with the Wilcoxon ranksum test results for the HRF
    % parameters with non-normal distribution from different groups in the 
    % Thr condition as well as its headers
    wilcoxon_thr_table = array2table(wilcoxon_thr,'VariableNames', wilcoxon_thr_headers);       

    % Merges the previous table with the statistical test's parameter names 
    wilcoxon_thr_table = [test_values wilcoxon_thr_table];
end


if isempty(wilcoxon_sub)==0                                                                     % when a HRF parameter from different groups in the Sub condition has non-normal distribution  

    wilcoxon_sub_headers = cell(1,size(wilcoxon_sub,2));                                        

    for l = 1:length(wilcoxon_indexes_sub)
        wilcoxon_sub_headers(l) = headers(wilcoxon_indexes_sub(l));                             % gets the headers of non-normal Sub HRF parameters
    end

    % Forms tables with the Wilcoxon ranksum test results for the HRF
    % parameters with non-normal distribution from different groups in the 
    % Sub condition as well as its headers
    wilcoxon_sub_table = array2table(wilcoxon_sub,'VariableNames', wilcoxon_sub_headers);       

    % Merges the previous table with the statistical test's parameter names 
    wilcoxon_sub_table = [test_values wilcoxon_sub_table];
end
    
    

% ************************** HRF Parameters ******************************
    

% «««««««««««««««««««« Average and standard deviation »»»»»»»»»»»»»»»»»»»»


% -------------------------------- CNT -----------------------------------

if isempty(avg_parameters_CNT_thr)==0                                                                   % when a CNT Thr HRF parameter has normal distribution --> average and standard deviation of the parameter per group and condition in each ROI

    CNT_thr_avg_headers = cell(1,size(avg_parameters_CNT_thr,2));                                       % headers for the table

    for l = 1:length(avg_std_indexes_CNT_thr)
        CNT_thr_avg_headers(l) = headers(avg_std_indexes_CNT_thr(l));                                   % gets the headers of normal CNT Thr HRF parameters
    end

    % Forms tables with the average and standard deviation for the CNT Thr
    % HRF parameters with normal distribution as well as its headers
    CNT_thr_avg_table = array2table(avg_parameters_CNT_thr,'VariableNames', CNT_thr_avg_headers);       
    CNT_thr_std_table = array2table(std_parameters_CNT_thr,'VariableNames', CNT_thr_avg_headers);              
end


if isempty(avg_parameters_CNT_sub)==0                                                                   % when a CNT Sub HRF parameter has normal distribution --> average and standard deviation of the parameter per group and condition in each ROI      

    CNT_sub_avg_headers = cell(1,size(avg_parameters_CNT_sub,2));                                         

    for l = 1:length(avg_std_indexes_CNT_sub)
        CNT_sub_avg_headers(l) = headers(avg_std_indexes_CNT_sub(l));                                   % gets the headers of normal CNT Sub HRF parameters
    end

    % Forms tables with the average and standard deviation for the CNT Sub
    % HRF parameters with normal distribution as well as its headers
    CNT_sub_avg_table = array2table(avg_parameters_CNT_sub,'VariableNames', CNT_sub_avg_headers);       
    CNT_sub_std_table = array2table(std_parameters_CNT_sub,'VariableNames', CNT_sub_avg_headers);         
end
    
    
% -------------------------------- T2DM ----------------------------------

if isempty(avg_parameters_T2DM_thr)==0                                                                  % when a T2DM Thr HRF parameter has normal distribution --> average and standard deviation of the parameter per group and condition in each ROI       

    T2DM_thr_avg_headers = cell(1,size(avg_parameters_T2DM_thr,2));                                       

    for l = 1:length(avg_std_indexes_T2DM_thr)
        T2DM_thr_avg_headers(l) = headers(avg_std_indexes_T2DM_thr(l));                                 % gets the headers of normal T2DM Thr HRF parameters
    end
    
    % Forms tables with the average and standard deviation for the T2DM Thr
    % HRF parameters with normal distribution as well as its headers
    T2DM_thr_avg_table = array2table(avg_parameters_T2DM_thr,'VariableNames', T2DM_thr_avg_headers);    
    T2DM_thr_std_table = array2table(std_parameters_T2DM_thr,'VariableNames', T2DM_thr_avg_headers);      
end


if isempty(avg_parameters_T2DM_sub)==0                                                                  % when a T2DM Sub HRF parameter has normal distribution --> average and standard deviation of the parameter per group and condition in each ROI      

    T2DM_sub_avg_headers = cell(1,size(avg_parameters_T2DM_sub,2));                                       

    for l = 1:length(avg_std_indexes_T2DM_sub)
        T2DM_sub_avg_headers(l) = headers(avg_std_indexes_T2DM_sub(l));                                 % gets the headers of normal T2DM Sub HRF parameters
    end

    % Forms tables with the average and standard deviation for the T2DM Sub
    % HRF parameters with normal distribution as well as its headers
    T2DM_sub_avg_table = array2table(avg_parameters_T2DM_sub,'VariableNames', T2DM_sub_avg_headers);    
    T2DM_sub_std_table = array2table(std_parameters_T2DM_sub,'VariableNames', T2DM_sub_avg_headers);      
end



% ««««««««««««««««««« Median and inter-quartile range »»»»»»»»»»»»»»»»»»»»


% -------------------------------- CNT -----------------------------------

if isempty(median_parameters_CNT_thr)==0                                                                            % when a CNT Thr HRF parameter has non-normal distribution --> median and interquartile range of the parameter per group and condition in each ROI    

    CNT_thr_median_headers = cell(1,size(median_parameters_CNT_thr,2));                                             % headers for the table

    for l = 1:length(median_iqr_indexes_CNT_thr)
        CNT_thr_median_headers(l) = headers(median_iqr_indexes_CNT_thr(l));                                         % gets the headers of non-normal CNT Thr HRF parameters
    end
    
    % Forms tables with the median and interquartile range for the CNT Thr
    % HRF parameters with non-normal distribution as well as its headers
    CNT_thr_median_table = array2table(median_parameters_CNT_thr,'VariableNames', CNT_thr_median_headers);          
    CNT_thr_iqr_table = array2table(iqr_parameters_CNT_thr,'VariableNames', CNT_thr_median_headers);                            
end

if isempty(median_parameters_CNT_sub)==0                                                                            % when a CNT Sub HRF parameter has non-normal distribution --> median and interquartile range of the parameter per group and condition in each ROI              

    CNT_sub_median_headers = cell(1,size(median_parameters_CNT_sub,2));                                            

    for l = 1:length(median_iqr_indexes_CNT_sub)
        CNT_sub_median_headers(l) = headers(median_iqr_indexes_CNT_sub(l));                                         % gets the headers of non-normal CNT Sub HRF parameters
    end
    
    % Forms tables with the median and interquartile range for the CNT Sub
    % HRF parameters with non-normal distribution as well as its headers
    CNT_sub_median_table = array2table(median_parameters_CNT_sub,'VariableNames', CNT_sub_median_headers);          
    CNT_sub_iqr_table = array2table(iqr_parameters_CNT_sub,'VariableNames', CNT_sub_median_headers);               
end


% -------------------------------- T2DM ----------------------------------


if isempty(median_parameters_T2DM_thr)==0                                                                           % when a T2DM Thr HRF parameter has non-normal distribution --> median and interquartile range of the parameter per group and condition in each ROI        

    T2DM_thr_median_headers = cell(1,size(median_parameters_T2DM_thr,2));                                          

    for l = 1:length(median_iqr_indexes_T2DM_thr)
        T2DM_thr_median_headers(l) = headers(median_iqr_indexes_T2DM_thr(l));                                       % gets the headers of non-normal T2DM Thr HRF parameters
    end

    % Forms tables with the median and interquartile range for the T2DM Thr
    % HRF parameters with non-normal distribution as well as its headers
    T2DM_thr_median_table = array2table(median_parameters_T2DM_thr,'VariableNames', T2DM_thr_median_headers);       
    T2DM_thr_iqr_table = array2table(iqr_parameters_T2DM_thr,'VariableNames', T2DM_thr_median_headers);       
end

if isempty(median_parameters_T2DM_sub)==0                                                                           % when a T2DM Sub HRF parameter has non-normal distribution --> median and interquartile range of the parameter per group and condition in each ROI        

    T2DM_sub_median_headers = cell(1,size(median_parameters_T2DM_sub,2));                                          

    for l = 1:length(median_iqr_indexes_T2DM_sub)
        T2DM_sub_median_headers(l) = headers(median_iqr_indexes_T2DM_sub(l));                                       % gets the headers of non-normal T2DM Sub HRF parameters
    end

    % Forms tables with the median and interquartile range for the T2DM Sub
    % HRF parameters with non-normal distribution as well as its headers
    T2DM_sub_median_table = array2table(median_parameters_T2DM_sub,'VariableNames', T2DM_sub_median_headers);       
    T2DM_sub_iqr_table = array2table(iqr_parameters_T2DM_sub,'VariableNames', T2DM_sub_median_headers);       
end