clear
close all
clc

%% Specify the model version: 'BNT', 'CDR', 'min', or 'max'
version = 'max';
% If this flag is set to one, it will regress out age, sex, and education effects, if not, it won't
flag_regress_out = 0;

%%
cd /data/dadmah/metame/NIFD_PLSpaper
addpath(genpath('/data/dadmah/resources/Pls'));
addpath(genpath('/data/datamah/in_vivo_Datasets/NIFD'));
addpath(genpath('/data/datamah/in_vivo_Datasets/NIFD/PLS'));

%% Read the descriptions in /data/dadmah/resources/Pls/plscmd/pls_analysis.m to see the details on PLS options
%figure
%set(gcf, 'PaperPosition', [0 0 10 12]);

%% Reading the data and making an overall FTD group with all 3 subgroups in it
Info=readtable('NIFD_MRI_Measures.csv','Delimiter',',');
%Info=readtable('Info_agematched.csv','Delimiter',',');
%Info(:,'sum_PPVT')=[];

Info.DXg(strcmp(Info.DX,'CON'))={'NC'};
Info.DXg(strcmp(Info.DX,'SV'))={'FTD'};
Info.DXg(strcmp(Info.DX,'BV'))={'FTD'};
Info.DXg(strcmp(Info.DX,'PNFA'))={'FTD'};

%% Remove participant with inconsistent data
Info(strcmp(Info.visit_id, '1_S_0020_1'), :) = [];

%% Imputing education variables
Info.EDUCATION(isnan(Info.EDUCATION)&strcmp(Info.DX,'CON'))=nanmean(Info.EDUCATION(strcmp(Info.DX,'CON')));
Info.EDUCATION(isnan(Info.EDUCATION)&strcmp(Info.DX,'SV'))=nanmean(Info.EDUCATION(strcmp(Info.DX,'SV')));
Info.EDUCATION(isnan(Info.EDUCATION)&strcmp(Info.DX,'BV'))=nanmean(Info.EDUCATION(strcmp(Info.DX,'BV')));
Info.EDUCATION(isnan(Info.EDUCATION)&strcmp(Info.DX,'PNFA'))=nanmean(Info.EDUCATION(strcmp(Info.DX,'PNFA')));

%% add total PPVT score
PPVT_indices = [35:38];
PPVTsum = Info{:,PPVT_indices};
sumPPVT = sum(PPVTsum, 2, 'includenan');
Info.sum_PPVT = sumPPVT;

%% pick behavioral variables to include
%13,14,16 CDR
%18 MMSE
%19,20,22,24 CVLT
%25,26 TMT
%28,29 digit span
%30,31 verbal fluency
%32 BNT
%35-38/253 PPVT

if strcmp(version, 'BNT')
    ind_beh = [32, 8:10];  % BNT only model
elseif strcmp(version, 'CDR')
    ind_beh = [13, 14, 16, 8:10];  % CDR only model
elseif strcmp(version, 'min')
    ind_beh = [13, 14, 16, 32, 8:10];  % Minimal model (CDR + BNT)
elseif strcmp(version, 'max')
    ind_beh = [13, 14, 16, 18, 19, 20, 22, 24:26, 28:32, 253, 8:10];  % Maximal model
else
    error('Invalid version selected. Choose ''BNT'', ''CDR'', ''min'', or ''max''.');
end

%% change negative behavioural values to NaN
data_subset = table2array(Info(:, ind_beh));
data_subset(data_subset < 0) = NaN;
Info(:, ind_beh) = array2table(data_subset);

Info(sum(isnan(table2array(Info(:, ind_beh))),2)>3,:)=[];
%% We're running the PLS on cross-sectional data, so here we just select and use one timepoint per participant
%pick baseline vs all longitudinal visits
Info.CLINICAL_LINKDATE = datetime(Info.CLINICAL_LINKDATE, 'InputFormat', 'MM-dd-yy');
Info = sortrows(Info, {'LONI_ID', 'CLINICAL_LINKDATE'});

[ia,ib,ic]=unique(Info.LONI_ID);
Info_Lng=Info;
Info=Info(ib,:);

[ia,ib,ic]=unique(Info.LONI_ID);
dx_grp=1;
data=Info(ib,:);
dx_grps=[{'FTD'}];
dx_grps_c=[{'k'},{'b'},{'r'}];
data.Sex = ismember(data.GENDER,2);

%separate controls
data_control=data(strcmp(data.DXg,'NC'),:);
data = data(ismember(data.DXg,dx_grps(dx_grp)),:);

%separate table for only chosen behavioural variables (controls vs FTD)
data_beh = data(:,ind_beh);
data_beh_mat = table2array(data_beh);

data_beh_control = data_control(:,ind_beh);
data_beh_mat_control = table2array(data_beh_control);

%% pick DBM/brain atrophy variables
ind_brain = [131:232];
data_brain = data(:,ind_brain);
data_brain_mat = table2array(data_brain);

%% remove NaNs
[r1,c1]=find(isnan(data_beh_mat));
[r2,c2]=find(isnan(data_brain_mat));
r=[r1;r2];

dx_included=data.DX(setdiff(1:size(data_brain_mat,1),r));
ids_included=data.LONI_ID(setdiff(1:size(data_brain_mat,1),r));
tabletrain=data(setdiff(1:size(data_brain_mat,1),r),:);

%remove NaN cases
data_beh_mat = data_beh_mat(setdiff(1:size(data_beh_mat,1),r),:);
data_brain_mat = data_brain_mat(setdiff(1:size(data_brain_mat,1),r),:);

%% Regress out effects of age, sex, and education based on healthy controls (if chosen above)
if flag_regress_out == 1
    %% Regress out effects of age, sex, and education based on healthy controls
    for i = 1:size(data_beh_mat, 2) - 3
        % Calculate the regression weights
        w = regress(data_beh_mat_control(:, i), ...
            [data_beh_mat_control(:, end-2:end) ones(size(data_beh_mat_control, 1), 1)]);
        % Regress out the effects from the full dataset
        data_beh_mat(:, i) = data_beh_mat(:, i) - ...
            [data_beh_mat(:, end-2:end) ones(size(data_beh_mat, 1), 1)] * w;
    end
    disp('Age, sex, and education effects have been regressed out.');
else
    disp('Skipping regression of age, sex, and education effects.');
end

%% Extracting the names of the variables

if strcmp(version, 'BNT')
    labels_x = {'BNT'; 'Age'; 'Sex';'Education'};
elseif strcmp(version, 'CDR')
    labels_x = {'CDR language'; 'CDR behavior'; 'CDR sum of boxes'; 'Age'; 'Sex';'Education'};
elseif strcmp(version, 'min')
    labels_x = {'CDR language'; 'CDR behavior'; 'CDR sum of boxes'; 'BNT'; 'Age'; 'Sex';'Education'};
elseif strcmp(version, 'max')
    labels_x = {'CDR language'; 'CDR behavior'; 'CDR sum of boxes'; 'MMSE'; 'CVLT total recall'; 'CVLT 30s delay'; 'CVLT 10m delay'; 'CVLT recognition'; 'MTMT time'; 'MTMT correct lines'; 'Forward digit span'; 'Backward digit span'; 'Verbal fluency (phon.)'; 'Verbal fluency (sem.)'; 'BNT'; 'PPVT'; 'Age'; 'Sex';'Education'};
else
    error('Invalid version selected. Choose ''BNT'', ''CDR'', ''min'', or ''max''.');
end

labels_y = data.Properties.VariableNames(ind_brain)';labels_y = strrep(labels_y,'_',' ')';

%% Zscoring the data
x = zscore(data_beh_mat);
y = zscore(data_brain_mat);

%% Partial Least Squares
%% Selecting PLS options and running analysis, check against the PLS descriptions

num_subj_lst = size(x,1);
num_cond = 1;
option.num_boot = 500;
option.num_perm = 500;
option.method = 3;
rng(12)
datamat_lst{1} = x;
option.stacked_behavdata = y;
result = pls_analysis(datamat_lst, num_subj_lst, num_cond, option);
datamat_lst{1} = y;
option.stacked_behavdata = x;
result_reverse = pls_analysis(datamat_lst, num_subj_lst, num_cond, option);

%% save results as a matlab file
 if strcmp(version, 'BNT')
    save results_PLS_BNT
 elseif strcmp(version, 'CDR')
     save results_PLS_CDR
 elseif strcmp(version, 'min')
     save results_PLS_min
 elseif strcmp(version, 'max')
     save results_PLS_max
 else

 end

%% print amount of covariance explained of each LV

(result.s.^2)./sum(result.s.^2)

%% save csv files for LVs
% save csv file with all significantly contributing brain regions for each LV
if strcmp(version, 'max')

    for comp=1:4

        U_values = result.boot_result.ulcorr(:, comp);
        L_values = result.boot_result.llcorr(:, comp);
        C_values = double(result.boot_result.orig_corr(:, comp));
        labels = labels_y;

        test = table(U_values, L_values, C_values);

        test.U = U_values-C_values;
        test.L = test.C_values-L_values;
        test.B = data.Properties.VariableNames(ind_brain)';labels_y = strrep(labels_y,'_',' ')';

        test.X = -1*C_values+test.U;
        test.Y = -1*C_values-test.L;

        condition_1 = test.X> 0 & test.Y > 0;
        condition_2 = test.X < 0 & test.Y < 0;
        new_column = zeros(size(test.X));
        new_column(condition_1) = 1;
        new_column(condition_2) = 2;
        test = addvars(test, new_column, 'Before', 1, 'NewVariableName', 'sig');

        test.mean=(test.X + test.Y)/2;

        test = test(test.sig~=0,:);

        switch comp
            case 1
                writetable(test,'CI_brain1.csv');
            case 2
                writetable(test,'CI_brain2.csv');
            case 3
                writetable(test,'CI_brain3.csv');
            case 4
                writetable(test,'CI_brain4.csv');
        end
    end

else

end

% save csv file with all significantly contributing cognition scores for each LV
if strcmp(version, 'max')
    for comp=1:4

        U_values = result_reverse.boot_result.ulcorr(:, comp);
        L_values = result_reverse.boot_result.llcorr(:, comp);
        C_values = double(result_reverse.boot_result.orig_corr(:, comp));

        test = table(U_values, L_values, C_values, labels_x);

        test.U = U_values-C_values;
        test.L = test.C_values-L_values;

        test.upper = -1*test.C_values+test.U;
        test.lower = -1*test.C_values-test.L;

        test.beh = (sign(result_reverse.boot_result.llcorr(:,comp).*result_reverse.boot_result.ulcorr(:,comp)))>0;

        condition_1 = test.upper> 0 & test.lower > 0;
        condition_2 = test.upper < 0 & test.lower < 0;
        new_column = zeros(size(test.upper));
        new_column(condition_1) = 1;
        new_column(condition_2) = 2;
        test = addvars(test, new_column, 'Before', 1, 'NewVariableName', 'sig');

        test.mean=(test.upper + test.lower)/2;

        switch comp
            case 1
                writetable(test,'CI_beh1.csv');
            case 2
                writetable(test,'CI_beh2.csv');
            case 3
                writetable(test,'CI_beh3.csv');
            case 4
                writetable(test,'CI_beh4.csv');
        end
    end

else

end