%% this is where we run the FTD subtype classification model
%% in a cross-validated loop
%% the resulting table includes accuracies, sensitivity, specificity, and balanced accuracies
%% for all model inputs that are selected (brain scores only, cognition scores only, or both)
%% and for longitudinal validation

%% Specify the model version: 'max', 'min', 'BNT', 'CDR'

    version = 'max';

%% pick which inputs you want to run (both = brain + cognition patters, neur = only brain, cogn = only cognition)

modelList = {'both','neur','cogn'};

%set seed for reproducibility
rng(82);

if strcmp(version, 'min')
    load results_PLS_min_lng
elseif strcmp(version, 'max')
    load results_PLS_max_lng
elseif strcmp(version, 'CDR')
    load results_PLS_CDR_lng
elseif strcmp(version, 'BNT')
    load results_PLS_BNT_lng
else
end

n_comp = 4; % number of components to use (based on how many were significant)

%% add table with all participants that were included in the maximal model (to compare performance in that sample)

tablemax=readtable('table_maxsample.csv');

%% crossvalidated prediction analysis

temp = cell(numel(modelList), 1);

for i = 1:numel(modelList)  
    models = modelList{i};  % Select the current model

    % Select data based on the model type
    switch models
        case 'both'
            X_all = [result.usc(:, 1:n_comp), result.vsc(:, 1:n_comp)];
            X_cs = table2array(cs(:, 3:end));
            X_lng = table2array(lng(:, 3:end));
            
        case 'neur'
            X_all = result.vsc(:, 1:n_comp);
            X_cs = table2array(cs(:, 7:end));
            X_lng = table2array(lng(:, 7:end));
            
        case 'cogn'
            X_all = result.usc(:, 1:n_comp);
            X_cs = table2array(cs(:, 3:6));
            X_lng = table2array(lng(:, 3:6));
            
        otherwise
            error('Invalid model type: %s', models);
    end

    %index DX groups
    Y_all (strcmp(dx_included,'BV'),1)=1;
    Y_all (strcmp(dx_included,'SV'),1)=2;
    Y_all (strcmp(dx_included,'PNFA'),1)=3;

    Y_cs (strcmp(cs.DX,'BV'),1)=1;
    Y_cs (strcmp(cs.DX,'SV'),1)=2;
    Y_cs (strcmp(cs.DX,'PNFA'),1)=3;

    Y_lng (strcmp(lng.DX,'BV'),1)=1;
    Y_lng (strcmp(lng.DX,'SV'),1)=2;
    Y_lng (strcmp(lng.DX,'PNFA'),1)=3;

    %index for CDR matched (lower severity) participants
    ind_cdr_lt1 = tabletrain.CDR_TOT<1.5;

    %index for participants in max model (to compare with min model)
    [~, idx1] = ismember(tablemax.LONI_ID, tabletrain.LONI_ID);idx1=idx1(idx1>0);

    %prediction loop

    for r = 1:100 % do it 100 times to get an accurate estimate
        r
        indices = crossvalind('Kfold',size(X_all,1),10); % 10-fold cross validation indices
        for k = 1:10
            X_train = X_all(indices ~= k,:);
            X_test =  X_all(indices == k,:);
            Y_train = Y_all(indices ~= k);
            Y_test = Y_all(indices == k);

            Mdl = fitcensemble(X_train,Y_train,'Method','Bag','Learners','discriminant');
            Y_hat(indices == k,1) = predict(Mdl,X_test); % train on 90%, test on the rest, y_hat is the vector of model predictions
            ind_lng=find(ismember(cs.LONI_ID,tabletrain.LONI_ID(indices == k)));
            Y_hat_cs(ind_lng,1)=predict(Mdl,X_cs(ind_lng,:));
            Y_hat_lng(ind_lng,1)=predict(Mdl,X_lng(ind_lng,:));

        end

        % save prediction accuracy

        Accuracy(r,1) = sum(Y_all==Y_hat)/size(Y_all,1);
        Accuracy_cs(r,1) = sum(Y_cs==Y_hat_cs)/size(Y_cs,1);
        Accuracy_lng(r,1) = sum(Y_lng==Y_hat_lng)/size(Y_lng,1);
        Accuracy_CDR_lt1(r,1) = sum(Y_all(ind_cdr_lt1) == Y_hat(ind_cdr_lt1))/size(Y_all(ind_cdr_lt1),1);
        Accuracy_CDR_max(r,1) = sum(Y_all(idx1) == Y_hat(idx1))/size(Y_all(idx1),1);

        Specificity_BV(r,1) = sum((Y_all~=1)&(Y_hat~=1))/sum(Y_all~=1);
        Specificity_SV(r,1) = sum((Y_all~=2)&(Y_hat~=2))/sum(Y_all~=2);
        Specificity_PNFA(r,1) = sum((Y_all~=3)&(Y_hat~=3))/sum(Y_all~=3);

        Sensitivity_BV(r,1) = sum((Y_all==1)&(Y_hat==1))/sum(Y_all==1);
        Sensitivity_SV(r,1) = sum((Y_all==2)&(Y_hat==2))/sum(Y_all==2);
        Sensitivity_PNFA(r,1) = sum((Y_all==3)&(Y_hat==3))/sum(Y_all==3);

        Specificity_BV_cs(r,1) = sum((Y_cs~=1)&(Y_hat_cs~=1))/sum(Y_cs~=1);
        Specificity_SV_cs(r,1) = sum((Y_cs~=2)&(Y_hat_cs~=2))/sum(Y_cs~=2);
        Specificity_PNFA_cs(r,1) = sum((Y_cs~=3)&(Y_hat_cs~=3))/sum(Y_cs~=3);

        Sensitivity_BV_cs(r,1) = sum((Y_cs==1)&(Y_hat_cs==1))/sum(Y_cs==1);
        Sensitivity_SV_cs(r,1) = sum((Y_cs==2)&(Y_hat_cs==2))/sum(Y_cs==2);
        Sensitivity_PNFA_cs(r,1) = sum((Y_cs==3)&(Y_hat_cs==3))/sum(Y_cs==3);

        Specificity_BV_lng(r,1) = sum((Y_lng~=1)&(Y_hat_lng~=1))/sum(Y_lng~=1);
        Specificity_SV_lng(r,1) = sum((Y_lng~=2)&(Y_hat_lng~=2))/sum(Y_lng~=2);
        Specificity_PNFA_lng(r,1) = sum((Y_lng~=3)&(Y_hat_lng~=3))/sum(Y_lng~=3);

        Sensitivity_BV_lng(r,1) = sum((Y_lng==1)&(Y_hat_lng==1))/sum(Y_lng==1);
        Sensitivity_SV_lng(r,1) = sum((Y_lng==2)&(Y_hat_lng==2))/sum(Y_lng==2);
        Sensitivity_PNFA_lng(r,1) = sum((Y_lng==3)&(Y_hat_lng==3))/sum(Y_lng==3);

        Specificity_BV_lt1(r,1) = sum((Y_all(ind_cdr_lt1)~=1)&(Y_hat(ind_cdr_lt1)~=1))/sum(Y_all(ind_cdr_lt1)~=1);
        Specificity_SV_lt1(r,1) = sum((Y_all(ind_cdr_lt1)~=2)&(Y_hat(ind_cdr_lt1)~=2))/sum(Y_all(ind_cdr_lt1)~=2);
        Specificity_PNFA_lt1(r,1) = sum((Y_all(ind_cdr_lt1)~=3)&(Y_hat(ind_cdr_lt1)~=3))/sum(Y_all(ind_cdr_lt1)~=3);

        Sensitivity_BV_lt1(r,1) = sum((Y_all(ind_cdr_lt1)==1)&(Y_hat(ind_cdr_lt1)==1))/sum(Y_all(ind_cdr_lt1)==1);
        Sensitivity_SV_lt1(r,1) = sum((Y_all(ind_cdr_lt1)==2)&(Y_hat(ind_cdr_lt1)==2))/sum(Y_all(ind_cdr_lt1)==2);
        Sensitivity_PNFA_lt1(r,1) = sum((Y_all(ind_cdr_lt1)==3)&(Y_hat(ind_cdr_lt1)==3))/sum(Y_all(ind_cdr_lt1)==3);

        Specificity_BV_max(r,1) = sum((Y_all(idx1)~=1)&(Y_hat(idx1)~=1))/sum(Y_all(idx1)~=1);
        Specificity_SV_max(r,1) = sum((Y_all(idx1)~=2)&(Y_hat(idx1)~=2))/sum(Y_all(idx1)~=2);
        Specificity_PNFA_max(r,1) = sum((Y_all(idx1)~=3)&(Y_hat(idx1)~=3))/sum(Y_all(idx1)~=3);

        Sensitivity_BV_max(r,1) = sum((Y_all(idx1)==1)&(Y_hat(idx1)==1))/sum(Y_all(idx1)==1);
        Sensitivity_SV_max(r,1) = sum((Y_all(idx1)==2)&(Y_hat(idx1)==2))/sum(Y_all(idx1)==2);
        Sensitivity_PNFA_max(r,1) = sum((Y_all(idx1)==3)&(Y_hat(idx1)==3))/sum(Y_all(idx1)==3);

        %display accuracy in current loop
        sum(Y_all==Y_hat)/size(Y_all,1)

    end

   % save all results in a table
Results = table;   
Results.model = models;

% Overall accuracy metrics (cross-sectional)
Results.Overall_Accuracy = [round(mean(Accuracy) * 100, 2)];
Results.Overall_sd = [round(std(Accuracy) * 100, 2)];% Convert to percentage
Results.Sensitivity = struct('bvFTD', [round(mean(Sensitivity_BV) * 100, 2), round(std(Sensitivity_BV) * 100, 2)], ...
                             'svPPA', [round(mean(Sensitivity_SV) * 100, 2), round(std(Sensitivity_SV) * 100, 2)], ...
                             'nfvPPA', [round(mean(Sensitivity_PNFA) * 100, 2), round(std(Sensitivity_PNFA) * 100, 2)]);

Results.Specificity = struct('bvFTD', [round(mean(Specificity_BV) * 100, 2), round(std(Specificity_BV) * 100, 2)], ...
                             'svPPA', [round(mean(Specificity_SV) * 100, 2), round(std(Specificity_SV) * 100, 2)], ...
                             'nfvPPA', [round(mean(Specificity_PNFA) * 100, 2), round(std(Specificity_PNFA) * 100, 2)]);

% Balanced accuracies
Results.Balanced_Accuracy = struct('bvFTD', round(((mean(Specificity_BV) + mean(Sensitivity_BV)) / 2) * 100, 2), ...
                                   'svPPA', round(((mean(Specificity_SV) + mean(Sensitivity_SV)) / 2) * 100, 2), ...
                                   'nfvPPA', round(((mean(Specificity_PNFA) + mean(Sensitivity_PNFA)) / 2) * 100, 2));

% Longitudinal accuracy metrics
Results.Longitudinal_Accuracy = [round(mean(Accuracy_lng) * 100, 2)]; % Convert to percentage
Results.Longitudinal_sd = [round(std(Accuracy_lng) * 100, 2)]; % Convert to percentage

Results.Longitudinal_Sensitivity = struct(...
    'bvFTD', [round(mean(Sensitivity_BV_lng) * 100, 2), round(std(Sensitivity_BV_lng) * 100, 2)], ...
    'svPPA', [round(mean(Sensitivity_SV_lng) * 100, 2), round(std(Sensitivity_SV_lng) * 100, 2)], ...
    'nfvPPA', [round(mean(Sensitivity_PNFA_lng) * 100, 2), round(std(Sensitivity_PNFA_lng) * 100, 2)]);

Results.Longitudinal_Specificity = struct(...
    'bvFTD', [round(mean(Specificity_BV_lng) * 100, 2), round(std(Specificity_BV_lng) * 100, 2)], ...
    'svPPA', [round(mean(Specificity_SV_lng) * 100, 2), round(std(Specificity_SV_lng) * 100, 2)], ...
    'nfvPPA', [round(mean(Specificity_PNFA_lng) * 100, 2), round(std(Specificity_PNFA_lng) * 100, 2)]);

% Balanced accuracies
Results.Longitudinal_Balanced_Accuracy = struct(...
    'bvFTD', round(((mean(Specificity_BV_lng) + mean(Sensitivity_BV_lng)) / 2) * 100, 2), ...
    'svPPA', round(((mean(Specificity_SV_lng) + mean(Sensitivity_SV_lng)) / 2) * 100, 2), ...
    'nfvPPA', round(((mean(Specificity_PNFA_lng) + mean(Sensitivity_PNFA_lng)) / 2) * 100, 2));

% CDR matched sample accuracy
Results.CDR_Matched_Accuracy = [round(mean(Accuracy_CDR_lt1) * 100, 2)];
Results.CDR_Matched_sd = [round(std(Accuracy_CDR_lt1) * 100, 2)];

Results.CDR_Matched_Sensitivity = struct(...
    'bvFTD', [round(mean(Sensitivity_BV_lt1) * 100, 2), round(std(Sensitivity_BV_lt1) * 100, 2)], ...
    'svPPA', [round(mean(Sensitivity_SV_lt1) * 100, 2), round(std(Sensitivity_SV_lt1) * 100, 2)], ...
    'nfvPPA', [round(mean(Sensitivity_PNFA_lt1) * 100, 2), round(std(Sensitivity_PNFA_lt1) * 100, 2)]);

Results.CDR_Matched_Specificity = struct(...
    'bvFTD', [round(mean(Specificity_BV_lt1) * 100, 2), round(std(Specificity_BV_lt1) * 100, 2)], ...
    'svPPA', [round(mean(Specificity_SV_lt1) * 100, 2), round(std(Specificity_SV_lt1) * 100, 2)], ...
    'nfvPPA', [round(mean(Specificity_PNFA_lt1) * 100, 2), round(std(Specificity_PNFA_lt1) * 100, 2)]);

Results.CDR_Matched_Balanced_Accuracy = struct(...
    'bvFTD', round(((mean(Specificity_BV_lt1) + mean(Sensitivity_BV_lt1)) / 2) * 100, 2), ...
    'svPPA', round(((mean(Specificity_SV_lt1) + mean(Sensitivity_SV_lt1)) / 2) * 100, 2), ...
    'nfvPPA', round(((mean(Specificity_PNFA_lt1) + mean(Sensitivity_PNFA_lt1)) / 2) * 100, 2));

% Max sample accuracy
Results.Max_Sample_Accuracy = [round(mean(Accuracy_CDR_max) * 100, 2)];
Results.Max_Sample_sd = [round(std(Accuracy_CDR_max) * 100, 2)];

Results.Max_Sample_Sensitivity = struct(...
    'bvFTD', [round(mean(Sensitivity_BV_max) * 100, 2), round(std(Sensitivity_BV_max) * 100, 2)], ...
    'svPPA', [round(mean(Sensitivity_SV_max) * 100, 2), round(std(Sensitivity_SV_max) * 100, 2)], ...
    'nfvPPA', [round(mean(Sensitivity_PNFA_max) * 100, 2), round(std(Sensitivity_PNFA_max) * 100, 2)]);

Results.Max_Sample_Specificity = struct(...
    'bvFTD', [round(mean(Specificity_BV_max) * 100, 2), round(std(Specificity_BV_max) * 100, 2)], ...
    'svPPA', [round(mean(Specificity_SV_max) * 100, 2), round(std(Specificity_SV_max) * 100, 2)], ...
    'nfvPPA', [round(mean(Specificity_PNFA_max) * 100, 2), round(std(Specificity_PNFA_max) * 100, 2)]);

Results.Max_Sample_Balanced_Accuracy = struct(...
    'bvFTD', round(((mean(Specificity_BV_max) + mean(Sensitivity_BV_max)) / 2) * 100, 2), ...
    'svPPA', round(((mean(Specificity_SV_max) + mean(Sensitivity_SV_max)) / 2) * 100, 2), ...
    'nfvPPA', round(((mean(Specificity_PNFA_max) + mean(Sensitivity_PNFA_max)) / 2) * 100, 2));

temp{i} = Results;

end

% create more readable table

table1 = temp{1};
table2 = temp{2};
table3 = temp{3};
mergedTable = vertcat(table1, table2, table3);

summary = table();
for i = 1:height(mergedTable)
    tempStruct = struct();
    fields = mergedTable.Properties.VariableNames; 
    for j = 1:length(fields)
        fieldName = fields{j};
        if isstruct(mergedTable.(fieldName)(i)) 
            subFields = fieldnames(mergedTable.(fieldName)(i)); 
            for k = 1:length(subFields)
                subFieldName = subFields{k};
               tempStruct.([fieldName '_' subFieldName]) = mergedTable.(fieldName)(i).(subFieldName);
            end
        else
            tempStruct.(fieldName) = mergedTable.(fieldName)(i); 
        end
    end
    summary = [summary; struct2table(tempStruct)];
end

summary.cohort(1:3,:)="UCSF";
summary = movevars(summary, "cohort", "Before", "model");
summary.version(1:3,:)="min";
summary = movevars(summary, "version", "Before", "cohort");

% save results
% allresults=readtable('predictionresults.csv');
% writetable(summary, 'temp.csv');
% temp=readtable('temp.csv');
% temp=[allresults;temp];
% writetable(temp, 'predictionresults.csv');

%% print numbers of subjects in each model

[unique_vals, ~, idx] = unique(dx_included);
accumarray(idx, 1)

[unique_vals, ~, idx] = unique(dxlng_included);
accumarray(idx, 1)
