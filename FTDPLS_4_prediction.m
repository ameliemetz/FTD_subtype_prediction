clear
close all
clc

%% Specify the model version: 'BNT', 'CDR', 'min', or 'max'
version = 'max';

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

n_comp = 4; % number of components to use
tablemax=readtable('tabletrain.csv');

modelList = {'both','neur','cogn'};

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
        % prediction accuracy, the %of correct predictions
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

        sum(Y_all==Y_hat)/size(Y_all,1)
    end

    % %% overall accuracy, crosssectional
    % [mean(Accuracy) std(Accuracy)]
    % [mean(Sensitivity_BV) std(Sensitivity_BV)]
    % [mean(Sensitivity_SV) std(Sensitivity_SV)]
    % [mean(Sensitivity_PNFA) std(Sensitivity_PNFA)]
    %
    % [mean(Specificity_BV) std(Specificity_BV)]
    % [mean(Specificity_SV) std(Specificity_SV)]
    % [mean(Specificity_PNFA) std(Specificity_PNFA)]
    %
    % %balanced accuracies
    % (mean(Specificity_BV) + mean(Sensitivity_BV))/2
    % (mean(Specificity_SV) + mean(Sensitivity_SV))/2
    % (mean(Specificity_PNFA) + mean(Sensitivity_PNFA))/2
    %
    % %% accuracy crosssectional, sample with lng data
    % [mean(Accuracy_cs) std(Accuracy_cs)]
    % [mean(Sensitivity_BV_cs) std(Sensitivity_BV_cs)]
    % [mean(Sensitivity_SV_cs) std(Sensitivity_SV_cs)]
    % [mean(Sensitivity_PNFA_cs) std(Sensitivity_PNFA_cs)]
    %
    % [mean(Specificity_BV_cs) std(Specificity_BV_cs)]
    % [mean(Specificity_SV_cs) std(Specificity_SV_cs)]
    % [mean(Specificity_PNFA_cs) std(Specificity_PNFA_cs)]
    %
    % %% accuracy longitudinal
    % [mean(Accuracy_lng) std(Accuracy_lng)]
    %
    % [mean(Sensitivity_BV_lng) std(Sensitivity_BV_lng)]
    % [mean(Sensitivity_SV_lng) std(Sensitivity_SV_lng)]
    % [mean(Sensitivity_PNFA_lng) std(Sensitivity_PNFA_lng)]
    %
    % [mean(Specificity_BV_lng) std(Specificity_BV_lng)]
    % [mean(Specificity_SV_lng) std(Specificity_SV_lng)]
    % [mean(Specificity_PNFA_lng) std(Specificity_PNFA_lng)]
    %
    % %% accuarcy in CDR matched sample
    % [mean(Accuracy_CDR_lt1) std(Accuracy_CDR_lt1)]
    %
    % [mean(Sensitivity_BV_lt1) std(Sensitivity_BV_lt1)]
    % [mean(Sensitivity_SV_lt1) std(Sensitivity_SV_lt1)]
    % [mean(Sensitivity_PNFA_lt1) std(Sensitivity_PNFA_lt1)]
    %
    % [mean(Specificity_BV_lt1) std(Specificity_BV_lt1)]
    % [mean(Specificity_SV_lt1) std(Specificity_SV_lt1)]
    % [mean(Specificity_PNFA_lt1) std(Specificity_PNFA_lt1)]
    %
    % %% accuracy min in sample of max model
    % [mean(Accuracy_CDR_max) std(Accuracy_CDR_max)]
    %
    % [mean(Sensitivity_BV_max) std(Sensitivity_BV_max)]
    % [mean(Sensitivity_SV_max) std(Sensitivity_SV_max)]
    % [mean(Sensitivity_PNFA_max) std(Sensitivity_PNFA_max)]
    %
    % [mean(Specificity_BV_max) std(Specificity_BV_max)]
    % [mean(Specificity_SV_max) std(Specificity_SV_max)]
    % [mean(Specificity_PNFA_max) std(Specificity_PNFA_max)]

    % Calculate overall accuracy and relevant metrics
    Results = table;   
    Results.model = models;

    % Overall accuracy (cross-sectional)
    Results.Overall_Accuracy = [mean(Accuracy), std(Accuracy)];
    Results.Sensitivity_BV = [mean(Sensitivity_BV), std(Sensitivity_BV)];
    Results.Sensitivity_SV = [mean(Sensitivity_SV), std(Sensitivity_SV)];
    Results.Sensitivity_PNFA = [mean(Sensitivity_PNFA), std(Sensitivity_PNFA)];
    Results.Specificity_BV = [mean(Specificity_BV), std(Specificity_BV)];
    Results.Specificity_SV = [mean(Specificity_SV), std(Specificity_SV)];
    Results.Specificity_PNFA = [mean(Specificity_PNFA), std(Specificity_PNFA)];

    % Balanced accuracies
    Results.Balanced_Accuracy_BV = (mean(Specificity_BV) + mean(Sensitivity_BV)) / 2;
    Results.Balanced_Accuracy_SV = (mean(Specificity_SV) + mean(Sensitivity_SV)) / 2;
    Results.Balanced_Accuracy_PNFA = (mean(Specificity_PNFA) + mean(Sensitivity_PNFA)) / 2;

    % Accuracy (cross-sectional, sample with longitudinal data)
    Results.Accuracy_cs = [mean(Accuracy_cs), std(Accuracy_cs)];
    Results.Sensitivity_BV_cs = [mean(Sensitivity_BV_cs), std(Sensitivity_BV_cs)];
    Results.Sensitivity_SV_cs = [mean(Sensitivity_SV_cs), std(Sensitivity_SV_cs)];
    Results.Sensitivity_PNFA_cs = [mean(Sensitivity_PNFA_cs), std(Sensitivity_PNFA_cs)];
    Results.Specificity_BV_cs = [mean(Specificity_BV_cs), std(Specificity_BV_cs)];
    Results.Specificity_SV_cs = [mean(Specificity_SV_cs), std(Specificity_SV_cs)];
    Results.Specificity_PNFA_cs = [mean(Specificity_PNFA_cs), std(Specificity_PNFA_cs)];

    % Accuracy (longitudinal)
    Results.Accuracy_lng = [mean(Accuracy_lng), std(Accuracy_lng)];
    Results.Sensitivity_BV_lng = [mean(Sensitivity_BV_lng), std(Sensitivity_BV_lng)];
    Results.Sensitivity_SV_lng = [mean(Sensitivity_SV_lng), std(Sensitivity_SV_lng)];
    Results.Sensitivity_PNFA_lng = [mean(Sensitivity_PNFA_lng), std(Sensitivity_PNFA_lng)];
    Results.Specificity_BV_lng = [mean(Specificity_BV_lng), std(Specificity_BV_lng)];
    Results.Specificity_SV_lng = [mean(Specificity_SV_lng), std(Specificity_SV_lng)];
    Results.Specificity_PNFA_lng = [mean(Specificity_PNFA_lng), std(Specificity_PNFA_lng)];

    % Accuracy (CDR matched sample)
    Results.Accuracy_CDR_lt1 = [mean(Accuracy_CDR_lt1), std(Accuracy_CDR_lt1)];
    Results.Sensitivity_BV_lt1 = [mean(Sensitivity_BV_lt1), std(Sensitivity_BV_lt1)];
    Results.Sensitivity_SV_lt1 = [mean(Sensitivity_SV_lt1), std(Sensitivity_SV_lt1)];
    Results.Sensitivity_PNFA_lt1 = [mean(Sensitivity_PNFA_lt1), std(Sensitivity_PNFA_lt1)];
    Results.Specificity_BV_lt1 = [mean(Specificity_BV_lt1), std(Specificity_BV_lt1)];
    Results.Specificity_SV_lt1 = [mean(Specificity_SV_lt1), std(Specificity_SV_lt1)];
    Results.Specificity_PNFA_lt1 = [mean(Specificity_PNFA_lt1), std(Specificity_PNFA_lt1)];

    % Accuracy (max sample)
    Results.Accuracy_CDR_max = [mean(Accuracy_CDR_max), std(Accuracy_CDR_max)];
    Results.Sensitivity_BV_max = [mean(Sensitivity_BV_max), std(Sensitivity_BV_max)];
    Results.Sensitivity_SV_max = [mean(Sensitivity_SV_max), std(Sensitivity_SV_max)];
    Results.Sensitivity_PNFA_max = [mean(Sensitivity_PNFA_max), std(Sensitivity_PNFA_max)];
    Results.Specificity_BV_max = [mean(Specificity_BV_max), std(Specificity_BV_max)];
    Results.Specificity_SV_max = [mean(Specificity_SV_max), std(Specificity_SV_max)];
    Results.Specificity_PNFA_max = [mean(Specificity_PNFA_max), std(Specificity_PNFA_max)];

    %% Save the results
     temp{i} = Results;
end

finalTable = temp{1};
for r = 2:length(temp)
    finalTable = outerjoin(finalTable, temp{r}, 'MergeKeys', true);
end

switch version
    case 'max'
        writetable(summary, 'predictionresults_max.csv');
    case 'min'
        writetable(summary, 'predictionresults_min.csv');
    case 'CDR'
        writetable(summary, 'predictionresults_CDR.csv');
    case 'BNT'
        writetable(summary, 'predictionresults_BNT.csv');
end

%% print numbers of subjects in each model
[unique_vals, ~, idx] = unique(dx_included);
accumarray(idx, 1)

[unique_vals, ~, idx] = unique(dxlng_included);
accumarray(idx, 1)