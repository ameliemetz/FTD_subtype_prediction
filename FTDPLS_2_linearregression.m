%% this code is used to test whether the relationship in brain patterns and associated cognition patterns
%% based on latent variables (resulting from PLS analysis)
%% differs between FTD subtypes

%% we are running a linear regression model for each significant latent variable
%% cognition ~ brain * DX

%% Specify the model version: 'max', 'min', 'BNT', 'CDR'

    version = 'max';

%either run PLS script right before or load saved results and pick number of components to analyse

if strcmp(version, 'min')
load results_PLS_min
c3=1:3;
else
load results_PLS_max
c3=1:4;
end

%% run linear regression model relating brain and cognition scores in each FTD subtype
% save results for intercepts and interactions

results_table = table;

for comp=c3

    new=table;
    new.usc=result.usc(:,comp);
    new.vsc=result.vsc(:,comp);
    new.dx=dx_included;
%vs BV
    new = sortrows(new, {'dx'}, "ascend");
    mylm1=fitlm(new,'usc~vsc*dx')
%vs SV
    new = sortrows(new, {'dx'}, "descend");
    mylm2=fitlm(new,'usc~vsc*dx')

    stats1 = mylm1.Coefficients;
    stats2 = mylm2.Coefficients;

    temp_table = table;
    temp_table.comp = repmat(comp, size(stats1, 1) + size(stats2, 1), 1);  
    temp_table.Term = [stats1.Properties.RowNames; stats2.Properties.RowNames];
    temp_table.Estimate = [stats1.Estimate; stats2.Estimate];
    temp_table.tStat = [stats1.tStat; stats2.tStat];
    temp_table.pValue = [stats1.pValue; stats2.pValue];

    results_table = [results_table; temp_table];

end

if strcmp(version, 'min')
writetable(results_table,'Regression_results_min.csv');
else
writetable(results_table,'Regression_results_max.csv');
end
