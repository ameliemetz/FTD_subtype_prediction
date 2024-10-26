clear
close all
clc

%% longitudinal progression of brain and behavior patterns

load results_PLS_max_lng

%get time between baseline and followup
time = tablelng.TimefromBaseline(idx2,:);
%calculate delta (followup-baseline)/time
delta = table((Xlng(:,1:8)- Xbase(:,1:8)));
delta = delta.Var1(:,1:8) ./ time;

tbl = table(tablelng.LONI_ID(idx2, :), dxlng_included(idx2, :),Xbase(:,1:8),Xlng(:,1:8), time, delta(:,1:8), tablelng.Sex(idx2, :),tablelng.AGE(idx2, :), ...
    'VariableNames',{'ID','DX','base','lng','time','delta','sex','age'});

%% delta between DX (ttests)
for i=1:8
    i
    tbl.d = tbl.delta(:,i);
    x = tbl.d(strcmp(tbl.DX,'BV'),:);
    y = tbl.d(strcmp(tbl.DX,'PNFA'),:);
    z = tbl.d(strcmp(tbl.DX,'SV'),:);
    [h,p,ci,stats] = ttest2(x,y)
    [h,p,ci,stats] = ttest2(x,z)
    [h,p,ci,stats] = ttest2(y,z)
end

%% ttest to compare each LV in all participants
[h,p] = ttest(Xbase(:,1:4),Xlng(:,1:4)) %compare base and follow up behavioral scores
[h,p] = ttest(Xbase(:,5:8),Xlng(:,5:8)) %compare base and follow up brain scores

%% write tables to put them into R for the boxplots
writetable(tbl,'delta.csv')
writetable(cs,'compare_base.csv')
writetable(lng,'compare_lng.csv')