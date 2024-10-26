clear
close all
clc

version = 'max';

if strcmp(version, 'min')
    load results_PLS_min
elseif strcmp(version, 'max')
    load results_PLS_max
elseif strcmp(version, 'CDR')
    load results_PLS_CDR
elseif strcmp(version, 'BNT')
    load results_PLS_BNT
else
end

%% pick timepoint closest to one year

for i=1:size(Info_Lng,1)
    tmp=cell2mat(Info_Lng.DOB(i));
    Info_Lng.Age_exact(i,1)=(datenum(Info_Lng.CLINICAL_LINKDATE(i))-datenum(strcat(tmp(1:5),'-19',tmp(7:8)),'mm-dd-yy'))/365;
end

ids = unique(Info_Lng.LONI_ID);
for i = 1:numel(ids)
    ind = find(strcmp(ids(i), Info_Lng.LONI_ID));
    Info_Lng.TimefromBaseline(ind) = Info_Lng.Age_exact(ind) - min(Info_Lng.Age_exact(ind));
    Info_Lng.N_visits(ind) = numel(ind);
end

% determine which timepoint has all data available
Lng = Info_Lng(Info_Lng.TimefromBaseline~=0,:);
rows_with_nan = false(height(Lng), 1);
for col = ind_beh
    rows_with_nan = rows_with_nan | isnan(Lng{:, col});
end
Lng(rows_with_nan, :) = [];

Lng.Keep = zeros(height(Lng), 1);

%pick closest to one year
unique_ids = unique(Lng.LONI_ID);

for i = 1:length(unique_ids)
    ind = find(strcmp(Lng.LONI_ID, unique_ids{i}));
    if length(ind) == 1
        Lng.Keep(ind) = 1;
    else
        [~, min_idx] = min(abs(Lng.TimefromBaseline(ind) - 1));
        Lng.Keep(ind(min_idx)) = 1;
    end
end

datalng=Lng(Lng.Keep == 1,:);
database=Info;
datalng.Sex = ismember(datalng.GENDER,2);

datalng_control=datalng(strcmp(datalng.DXg,'NC'),:);
datalng = datalng(ismember(datalng.DXg,dx_grps(dx_grp)),:);

%% create individual brain and cognitive scores for each longitudinal visit

datalng_beh = datalng(:,ind_beh);
datalng_beh_mat = table2array(datalng_beh);

ind_brain = [131:232];
datalng_brain = datalng(:,ind_brain); % indices for the brain data you want to include (check the table)
datalng_brain_mat = table2array(datalng_brain);

%% remove the NaN cases
[r1,c1]=find(isnan(datalng_beh_mat));
[r2,c2]=find(isnan(datalng_brain_mat));
r=[r1;r2];
dxlng_included=datalng.DX(setdiff(1:size(datalng_brain_mat,1),r));
tablelng=datalng(setdiff(1:size(datalng_brain_mat,1),r),:);
datalng_beh_mat = datalng_beh_mat(setdiff(1:size(datalng_beh_mat,1),r),:);
datalng_brain_mat = datalng_brain_mat(setdiff(1:size(datalng_brain_mat,1),r),:);

%% Zscoring the data
%% this is equivalent of zscoring data_beh_mat and data_brain_mat
x = data_beh_mat;y = data_brain_mat;
means_x=mean(data_beh_mat);stds_x=std(data_beh_mat);
means_y=mean(data_brain_mat);stds_y=std(data_brain_mat);
x=(x-repmat(means_x,size(x,1),1))./repmat(stds_x,size(x,1),1);
y=(y-repmat(means_y,size(y,1),1))./repmat(stds_y,size(y,1),1);
%% this is equivalent of zscoring datalng_beh_mat and datalng_brain_mat based on means of data_beh_mat and data_brain_mat
xlng = datalng_beh_mat;ylng = datalng_brain_mat;
xlng=(xlng-repmat(means_x,size(xlng,1),1))./repmat(stds_x,size(xlng,1),1);
ylng=(ylng-repmat(means_y,size(ylng,1),1))./repmat(stds_y,size(ylng,1),1);

%% tables for CVprediction

%get subjects included in lng projection
[~, idx1] = ismember(tablelng.LONI_ID, tabletrain.LONI_ID);idx1=idx1(idx1>0);
[~, idx2] = ismember(tabletrain.LONI_ID, tablelng.LONI_ID);idx2=idx2(idx2>0);

X_test = [xlng*result_reverse.v(:,1:4),ylng*result_reverse.u(:,1:4)] ;
X_all = [result.usc(:,1:4),result.vsc(:,1:4)];
Xlng = X_test(idx2, :);
Xbase = X_all(idx1, :);

cs = table(tablelng.LONI_ID(idx2, :), dxlng_included(idx2, :), Xbase(:,1),Xbase(:,2),Xbase(:,3),Xbase(:,4),Xbase(:,5),Xbase(:,6),Xbase(:,7),Xbase(:,8),'VariableNames',{'LONI_ID','DX','cog1','cog2','cog3','cog4','bra1','bra2','bra3','bra4'});
lng = table(tablelng.LONI_ID(idx2, :), dxlng_included(idx2, :), X_test(idx2,1),X_test(idx2,2),X_test(idx2,3),X_test(idx2,4),X_test(idx2,5),X_test(idx2,6),X_test(idx2,7),X_test(idx2,8),'VariableNames',{'LONI_ID','DX','cog1','cog2','cog3','cog4','bra1','bra2','bra3','bra4'});

if strcmp(version, 'max')
    save results_PLS_max_lng
elseif strcmp(version, 'min')
    save results_PLS_min_lng
elseif strcmp(version, 'CDR')
    save results_PLS_CDR_lng
elseif strcmp(version, 'BNT')
    save results_PLS_BNT_lng
else
end