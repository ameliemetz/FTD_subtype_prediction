%% PLOTS

%this code produces plots for PLS analysis results (brain patterns)
%other plots produced in R

addpath(genpath('/data/dadmah/resources/Pls'));
addpath(genpath('/data/datamah/in_vivo_Datasets/NIFD'));
addpath(genpath('/data/datamah/in_vivo_Datasets/NIFD/PLS'));

%% Specify the model version: 'max', 'min'

    version = 'max';

if strcmp(version, 'min')
    load results_PLS_min
    c3=1:3;
else
    load results_PLS_max
    c3=1:4;
end

%% plotting covariance explained and p values for the latent variables

figure
colororder({'#D55E00'});
plot((result.s.^2)./sum(result.s.^2),'MarkerSize',14, 'Color' ,'#0072B2', 'Marker', '.', 'LineStyle','none');hold on;xlabel('Latent Variable');ylabel('% Covariance');hold on
plot(result.perm_result.sprob,'*','MarkerSize',6, 'Color','#D55E00');hline = refline([0 0.05]);hline.Color = '#D55E00';hline.LineStyle = "--";
title('Covariance explained and permutation p-values');legend({'covariance','p value'})
yyaxis right
ylabel('p values', 'Color','#D55E00');

if strcmp(version, 'max')
    print(strcat('plot_signLVs_max',cell2mat(dx_grps(dx_grp))),'-dpng');
else
    print(strcat('plot_signLVs_min',cell2mat(dx_grps(dx_grp))),'-dpng');
end

%% thresholded brain maps (we are using this code: https://github.com/VANDAlab/Plotting-Tools/blob/main/Plotting_Atlas.m)

for comp=c3
    cerebra_regions=double(result.boot_result.orig_corr(:,comp));
    pass_region = (sign(result.boot_result.llcorr(:,comp).*result.boot_result.ulcorr(:,comp)))>0;
    cerebra_regions = cerebra_regions.*pass_region;
    figure;
    plot_number=1;
    for j=[50,80,90,100]
        Plotting_Atlas(ones(102,1),cerebra_regions,'axial',j,'',2,2,plot_number);plot_number=plot_number+1;
    end
end
    if strcmp(version, 'max')
        switch comp
            case 1
                print(strcat('plot_brainLV1_max',cell2mat(dx_grps(dx_grp))),'-dpng');
            case 2
                print(strcat('plot_brainLV2_max',cell2mat(dx_grps(dx_grp))),'-dpng');
            case 3
                print(strcat('plot_brainLV3_max',cell2mat(dx_grps(dx_grp))),'-dpng');
            case 4
                print(strcat('plot_brainLV4_max',cell2mat(dx_grps(dx_grp))),'-dpng');
        end
    else
        switch comp
            case 1
                print(strcat('plot_brainLV1_min',cell2mat(dx_grps(dx_grp))),'-dpng');
            case 2
                print(strcat('plot_brainLV2_min',cell2mat(dx_grps(dx_grp))),'-dpng');
            case 3
                print(strcat('plot_brainLV3_min',cell2mat(dx_grps(dx_grp))),'-dpng');
            end
    end

end

cmap = colormap("jet"); %get current colormap
cmap=cmap([0 0.3],:); % set your range here
colormap(cmap); % apply new colormap
colorbar()

%% thresholded surface brain maps

brainregions=readtable('cerebraindices - Sheet2.csv','Delimiter',',');

for comp=c3

    cerebra_regions=double(result.boot_result.orig_corr(:,comp));
    pass_region = (sign(result.boot_result.llcorr(:,comp).*result.boot_result.ulcorr(:,comp)))>0;
    cerebra_regions = cerebra_regions.*pass_region;

    tools = '/data/zeiyas/tools/';
    addpath(genpath(strcat(tools,'matlab_toolboxes/surfstat')));
    addpath(genpath(strcat(tools,'useful_scripts')));
    % read surface template
    model_folder_s = strcat(tools,'parcellations_atlas_mni/models/');
    s = SurfStatReadSurf( {strcat(model_folder_s,'surface_models/icbm_avg_mid_sym_mc_left.obj'),...
        strcat(model_folder_s,'surface_models/icbm_avg_mid_sym_mc_right.obj')} );
    % make Cerebra surface
    atlases = load(strcat(tools,'parcellations_atlas_mni/yz_atlasses/all_atlases_vectorized.mat'));
    Cerebra_s = ICBM_volume_to_surface_map(atlases.all_atlases_volume.Cerebra_gm_h);

    Cerebra_s_stat = tabulate(Cerebra_s);
    remove_threshold = 10;remove_regions = Cerebra_s_stat((Cerebra_s_stat(:,2)<remove_threshold),1);
    Cerebra_s(ismember(Cerebra_s,remove_regions))=0;

    complementary_atlas= atlases.all_atlases_surface.Schaefer_7_1000;
    Cerebra_s(complementary_atlas==0)=0;
    ind_miss = find(Cerebra_s==0 & complementary_atlas~=0);

    for i=1:size(ind_miss)
        temp_cerebra = Cerebra_s(complementary_atlas==complementary_atlas(ind_miss(i)));
        temp_val(i,1)= mode(temp_cerebra(temp_cerebra~=0));
    end
    Cerebra_s(ind_miss)=temp_val;clear temp_val
    Cerebra_s(isnan(Cerebra_s))=0;

    u_regions = unique(Cerebra_s);u_regions(u_regions==0)=[];
    out = cerebra_regions;out = out(u_regions);
    out_s =region_to_atlas(out,Cerebra_s);

    figure(comp);SurfStatViewData_yz_22(out_s,s,[-.6,.6],strcat('bootstrap ratio'));

    if strcmp(version, 'max')
        switch comp
            case 1
                print(strcat('plot_surfaceLV1_max',cell2mat(dx_grps(dx_grp))),'-dpng');
            case 2
                print(strcat('plot_surfaceLV2_max',cell2mat(dx_grps(dx_grp))),'-dpng');
            case 3
                print(strcat('plot_surfaceLV3_max',cell2mat(dx_grps(dx_grp))),'-dpng');
            case 4
                print(strcat('plot_surfaceLV4_max',cell2mat(dx_grps(dx_grp))),'-dpng');
        end
    else
        switch comp
            case 1
                print(strcat('plot_surfaceLV1_min',cell2mat(dx_grps(dx_grp))),'-dpng');
            case 2
                print(strcat('plot_surfaceLV2_min',cell2mat(dx_grps(dx_grp))),'-dpng');
            case 3
                print(strcat('plot_surfaceLV3_min',cell2mat(dx_grps(dx_grp))),'-dpng');
           end
    end

    %print 20 brain regions with highest bootstrap ratios
    [a,b]=sort(cerebra_regions, "ascend");
    b(1:10,:)
    a(1:10,:)
    brainregions(ismember(brainregions.Var2, b(1:10,:)),:)

    [a,b]=sort(cerebra_regions, "descend");
    b(1:10,:)
    a(1:10,:)
    brainregions(ismember(brainregions.Var2, b(1:10,:)),:)

end

%% plot relationship between brain and cognition scores

if strcmp(version,'max')
    for comp = c3
        figure;
        scatter(result.vsc(:,comp), result.usc(:,comp), [], data_beh_mat(:,4), 'filled'),lsline()
        xlabel('Brain Score'); ylabel('Cognitive Score');
        colormap("copper")
        colorbar;
        c = colorbar; % Creating colorbar object
        c.Label.String = 'MMSE Score'; % Adding label to colorbar

        switch comp
            case 1
                title('Latent variable I');
                print(strcat('plot_surfaceLV1_max',cell2mat(dx_grps(dx_grp))),'-dpng');
            case 2
                title('Latent variable II');
                print(strcat('plot_surfaceLV2_max',cell2mat(dx_grps(dx_grp))),'-dpng');
            case 3
                title('Latent variable III');
                print(strcat('plot_surfaceLV3_max',cell2mat(dx_grps(dx_grp))),'-dpng');
            case 4
                title('Latent variable IV');
                print(strcat('plot_surfaceLV4_max',cell2mat(dx_grps(dx_grp))),'-dpng');
        end
    end
else
end
