function CompareDatasets(data,psth)
%Camden - timeless
%Creates a collection of supplemental figures comparing the results when
%performing RRR on PSTH-subtracted or PSTH data. 


%As state in the R2R

% these findings, more predictive dimensions were needed to predict the full activity than the 
% residual variability. Similarly, The hierarchy of regions was maintained, and similar subspace 
% networks were observed. Analyzing the ‘full’ data without removing the average evoked activity 
% just increased the relative number of both local and interregional dimensions. This is expected 
% since the full activity is composed of both the evoked component and residual variability. 


if nargin <1; data = LoadSubspaceData('in_grouped'); end
if nargin <2; psth = LoadSubspaceData('in_grouped_psth'); end

%% remake the Fig 1 dimensionality plots  

%important. need to go into the function and uncomment for psth
PlotSubspaceDimensionality_v2(psth,'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\KeepPSTH');

%% FIG 1 findings: compare ridge regression performance
%ridge repression performance was similar in both datasets
RidgeRegressionCombinedFigures(data)
title('NoPSTH')
RidgeRegressionCombinedFigures(psth)
title('PSTH')
[x,~] = LoadVariable(data,'ridge_performance',[]);
[y,~] = LoadVariable(psth,'ridge_performance',[]);
simcompare(x,y)
title({'RRR Explains more variance','with PSTH included'},'fontweight','normal')
plot([0 0.5],[0 0.5],'color','r','linewidth',1,'linestyle','--')
xlim([0 0.5])
ylim([0 0.5])
[~,stats_modelperformance] = pairedBootstrap([x(:),y(:)],@nanmean);

%% Get the number of dimensions needed to explain the variance
% more predictive dimensions were needed to predict the full activity than the 
% residual variability. 
% Analyzing the ‘full’ data without removing the average evoked activity just increased the 
% relative number of both local and interregional dimensions. This is expected since the full 
% activity is composed of both the evoked component and residual variability. 

[x,~] = LoadVariable(data,'rrr_dim',[],0.8);
[y,~] = LoadVariable(psth,'rrr_dim',[],0.8);
simcompare(x,y)
title({'More dimensions are needed to','predict full activity than residual variability'},'fontweight','normal')
plot([0 30],[0 30],'color','r','linewidth',1,'linestyle','--')
[~,stats_numdimensions] = pairedBootstrap([x(:),y(:)],@nanmean);
save('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\KeepPSTH\dimensionsStats.mat','stats_numdimensions','stats_modelperformance')
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'similarity_dimensions_psth','Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\KeepPSTH\',0); close all


%% Plot the visual subspace example
PlotExampleSubspace(data,'VIS',1:14,1)
PlotExampleSubspace(psth,'VIS',1:14,1)

%% Fig 2 findings: 
% In addition, similar subspace networks were observed using both preprocessing approaches.
% For example, Fig. SX shows the same AUC plot as in the original main text. This generalized
% to the entire dataset, and show the correlation in structure of the subspace networks throughout. 
fp = fig_params_cortdynamics;
%get the correlation between the aucs
[rho_data,rho_perm] = getCorr(data,psth);
%%
%plot per area and dimension
col = arrayfun(@(n) fp.c_area(n,:),1:size(fp.c_area,1),'UniformOutput',0);
figure; hold on; 
t = tiledlayout(3,3,'TileSpacing','tight');
for i = 1:8
    nexttile
    x = squeeze(rho_data(:,i,:));
    xboot = pairedBootstrap(x',@nanmean);
    CompareViolins(xboot',fp,'plotspread',0,'connectline',col{i},'plotspread',0,'divfactor',.5,'sidebyside',0,'distWidth',0.9,'col',repmat(col(i),1,10));
    fp.FormatAxes(gca);
    ylim([0 1])
    box on; grid off    
end
fp.FigureSizing(gcf,[3 2 8 8],[5 5 15 15]) 
title(t,'Subspace similarity with/without PSTH')
ylabel(t,'Rho')
xlabel(t,'Dimension')

%permuted data
col = arrayfun(@(n) fp.c_area(n,:),1:size(fp.c_area,1),'UniformOutput',0);
figure; hold on; 
t = tiledlayout(3,3,'TileSpacing','tight');
for i = 1:8
    nexttile
    x = squeeze(rho_perm(:,i,:));
    xboot = pairedBootstrap(x',@nanmean);
    CompareViolins(xboot',fp,'plotspread',0,'connectline',col{i},'plotspread',0,'divfactor',.5,'sidebyside',0,'distWidth',0.9,'col',repmat(col(i),1,10));
    fp.FormatAxes(gca);
    ylim([-1 1])
    box on; grid off    
end
fp.FigureSizing(gcf,[3 2 8 8],[5 5 15 15]) 
title(t,'Subspace similarity with/without PSTH PERMUTED DATA')
ylabel(t,'Rho')
xlabel(t,'Dimension')

%% Create a plot just showing the first dimension for each region | then include stats showing all dimensions
col = arrayfun(@(n) fp.c_area(n,:),1:size(fp.c_area,1),'UniformOutput',0);
figure; hold on; 
[xboot,stats_similarity_by_region] = pairedBootstrap(squeeze(rho_data(1,:,:))',@nanmean);

%reorder by decreasing similarity for plotting
[~,idx] = sort(stats_similarity_by_region.mean,'descend'); 
xboot = xboot(:,idx);
col = col(idx);

CompareViolins(xboot',fp,'plotspread',0,'plotspread',0,'divfactor',.5,'sidebyside',0,'distWidth',0.9,'col',col);
fp.FormatAxes(gca);
ylim([0 1])
box on; grid on 

fp.FigureSizing(gcf,[3 2 3.75 3.75],[5 5 15 15]) 
title({'Subspace dimension 1','similarity with/without PSTH'},'FontWeight','normal')
ylabel('Rho')
xlabel('Region')
% mean and CI collectively across all regions and all dimensions
[~,stats_similarity_allregions_alldim] = pairedBootstrap(rho_data(:),@nanmean);
save('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\KeepPSTH\similaritystats.mat','stats_similarity_by_region','stats_similarity_allregions_alldim')

saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'similarity_subspaces','Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\KeepPSTH\',0); close all


end 



function simcompare(x,y,fp)
if nargin <3; fp = fig_params_cortdynamics; end 
figure; hold on; 
s1 = scatter(x(:),y(:),'o','MarkerEdgeAlpha',0.25);  
s1.SizeData = s1.SizeData/5; 
s1.MarkerEdgeColor = 'k';   
xlabel('Without PSTH'); ylabel('With PSTH');
fp.FormatAxes(gca);
box on; grid off
fp.FigureSizing(gcf,[3 2 3.75 3.75],[10 10 15 10]) 

end



function [rho_data,rho_perm] = getCorr(data,psth)
    
    [~, area_all] = LoadVariable(data,'rrr_beta','VIS',1); %load area list
    rho_data = NaN(10,numel(area_all),84);
    rho_perm = NaN(10,numel(area_all),84);
    for cur_d = 1:10 %use the top 10 dimensions
        for cur_area = 1:8 %loop through the areas
            %% do it for the data
            targ_area= area_all{cur_area};
            B = cell(6,14);
            AArea = cell(6,14);
            AreaLabel = cell(6,14);
    
            for cur_rec = 1:6
                for cur_m = 1:14        
                    [B{cur_rec,cur_m},AArea{cur_rec,cur_m},AreaLabel{cur_rec,cur_m}] = loadVisBeta(data,cur_d,cur_rec,cur_m,targ_area);
                end
            end
            
            %sweep fractions
            idx = 0:0.1:1;
            x = NaN(84,7,numel(idx));
            for i = 1:numel(idx)
                [~,y,~,~,datasetidx] = SubspaceUniformity(B,AArea,AreaLabel,data,targ_area,idx(i),1);
                x(1:size(y,1),1:size(y,2),i) = y;
            end
            
            %get the order of the AUCs for each dataset
            x = x/100;
            auc = NaN(size(x,1),size(x,2));        
            for i = 1:size(x,1)
                for j = 1:size(x,2)
                    xx = squeeze(x(i,j,:));
                    auc(i,j) = trapz(xx/numel(xx)); 
                end
            end

            %% do for the psth
            targ_area= area_all{cur_area};
            B = cell(6,14);
            AArea = cell(6,14);
            AreaLabel = cell(6,14);
    
            for cur_rec = 1:6
                for cur_m = 1:14        
                    [B{cur_rec,cur_m},AArea{cur_rec,cur_m},AreaLabel{cur_rec,cur_m}] = loadVisBeta(psth,cur_d,cur_rec,cur_m,targ_area);
                end
            end
            
            %sweep fractions
            idx = 0:0.1:1;
            x = NaN(84,7,numel(idx));
            for i = 1:numel(idx)
                [~,y,~,~,~] = SubspaceUniformity(B,AArea,AreaLabel,psth,targ_area,idx(i),1);
                x(1:size(y,1),1:size(y,2),i) = y;
            end
            
            %get the order of the AUCs for each dataset
            x = x/100;
            aucpsth = NaN(size(x,1),size(x,2));        
            for i = 1:size(x,1)
                for j = 1:size(x,2)
                    xx = squeeze(x(i,j,:));
                    aucpsth(i,j) = trapz(xx/numel(xx)); 
                end
            end

            %% now compare aucpsth and auc
            for i = 1:size(auc,1)            
                if  sum(~isnan(auc(i,:)))>0
                    rho_data(cur_d,cur_area,i) = corr(auc(i,:)',aucpsth(i,:)','type','Pearson');
                    a=auc(i,:)'; 
                    a = a(randperm(numel(a),numel(a)));
                    rho_perm(cur_d,cur_area,i) = corr(a,aucpsth(i,:)','type','Pearson');
                end
            end
        end
        fprintf('\n...done with %d of %d',cur_d,10)
    end

end





















