%Figure 1 
%Camden MacDowell

%% Probe tracing and visualization
fp = fig_params_cortdynamics;
%uncomment to run (these can take a bit and are older functions)
% AnatomicalCorrelellogram(EphysPath,probe_ccf_path,allen_atlas_path,mua_flag,fp,cur_probe)
% PlotCombinedProbeTrajectories(histo_wrkspace_path,spike_opts_path)

% Get a summary of the in/out correlations across recordings per area
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\ValidationProbeLocation';
Plot_AnatomicalCorrelationInOut(1,'general');
fp.FigureSizing(gcf,[3 2 4 4],[10 10 10 10])
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'In_OutCorrelations',savedir,0); close all



savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\ValidationProbeLocation';
for i = 1:6
    PlotExampleCorrellelogram(i)
    fp.FigureSizing(gcf,[3 2 4 4],[10 10 10 10])
    set(gca,'CLim',[0 0.025])    
end
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'ActivityCorrelelogram',savedir,0); close all


%% Plot example raw data 
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\ExampleData';
if ~exist(savedir,'dir'); mkdir(savedir); end

cur_rec = 1;
tp = [1200 2100];
PlotExampleEphys(1,1,tp)
fp.FigureSizing(gcf,[3 2 16 6],[2 10 30 15])
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'ExampleEphys',savedir,0); close all

PlotExampleImagingTrace(cur_rec,tp)
fp.FigureSizing(gcf,[3 2 16 3],[2 10 40 15]); box on
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'ExampleImaging',savedir,0); close all

%plot motif onsets
motif_num = [3,5,14];
PlotExampleMotifWeightings(cur_rec,motif_num,tp); 
fp.FigureSizing(gcf,[3 2 16 1],[2 10 40 15])
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'ExampleMotifOnsets',savedir,0); close all

%plot motifs
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\Motifs';
if ~exist(savedir,'dir'); mkdir(savedir); end
PlotBasisMotifs(savedir) 

%plot the average response of ephys areas to motifs
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\EphysResponse';
if ~exist(savedir,'dir'); mkdir(savedir); end
Plot_MotifEphysActivity(1,1:14,1,'mean',savedir)


%% plot the motifs clustering
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\MotifClustering';
if ~exist(savedir,'dir'); mkdir(savedir); end
load('Mousepermouse_basis_motifs331.mat','tcorr_mat','cluster_idx');
idx = find(cluster_idx==2);
tcorr_mat(idx,:)=[];
tcorr_mat(:,idx)=[];
cluster_idx(idx)=[];
cluster_idx(cluster_idx>1)=cluster_idx(cluster_idx>1)-1;
Plot_OrderedSimilarityMatrix(tcorr_mat,cluster_idx);
set(gca,'clim',[0.2 0.6])
load('Mousepermouse_basis_motifs332.mat','tcorr_mat','cluster_idx');
idx = find(cluster_idx==2);
tcorr_mat(idx,:)=[];
tcorr_mat(:,idx)=[];
cluster_idx(idx)=[];
cluster_idx(cluster_idx>1)=cluster_idx(cluster_idx>1)-1;
Plot_OrderedSimilarityMatrix(tcorr_mat,cluster_idx);
set(gca,'clim',[0.2 0.6])
Plot_OrderedSimilarityMatrix(tcorr_mat,cluster_idx);
load('Mousepermouse_basis_motifs334.mat','tcorr_mat','cluster_idx');
idx = find(cluster_idx==2);
tcorr_mat(idx,:)=[];
tcorr_mat(:,idx)=[];
cluster_idx(idx)=[];
cluster_idx(cluster_idx>1)=cluster_idx(cluster_idx>1)-1;
Plot_OrderedSimilarityMatrix(tcorr_mat,cluster_idx);
set(gca,'clim',[0.2 0.6])
Plot_OrderedSimilarityMatrix(tcorr_mat,cluster_idx);
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'ExampleMotifOnsets',savedir,0); close all


%% Load subspace data
% data = LoadSubspaceData('in_grouped');
data = LoadSubspaceData('in_grouped_full_mean');
dataout = LoadSubspaceData('out_grouped');

%% ridge regression figures
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\RidgeRegressionFigures';
if ~exist(savedir,'dir'); mkdir(savedir); end
%for IN models
RidgeRegressionExampleFigures(data,'VIS',[5,3])
RidgeRegressionCombinedFigures(data)
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'RegressionExample_IN',savedir,0); close all


%% get the average 
[r,~] = LoadVariable(data,'rel_performance',[]);
%5 dimensions in visual 
x = squeeze(r(:,:,8,5));
p = nanmean(x(:));
pci = bootci(1000,@nanmean,x(:));

%10 dimension in visual 
x = squeeze(r(:,:,8,10));
p = nanmean(x(:));
pci = bootci(1000,@nanmean,x(:));

%5 dimensions across everything
x = squeeze(r(:,:,:,5));
p = nanmean(x(:));
pci = bootci(1000,@nanmean,x(:));

%10 dimensions across everything
x = squeeze(r(:,:,:,10));
p = nanmean(x(:));
pci = bootci(1000,@nanmean,x(:));


%load the shuffled models | THIS WILL TAKE DAYS, reccomend to do on spock. 
% [fn,~] = GrabFiles('\w*shuf\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\RidgePermutations'});
% temp = [];
% for i = 1:numel(fn)
%     try
%         temp{i} = load(fn{i},'rsq');
%     catch
%     end
%     
% end
% x = cellfun(@(x) x.rsq,temp,'UniformOutput',0);
% x = cat(2,x{:});
% max(x(:))

%for OUT models
RidgeRegressionExampleFigures(dataout,'VIS',[5,3])
RidgeRegressionCombinedFigures(dataout)
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'RegressionExample_OUT',savedir,0); close all

%% Get the explained variance of each motif
[fn,~] = GrabFiles('\w*.mat',0,{'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\BasisMotifFits'}); 
data = cellfun(@(x) load(x,'stats_refit'),fn);
x = arrayfun(@(n) data(n).stats_refit.pev,1:numel(data));
nanmean(x)
bootci(1000,@nanmean,x)

%next to add: 
%plot showing the correlation between data in and data out... explain in
%text that this is not required and highlight it's difference in strength
%between different areas. 
%plot showing the decodability in average |trial-to-trial activity| across all neurons between each pair of motifs.  
%the permutation tests
%by end of next week | knock out figure 2 and have explored 3


%% Beta weight fits
folder = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\CrossValidationSimilarity\';
[fn,~] = GrabFiles('\w*.mat',0,{folder}); 





%%




















