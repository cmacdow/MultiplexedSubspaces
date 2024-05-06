function CompareEphyCalciumNetworks(data,data_all)
%Pull together AUC of the % of neurons vs |B|
%data is ephys | data = LoadSubspaceData('in_grouped');
%data_all is widefield | data_all = LoadCorticalNetworks(0);

[~, area_all] = LoadVariable(data,'rrr_beta','VIS',1); %load area list

%cortical area list | left hemi (right side img)
%MOs, RSP, SS, SSBFD, VIS to match the order above
% coords = [13,30; 36,30; 36,26; 27,17; 46,18]; %right hemi (left img) 
coords = [13,42; 36,42; 36,46; 25,53; 46,52];
cort_areas = [2,4,5,6,8]; %remove areas 1, 3, and 7
cort_names = [{'MOs'},{'RSP'},{'SS'},{'SSp-bfd'},{'VIS'}];

%remake so that it does it per dataset 

%% get all of the AUC orders
rho_data = NaN(10,numel(area_all),84);
% perm_idx = perms(1:4);
% rho_data_perm = NaN(10,numel(area_all),84,size(perm_idx,1)); %4! = 24 different perms


for cur_d = 1:10 %use the top 10 dimensions
    fprintf('\n working on dim %d... \n',cur_d);
    for cur_area = cort_areas %loop through the areas
        
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

        % loop through each and find corresponding value in img data
        rho_all = NaN(84,numel(cort_areas)-1);
        auc_all = NaN(84,numel(cort_areas)-1);
        for i = 1:size(auc,1)            
            cur_rec = datasetidx(i,1); %rec
            cur_m = datasetidx(i,2); %motif
                         
            %get cortical network
            cort = data_all{cur_rec};
            cort_idx = find([cort.cur_motif]==cur_m & [cort.cur_a] == cur_area);
            if ~isempty(cort_idx) && sum(~isnan(auc(i,:)))>0
                cort = cort(cort_idx).rho_all(:,:,cur_d);
                % get the cortical correlations
                rho = diag(cort(coords(:,1),coords(:,2))); 
                %remove the areas not recorded on the cortex from AUCs
                temp_areas = AreaLabel{cur_rec,cur_m};
                temp_auc = auc(i,:); 
                temp_auc(isnan(temp_auc))=[];
                temp_auc = temp_auc(ismember(temp_areas,area_all(cort_areas))); 
                %remove self and missing areas from rho
                rho(ismember(cort_names,area_all(~ismember(area_all,temp_areas))))=[];
                %store in big matrixsince not enough areas per dataset for corr
                rho_all(i,1:numel(rho))=rho;
                auc_all(i,1:numel(temp_auc))=temp_auc;
                %get the correlation
                rho_data(cur_d,cur_area,i) = corr(rho(:),temp_auc(:),'type','Spearman');
               
            end
        end
    end
    fprintf('\n...done with %d of %d',cur_d,10)
end
% rho_data_perm = reshape(rho_data_perm,10,numel(area_all),84*size(perm_idx,1));
%% Plot the correlation across dimensions
xboot = NaN(1000,10);
for i = 1:10 
    x = squeeze(rho_data(i,:,:)); 
    xboot(:,i) = pairedBootstrap(x(:),@nanmean);    
end
pval = 1-nansum(xboot>0)/size(xboot,1);
fp = fig_params_cortdynamics;
figure; hold on; 
CompareViolins(xboot',fp,'plotspread',0,'connectline',[0.5 0.5 0.5],'plotspread',0,'divfactor',.5,'sidebyside',0,'distWidth',0.5);
for i = 1:numel(pval)
    if pval(i)<0.05
        plot([i i],[-.15 -.15],'*','markersize',5,'color','k');
    end
end
plot([0 10],[0,0],'linewidth',1,'LineStyle','--','Color','k')
fp.FormatAxes(gca);
xlabel('Dimension');
ylabel(' Rho');
box on; grid off
fp.FigureSizing(gcf,[3 2 3.75 3.75],[10 10 10 10]) 

%get some Stats

%% Plot the correlation across everything
xboot = pairedBootstrap(rho_data(:),@nanmean);
figure; hold on;
histogram(xboot,'EdgeColor','none','FaceColor','k','FaceAlpha',0.3,'EdgeAlpha',0.3,'NumBins',15); 
plot([0 0],get(gca,'ylim'),'color','r','LineStyle','--','LineWidth',1);
xlim([-0.1 0.1])
fp.FormatAxes(gca);
xlabel('Rho');
title({'Ephys-to-Cortical','Network Similarity'},'Fontweight','normal')
ylabel('Samples')
box on; grid off
fp.FigureSizing(gcf,[3 2 2.25 2.25],[10 10 10 10]) 

%% Plot the correlation across everything then split into first 5 and second 5
xboot = pairedBootstrap(rho_data(:),@nanmean);
x = squeeze(rho_data(1:5,:,:)); 
xboot = [xboot,pairedBootstrap(x(:),@nanmean)]; 
x = squeeze(rho_data(6:10,:,:)); 
xboot = [xboot,pairedBootstrap(x(:),@nanmean)]; 
pval = 1-nansum(xboot>0)/size(xboot,1);
fp = fig_params_cortdynamics;
figure; hold on; 
CompareViolins(xboot',fp,'plotspread',0,'connectline',[],'plotspread',0,'divfactor',.5,'sidebyside',0,'distWidth',0.95);
for i = 1:numel(pval)
    if pval(i)<0.05
        plot([i i],[-.05 -.05],'*','markersize',5,'color','k');
    end
end
plot([-1 4.5],[0,0],'linewidth',1,'LineStyle','--','Color','k')
fp.FormatAxes(gca);
title({'Ephys-to-Cortical','Network Similarity'},'Fontweight','normal')
ylabel('Rho');
xlim([0 4])
set(gca,'XTickLabel',{'All','1-5','6-10'},'XTickLabelRotation',90)
xlabel(sprintf('Dimensions p=%0.2g',pval(3)))
box on; grid off
fp.FigureSizing(gcf,[3 2 2.5 3.75],[10 10 10 10]) 

%% Plot the correlation across the first dimension then all dimensions
[xboot,stats] = pairedBootstrap(rho_data(:),@nanmean);
x = squeeze(rho_data(1,:,:)); 
xboot = [xboot,pairedBootstrap(x(:),@nanmean)]; 
x = squeeze(rho_data(2,:,:)); 
xboot = [xboot,pairedBootstrap(x(:),@nanmean)]; 
pval = 1-nansum(xboot>0)/size(xboot,1);
fp = fig_params_cortdynamics;
figure; hold on; 
CompareViolins(xboot(:,2:end)',fp,'plotspread',0,'connectline',[],'plotspread',0,'divfactor',.5,'sidebyside',0,'distWidth',0.95);
for i = 2:numel(pval)
    if pval(i)<0.05
        plot([i-1 i-1],[-.05 -.05],'*','markersize',5,'color','k');
    end
end
plot([-1 3.5],[0,0],'linewidth',1,'LineStyle','--','Color','k')
fp.FormatAxes(gca);
title({'Ephys-to-Cortical','Network Similarity'},'Fontweight','normal')
ylabel('Rho');
xlim([0 3])
set(gca,'XTickLabel',{'1','2'},'XTickLabelRotation',90)
xlabel(sprintf('Dimensions'))
box on; grid off
fp.FigureSizing(gcf,[3 2 2.5 3.75],[10 10 10 10]) 

%%
savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Figures\EphysCorticalNetworksSimilarity';
saveCurFigs(get(groot, 'Children'),{'-dpng','-dsvg'},'Similarity',savedir,0); close all


% %% plot an example to confirm
% close; coords = [13,42; 36,42; 36,46; 25,53; 46,52];
% x = data_all{1,1}(57).rho_all;
% imagesc(x(:,:,10));
% hold on
% for i = 1:5
%     plot(coords(i,2),coords(i,1),'marker','x','color','r','MarkerSize',10); %flipped when plotting
%     rho = x(coords(i,1),coords(i,2),10); %can manually confirm with fig
% end


end %function end

function [rho] = CortRho(x,coords)
    rho = NaN(1,size(coords,1));
    for i = 1:size(coords,1)
       %get +/i 2 pixels in all directions
       r = 2;
       rho(i) = nanmean(x(coords(i,1)-r:coords(i,1)+r,coords(i,2)-r:coords(i,2)+r),'all');
    end
end








