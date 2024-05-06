function Local_Interregional_Dim(cur_rec,cur_motif,cur_area,data, num_areas,partitionmethod, computemethod, tempval, normtype)

%Add paths
if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\GenUtils'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Ephys'));
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'));
end

%load data
if nargin <4; data = LoadSubspaceData('in_grouped'); end
if nargin <5; num_areas = 7; end %use all of the n-1 areas as source areas
if nargin <6; partitionmethod = 'rand'; end %can be rand or split
if nargin <7; computemethod = 'rrr'; end %can be rrr or SVCA
if nargin <8; tempval = 0.8; end
if nargin <9; normtype = 'mean'; end


%%

all_areas = data{1}(1).area_label;

neu_range = 15:2:120; %the range of subselected neurons
n_subsample = 75; 
n = numel(neu_range);

x_dim = NaN(1,n,n_subsample);
x_dim_local = NaN(1,n,n_subsample);
rPEV= NaN(1,n,n_subsample); %the amount of reliable explained variance
rPEV_local = NaN(1,n,n_subsample); %the amount of reliable explained variance
curve_dim = NaN(n,100,n_subsample); %take more dims then expectedn
curve_dim_local = NaN(n,100,n_subsample);
num_neu = NaN(1,n,n_subsample);

idx = strcmp(data{cur_rec}(cur_motif).area_label,all_areas{cur_area});  

%get your local source population (randomly selected)
fr_loc = cat(1,data{cur_rec}(cur_motif).area_val{idx==1}); 
fr_loc = normalizeToBaseline(fr_loc,1:2,normtype);
fr_loc = fr_loc(:,3:end,:);                        
nn = floor(size(fr_loc,1)/2); %get the size of the local pop

%grab your interregional source population                                                
fr = cat(1,data{cur_rec}(cur_motif).area_val{idx==0}); 
fr = normalizeToBaseline(fr,1:2,normtype); 
fr = fr(:,3:end,:);
region_idx = data{cur_rec}(cur_motif).area_val(idx==0);
region_idx = arrayfun(@(n) n*ones(size(region_idx{n},1),1), 1:numel(region_idx),'UniformOutput',0);
region_idx = cat(1,region_idx{:});

%subsample the original population
if max(region_idx)>=num_areas %skip because you don't have requisite areas in dataset
for cur_subsample = 1:n_subsample
    for cur_neu_subsample = 1:n        
        cur_n = neu_range(cur_neu_subsample); 
        if cur_n<=nn %can't subsample more than you have
            rng(cur_subsample*cur_n); 
            if strcmp(partitionmethod,'rand') %randomly split local population into two pops 
                neuidx = [zeros(1,nn),ones(1,nn)];
                neuidx((cur_n+1):(numel(neuidx)-cur_n))=NaN; %subselect 
                neuidx = neuidx(randperm(numel(neuidx)));         
            elseif strcmp(partitionmethod,'split') %split local population into two groups based on position along probe        
                neuidx = [zeros(1,nn),ones(1,nn)];   
                %only grab the 'n_sub' #number from distal ends of the probe
                neuidx((cur_n+1):(numel(neuidx)-cur_n))=NaN;        
            else
                error('unknown partition method')
            end

            src = fr_loc(neuidx==1,:,:);
            trg = fr_loc(neuidx==0,:,:);                            

            %split trials into random halves
            tidx = randperm(size(fr_loc,3),floor(size(fr_loc,3)/2));    
            t = ones(size(fr_loc,2),size(fr_loc,3));
            t(:,tidx)=0;
            t = t(:);

            %if applicable, subselect the number of interregional sources
            %do this within the subsample loop to randomly select
            fr_int = fr(ismember(region_idx,randperm(max(region_idx),num_areas)),:,:);                      
        
            %get your distribution of the combined source populations
            full_dist = max(nanmean(cat(1,fr_int,src),[2,3]));
            edges = linspace(0,full_dist,5); 
            edges(end) = ceil(edges(end)); %camden you are doing this because the max value can get cut off by discretize... machine percision. 
            [srcdist,~] = discretize(nanmean(src,[2,3]),edges);  
        
            %discretize the interregional source population
            [intdisp] = discretize(nanmean(fr_int,[2,3]),edges);
        
            %build your interregional source population by
            %grabbing random set that matches the distribution
            %of srcdist
            intregional = NaN(size(src));
            for cur_bin = 1:4                                                                      
                tempidx = find(intdisp == cur_bin);
                srcidx = find(srcdist == cur_bin);
                if numel(srcidx)>numel(tempidx) %drop neurons in local source pop outside range of interreional pop (rare)
                    src(srcdist==cur_bin,:,:)=[];
                    intregional(srcdist==cur_bin,:,:)=[];
                    srcdist(srcdist==cur_bin)=[];            
                elseif ~isempty(srcidx)                 
                    tempidx = tempidx(randperm(numel(tempidx),sum(srcdist==cur_bin))); 
                    intregional(srcdist==cur_bin,:,:) = fr_int(tempidx,:,:);
                else
                end
            end                    
                
            if strcmp(computemethod,'rrr')
                [ndim,d,PEV] = LocalDimRRR(intregional,trg,t,tempval); 
                [ndim_local,d_loc,PEV_local] = LocalDimRRR(src,trg,t,tempval);              
            elseif strcmp(computemethod,'svca')
                [ndim,d,PEV] = SVCA_V2(src,intregional,t,tempval); 
                [ndim_local,d_loc,PEV_local] = SVCA_V2(src,trg,t,tempval);                             
            else
                error('unknown compute method')
            end
              
            %because we are so sparsely sampling the interregional population,
            %sometimes we get a bad draw with no reliable variance
            if isempty(ndim); ndim = 0; end
            if isempty(ndim_local); ndim_local=0; end
            
            rPEV(1,cur_neu_subsample,cur_subsample) = PEV;  
            rPEV_local(1,cur_neu_subsample,cur_subsample) = PEV_local; 
            x_dim(1,cur_neu_subsample,cur_subsample) = getLastD(d); %ndim;  
            x_dim_local(1,cur_neu_subsample,cur_subsample) = getLastD(d_loc); %ndim_local; 
            curve_dim(cur_neu_subsample,1:numel(d),cur_subsample) = d; 
            curve_dim_local(cur_neu_subsample,1:numel(d_loc),cur_subsample) = d_loc; 
            num_neu(1,cur_neu_subsample,cur_subsample) = size(trg,1);
        end        
    end    
    fprintf('\n %d %d %d and subsample %d of %d...',cur_rec, cur_motif, cur_area,cur_subsample,n_subsample);
end
end



%average across subsamples
% rPEV = squeeze(nanmean(rPEV,2));
% rPEV_local = squeeze(nanmean(rPEV_local,2));
% x_dim = squeeze(nanmean(x_dim,3));
% x_dim_local = squeeze(nanmean(x_dim_local,3));
% curve_dim = squeeze(nanmean(curve_dim,2));
% curve_dim_local = squeeze(nanmean(curve_dim_local,2));
% num_neu(1,cur_neu_subsample,cur_subsample) = size(trg,1);
%get number of reliable dimensions
% [d,dauc] = getLastD(d);
% [d_loc,dLocauc] = getLastD(d_loc);    



%%
%save off the data
if ispc
   savedir = 'Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\Local_Interregional_Dimensionality';   
%    savedir = [sprintf('%s_%s_%s_%dareas',savedir,partitionmethod,computemethod,num_areas),'\'];   
   savedir = [sprintf('%s_%s_%s_%dareas_keeppsth',savedir,partitionmethod,computemethod,num_areas),'\'];   
else
   savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/Local_Interregional_Dimensionality';   
%    savedir = [sprintf('%s_%s_%s_%dareas',savedir,partitionmethod,computemethod,num_areas),'/'];      
   savedir = [sprintf('%s_%s_%s_%dareas_keeppsth',savedir,partitionmethod,computemethod,num_areas),'/'];      
end

if ~exist(savedir,'dir')
    mkdir(savedir);
end

save([savedir,sprintf('rec%d_motif%d_area%d_thresh%g_norm%s.mat',cur_rec,cur_motif,cur_area,tempval,normtype)],'partitionmethod','computemethod',...
    'cur_area','cur_rec','cur_motif','tempval','normtype','n_subsample','x_dim','x_dim_local','curve_dim','curve_dim_local','num_neu','rPEV','rPEV_local','neu_range','num_areas');
fprintf('\nsaved %s',savedir)

%%

end %function end


%                 %get distribution of pop
%                 [bins,edges] = discretize(nanmean(fr_loc,[2,3]),4);
%                 %split each bin into two groups
%                 for i = 1:4
%                     tempn=floor(sum(bins==i)/2);
%                     if tempn>0
%                         neuidx = [zeros(1,tempn),ones(1,tempn)];
%                         binidx = find(bins==i,numel(neuidx)); %grab first x neurons in the list (e.g. will drop one neuron at the end if odd number)
% 
%                 end
                %just randomly split




