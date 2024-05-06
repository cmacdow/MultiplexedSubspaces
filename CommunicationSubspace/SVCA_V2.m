function [ndim,d,rPEV] = SVCA_V2(x,y,t,cutoff,npc)
%Camden MacDowell - timeless
%wrapper for stringer et., science 2019
if nargin <4; cutoff = []; end
if nargin <5; npc = floor(size(x,1)/2); end



npc = min(npc,floor(size(x,1)/2));

%subtract the psth
x = x-nanmean(x,3);
y = y-nanmean(y,3);

%concatentate across trials 
x = reshape(x,[size(x,1),size(x,2)*size(x,3)]);
y = reshape(y,[size(y,1),size(y,2)*size(y,3)]);

%get the timepoints
itrain = find(t==1);
itest = find(t==0);

xx = cat(1,x,y);
xx = xx-nanmean(xx,2); %mean center
ntrain = 1:size(x,1);
ntest = (size(x,1)+1):size(xx,1);

rPEV=0;
%% run SVCA
[sneur, varneur, ~, ~, ~, ~]=runSVCA(xx, npc, ntrain, ntest, itrain, itest);
% [sneur, varneur, ~, ~, s1, s2]=runSVCA(xx, npc, ntrain, ntrain, itrain, itrain); %this should be 100%
%%

if ~isempty(cutoff) %set a threshold in reliable variance
    d = cumsum(sneur./varneur);
    d = d/max(d); 
    ndim = find(d>=cutoff,1,'first');
else %get AUC
    d = sneur./varneur;
    %stop at first negative
    d(find(d<0,1,'first'):end)=[];
    d = cumsum(d)/max(cumsum(d));
    d = [0, d', 1];
    ndim = 1-trapz(d/numel(d));
end

d = cumsum(sneur./varneur);

end %function end


function [sneur, varneur, u, v, s1, s2] = runSVCA(Ff, npc, ntrain, ntest, itrain, itest)
%from stringer et., science 2019
cov = Ff(ntrain,itrain) * Ff(ntest,itrain)';
[u,~,v] = svdecon(cov);
u = u(:,1:npc);
v = v(:,1:npc);
s1 = u' * Ff(ntrain,itest);
s2 = v' * Ff(ntest,itest);
sneur = sum(s1 .* s2, 2);
varneur = sum(s1.^2 + s2.^2,2)/2;

end

% INPUTS: 
%     Ff (neurons x timepts)
%     npc (number of PCs to compute)
%     ntrain (one half of neurons)
%     ntest (other half of neurons)
%     itrain (one half of timepts)
%     itest (other half of timepts)
% OUTPUTS:
%     sneur (shared variance of each covariance component)
%     vneur (total variance of each covariance component)
%     u (left eigenvectors of covariance matrix btw ntrain and ntest on
%        itrain timepts)
%     v (right eigenvectors of covariance matrix btw ntrain and ntest on
%        itrain timepts)