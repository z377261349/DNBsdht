function [DG] = Selecting(data_s,bin,percentiage,alpha)
rn = size(data_s,1);
cn = size(data_s,2);
nk = floor(rn/bin);
vr = zeros(nk,cn);

for i=1:nk
    vr(i,:) = var(data_s((1+(i-1)*bin):(i*bin),:));  %%Calculate variance
end

mn = median(vr);
md = 1/icdf('norm',percentiage,0,1)*mad(vr,1);
vr_c = bsxfun(@rdivide,bsxfun(@minus,vr,mn),md);     %%Centralize and Standardize

ct = zeros(1,nk);
for i=1:cn
    [pks,locs] = findpeaks(vr_c(:,i));               %%peak detection
    m = size(pks,1);
    pv = 1-cdf('Normal',pks,0,1).^3;                 %%_pvalue
    induce = locs(pv<(alpha/m));                         
    ct(1,induce) = ct(1,induce) + 1;
end
[~,I] = max(ct);                                     %%find the bin with the most peak

DG = false(1,cn);
for i=1:cn
    [~,locs] = findpeaks(vr_c(:,i));
    DG(1,i) = logical(ismember(I,locs));             %%peak genes are DNB genes non-peak genes are non-DNB genes
end
end