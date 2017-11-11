function [Isc] = Calculate(data_s,bin,DG,type)
nk = floor(size(data_s,1)/bin);
varc = zeros(nk,1);
pcc_in = zeros(nk,1);
if(type==2)
    pcc_out = zeros(nk,1);
end
for i=1:nk
    varc(i,1) = mean(var(data_s((1+bin*(i-1)):(bin*i),DG)));
    middle = tril(abs(corr(data_s((1+bin*(i-1)):(bin*i),:))),-1);
    middle1 = middle(DG,DG);
    middle1 = middle1(middle1~=0);
    if(type==2)
        middle2 = middle(~DG,DG);
        middle2 = middle2(middle2~=0);
    end
    pcc_in(i,1) = nanmean(middle1(:));
    if(type==2)
        pcc_out(i,1) = nanmean(middle2(:));
    end
end
if(type==2)
    Isc = varc.*pcc_in./pcc_out;                     %%Iscore
else
    Isc = varc.*pcc_in;
end
end