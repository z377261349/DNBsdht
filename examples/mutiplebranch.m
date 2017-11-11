%%%%%%%% Load data %%%%%%%%
load('tribranch.mat');

%%%%%%%% Parameter %%%%%%%%
rng(1);
example = 3;
bin = 40;
percentiage = 0.75;
alpha = 0.05;

%%%%%%%% Preprocessing %%%%%%%%
trunk = 1200;
bran = 800;
slc1 = false(4*bran,1);
slc1(sort(randperm(4*bran,bran))) = true;
BR1 = zeros(1,4*bran);
BR1(slc1) = (1:bran)+trunk;

slc2 = false(2*bran,1);
slc2(sort(randperm(2*bran,bran))) = true;
BR2 = zeros(1,3*bran);
BR2(1:bran) = (1:bran) + trunk + bran;
BR2([false(bran,1);slc2]) = (1:bran) + trunk + 2*bran;
BR2([false(bran,1);~slc2]) = (1:bran) + trunk + 3*bran;
BR1(~slc1) = BR2;

branch = zeros(4400,1);
branch(1:trunk) = 1;
branch([false(trunk,1);slc1]) = 2;
branch1 = zeros(2400,1);
branch1(1:bran) = 3;
branch1([false(bran,1);slc2]) = 4;
branch1([false(bran,1);~slc2]) = 5;
branch([false(trunk,1);~slc1]) = branch1;

data_s = markerdata1([1:trunk,BR1],:);
%%%%%%%% Draw gene expression %%%%%%%%
figure;
for i=1:16
    subplot(4,4,i);
    scatter(1:4400,data_s(:,i),[],branch);
    hold on
    if(mod(i,4)==1)
    ylabel({'gene';'expression'})
    end
end
suptitle('Gene expression of DPT-ordering cells colored by branches')
saveas(gcf,'Figure_3_1','epsc')
%%%%%%%% Select DNB Group %%%%%%%%
DG = [true(1,8),false(1,8)];
%%%%%%%% Calculating Iscore %%%%%%%%
Isc = Calculate(data_s,bin,DG,2);
%%%%%%%% Hypothesis Testing %%%%%%%%
[Isc_c,p_value,locs] = Hypothesis_Testing(Isc,percentiage);
%%%%%%%% Visualize %%%%%%%%
pic = 3;
Draw(Isc_c,locs,p_value,alpha,example,pic,[' on All Branches']);

%%%%%%%% Preprocessing %%%%%%%%
data_s1 = data_s(branch==2,:);
%%%%%%%% Calculating Iscore %%%%%%%%
Isc1 = Calculate(data_s1,bin,DG,2);
%%%%%%%% Hypothesis Testing %%%%%%%%
[Isc1_c,p1_value,locs1] = Hypothesis_Testing(Isc1,percentiage);
%%%%%%%% Visualize %%%%%%%%
pic = 5;
Draw(Isc1_c,locs1,p1_value,alpha,example,pic,[' on Branch without Subbranch']);

%%%%%%%% Preprocessing %%%%%%%%
data_s2 = data_s(branch==3|branch==4|branch==5,:);
%%%%%%%% Calculating Iscore %%%%%%%%
Isc2 = Calculate(data_s2,bin,DG,2);
%%%%%%%% Hypothesis Testing %%%%%%%%
[Isc2_c,p2_value,locs2] = Hypothesis_Testing(Isc2,percentiage);
%%%%%%%% Visualize %%%%%%%%
pic = 7;
Draw(Isc2_c,locs2,p2_value,alpha,example,pic,[' on Branch with Subbranch']);