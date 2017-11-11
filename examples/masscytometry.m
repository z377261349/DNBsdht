%%%%%%%% Load data %%%%%%%%
load('trajectory.mat');
load('cluster.mat');
data = csvread('masscytometry.csv',1,1);
%%%%%%%% Parameter %%%%%%%%
example = 1;
bin = 40;
percentiage = 0.75;
alpha = 0.05;
%%%%%%%% Preprocessing %%%%%%%%
data_s = Preprocessing(data,trajectory);
%%%%%%%% Select DNB Group %%%%%%%%
DG = true(1,size(data_s,2));
%%%%%%%% Calculating Iscore %%%%%%%%
Isc = Calculate(data_s,bin,DG,1);
%%%%%%%% Hypothesis Testing %%%%%%%%
[Isc_c,p_value,locs] = Hypothesis_Testing(Isc,percentiage);
%%%%%%%% Visualize %%%%%%%%
pic = 1;
Draw(Isc_c,locs,p_value,alpha,example,pic,[]);
