%%%%%%%% Load data %%%%%%%%
load('ESC_qpcr_data.mat');
load('ESC_qpcr_DPT.mat');
%%%%%%%% Parameter %%%%%%%%
example = 4;
bin = 40;
percentiage = 0.75;
alpha = 0.05;
%%%%%%%% Preprocessing %%%%%%%%
data_s = Preprocessing(data,DPT);
%%%%%%%% Select DNB Group %%%%%%%%
DG = Selecting(data_s,bin,percentiage,alpha);
%%%%%%%% Calculating Iscore %%%%%%%%
Isc = Calculate(data_s,bin,DG,1);
%%%%%%%% Hypothesis Testing %%%%%%%%
[Isc_c,p_value,locs] = Hypothesis_Testing(Isc,percentiage);
%%%%%%%% Visualize %%%%%%%%
pic = 2;
Draw(Isc_c,locs,p_value,alpha,example,pic,[]);
