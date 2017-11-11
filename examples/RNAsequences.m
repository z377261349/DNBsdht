%%%%%%%% Load data %%%%%%%%
load('trajectory_RNAseq.mat');                 %%trajectory
load('branches_RNAseq.mat');                   %%branches
data = csvread('sample_scseq_data.csv',1,1);   %%sample
%%%%%%%% Parameter %%%%%%%%
example = 2;
bin = 20;
percentiage = 0.75;
alpha = 0.05;
%%%%%%%% Preprocessing %%%%%%%%
sigma = 10;
[T,phi0] = transition_matrix(data,'nn',[30,10]);
root = 3416;
branching = 1;
[M, tips] = dpt_input(T, phi0, branching, 'maxdptdist',root);
[Branch,DPT]=dpt_analyse(M,branching,tips);
[phi, lambda] = diffusionmap.eig_decompose_normalized(T,4);

BrcolorOrder=[0 0 0;
              0 0 1;
              1 0 1;
              1 0 0];
figure
subplot(1,2,1)
colormap jet
scatter(phi(:,2),phi(:,3),20,DPT,'fill');
xlabel('2nd diffusion component')
ylabel('3rd diffusion component')
title ('Colored by DPT Distance')
subplot(1,2,2)
scatter(phi(:,2),phi(:,3),20,Branch,'fill');
xlabel('2nd diffusion component')
ylabel('3rd diffusion component')
title('Colored by Branches')
saveas(gcf,'Figure_2_1','epsc')

data_s = Preprocessing(data,DPT);
%%%%%%%% Select DNB Group %%%%%%%%
DG = Selecting(data_s,bin,percentiage,alpha);
%%%%%%%% Calculating Iscore %%%%%%%%
Isc = Calculate(data_s,bin,DG,2);
%%%%%%%% Hypothesis Testing %%%%%%%%
[Isc_c,p_value,locs] = Hypothesis_Testing(Isc,percentiage);
%%%%%%%% Visualize %%%%%%%%
pic = 3;
Draw(Isc_c,locs,p_value,alpha,example,pic,[]);