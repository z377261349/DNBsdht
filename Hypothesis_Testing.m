function [seq_c,p_value,locs] = Hypothesis_Testing(seq,percentiage)
%%%%% Hypothesis Testing

%%%%% Centralize sequences
mean_hat = median(seq);
sigma_hat = 1/icdf('norm',percentiage,0,1)*mad(seq,1);
seq_c = (seq-mean_hat)/sigma_hat;

%%%%% Find peaks of sequences
[pks,locs] = findpeaks(seq_c);
p_value = 1-cdf('Normal',pks,0,1).^3;
end