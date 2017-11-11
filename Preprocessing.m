function [data_s] = Preprocessing(data,trajectory)
%%%%% sort cells
[~,I] = sort(trajectory);
data_s = data(I,:);
end