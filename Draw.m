function [] = Draw(seq_c,locs,p_value,alpha,example,pic,title)
%%%%% Draw pictures
    m = size(locs,1);
    [~,ind] = sort(p_value);
    [~,ind_c] = sort(ind);
    k = max(ind_c(p_value<(ind_c.*alpha/m)));
    if(isempty(k))
        rej = false(m,1);
    else
        rej = logical(ind_c<=k);
        %rej = logical(p_value<alpha);
    end
    
    color = zeros(size(p_value,1),3);
    color(p_value>=(alpha/m),3) = 1;
    color(p_value<(alpha/m),1) = 1;
    figure;
    subplot(2,1,1);
    plot(seq_c);
    hold on
    stem(locs(rej),seq_c(locs(rej)),':*r','BaseValue',-5);
    stem(locs(~rej),seq_c(locs(~rej)),':.b','BaseValue',-5);
    ylabel({'Normalized';'Iscore'})
    xlabel('(a)')
    subplot(2,1,2);
    stem(locs(~rej),p_value(~rej),':ob','filled','BaseValue',1);
    hold on
    stem(locs(rej),p_value(rej),':or','filled','BaseValue',1);
    ylabel('p-value')
    xlabel({'bin';'(b)'})
    suptitle(['Hypotheses Testing of Each Local Maximum',title])
    saveas(gcf,['Figure_',num2str(example),'_',num2str(pic)],'epsc')
end