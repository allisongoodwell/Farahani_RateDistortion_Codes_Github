% Differences between the entropy of the quantized forcing data 
% of Ta, U and VPD into N bins using fixed binning and Lloyd algorithm.

function [Hx_input,  Hx_input_Q] = diffentropy_function(inupt_variable) 

N=5;

[N1,edges] = histcounts(sort(inupt_variable),'Normalization','probability');

%%
for n_threshold = 1:N-1
    [pdf_fixedbin, Coords] = compute_pdf(inupt_variable,n_threshold+1); % add points showing means for N-fixed binnig method
  
    fprintf('working on threshold case: %d \n', n_threshold)
    cnt = ones(n_threshold,1) ;
    iteration = ones(n_threshold,1)*1e10 ;
    error = ones(n_threshold,1)*inf ;
    tol = ones(n_threshold,1)* .0001;
    
    t_final = zeros(n_threshold,1) ;
    diff = zeros(n_threshold,1) ;
    ave_Ex_T = zeros(n_threshold,1);
    t = zeros(n_threshold,1) ;
    mean_final = zeros(n_threshold+1,1) ;
    
    % Calculate the center of edges (i.e., average of each pair of edges)
    ave_edge = zeros(1) ;
    for it_edge = 1:size(edges,2) - 1
        ave_edge(1,it_edge) = mean([edges(it_edge) , edges(it_edge+1)]) ;
    end
    
    while any(error >= tol) && any(cnt <= iteration)
        
        equal_range = linspace(min(inupt_variable),max(inupt_variable),n_threshold+1) ;
        % Choose a randome threshold
        if cnt == 1
            for it_rand = 1:n_threshold
                t(it_rand,1) = rand(1)*(equal_range(it_rand+1) - equal_range(it_rand)) + ...
                    equal_range(it_rand) ;
            end
        else
            t =  ave_Ex_T ;
        end
        
        % Find values and relative probabilities
        edge = zeros(n_threshold+1,it_edge);
        for it_rand_t = 1:n_threshold+1
            if it_rand_t == 1
                edge(it_rand_t,:) = ave_edge < t(it_rand_t) ;
            elseif 1 < it_rand_t && it_rand_t < n_threshold+1
                edge(it_rand_t,:) = all([ave_edge < t(it_rand_t) ; ave_edge >= t(it_rand_t-1)],1) ;
            elseif it_rand_t == n_threshold+1
                edge(it_rand_t,:) = ave_edge >= t(it_rand_t-1) ;
            end
        end
        vals_edge = edge.*ave_edge;
        
        rel_prob_edge = (edge.*N1)./(sum(edge.*N1,2));
       
        
        % The expected input value
        Ex_T = sum(vals_edge.*rel_prob_edge,2);
        
        % new thresholds
        for i = 1:n_threshold
            ave_Ex_T(i,:) = (Ex_T(i+1,:)+Ex_T(i,:))./2;
        end
        % difference
        diff = abs(ave_Ex_T - t);
        
        if  diff <= error
            error = diff;
            t_final = t;
            pdf_input_Q = sum((edge.*N1),2);         
            mean_final = Ex_T;
            
        end
        
        cnt = cnt + 1;
    end
    
    %differences between the entropy of the quantized input data (fixed binning and Llyod algorithm).
    Hx_input(:,n_threshold) = compute_info_measures(pdf_fixedbin).Hx;
    Hx_input_Q(:,n_threshold) = compute_info_measures(pdf_input_Q).Hx;
 
end
end

