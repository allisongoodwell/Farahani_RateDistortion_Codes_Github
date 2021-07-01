function input_Q = Quantization_function(inupt_variable)
%22062020 N_level Qunatization
N=5;

[N1,edges] = histcounts(sort(inupt_variable),'Normalization','probability');

for n_threshold = 1:N-1 % Levels of quantization N = n_threshold+1
    
    fprintf('working on inupt_variable, threshold case: %d \n', n_threshold)
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
        % Choose a randome threshold: Set initial threshold guesses
        % (T_1,...,T_N-1)
        if cnt == 1
            for it_rand = 1:n_threshold
                t(it_rand,1) = rand(1)*(equal_range(it_rand+1) - equal_range(it_rand)) + ...
                    equal_range(it_rand) ;
            end
        else
            t =  ave_Ex_T ;
        end
        
        % Find values X_i and relative probabilities Pr(X_i|T_j-1<X_i<T_j);
        % for j=1,...,N
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
        rel_prob_edge = (edge.*N1)./(sum(edge.*N1,2)) ;
        
        % Calculate the expected values of X_i between two sccessive
        % thresholds (EV_j)
        Ex_T = sum(vals_edge.*rel_prob_edge,2);
        
        % Find new thresholds: T_j-1= (EV_j+EV_(j-1))/2
        for i = 1:n_threshold
            ave_Ex_T(i,:) = (Ex_T(i+1,:)+Ex_T(i,:))./2 ;
        end
        % Ecpected distortion: Error between new and old thresholds
        diff = abs(ave_Ex_T - t) ;
        
        if  diff <= error
            error = diff ;
            t_final = t ;  %find final thresholds
            
            mean_final = Ex_T; %find final corresponding representation levels
            
        end
        
        cnt = cnt + 1 ;
    end
    % simplify individual forcing time-series data to N levels of quantization
    for i = 1:size(inupt_variable)
        if n_threshold == 1
            if inupt_variable(i) <= t_final(n_threshold)
                inpt_variable_1(i,:) = mean_final(n_threshold);
            else
                inpt_variable_1(i,:) = mean_final(n_threshold+1);
            end

            
        elseif n_threshold ==2
            if inupt_variable(i) <= t_final(n_threshold-1)
                inpt_variable_2(i,:) = mean_final(n_threshold-1);
            elseif inupt_variable(i)> t_final(n_threshold-1) && inupt_variable(i)<t_final(n_threshold)
                inpt_variable_2(i,:) = mean_final(n_threshold);
            else
                inpt_variable_2(i,:) = mean_final(n_threshold+1);
            end

            
        elseif n_threshold == 3
            if inupt_variable(i) <= t_final(n_threshold-2)
                inpt_variable_3(i,:) = mean_final(n_threshold-2);
            elseif inupt_variable(i)> t_final(n_threshold-2) && inupt_variable(i)<t_final(n_threshold-1)
                inpt_variable_3(i,:) = mean_final(n_threshold-1);
            elseif inupt_variable(i)> t_final(n_threshold-1) && inupt_variable(i)<t_final(n_threshold)
                inpt_variable_3(i,:) = mean_final(n_threshold);
            else
                inpt_variable_3(i,:) = mean_final(n_threshold+1);
            end
            
        elseif n_threshold == 4
            if inupt_variable(i) <= t_final(n_threshold-3)
                inpt_variable_4(i,:) = mean_final(n_threshold-3);
            elseif inupt_variable(i)> t_final(n_threshold-3) && inupt_variable(i)<t_final(n_threshold-2)
                inpt_variable_4(i,:) = mean_final(n_threshold-2);
            elseif inupt_variable(i)> t_final(n_threshold-2) && inupt_variable(i)<t_final(n_threshold-1)
                inpt_variable_4(i,:) = mean_final(n_threshold-1);
            elseif inupt_variable(i)> t_final(n_threshold-1) && inupt_variable(i)<t_final(n_threshold)
                inpt_variable_4(i,:) = mean_final(n_threshold);
            else
                inpt_variable_4(i) = mean_final(n_threshold+1);
            end
            
        end
    end    
end
input_Q = [inpt_variable_1 inpt_variable_2 inpt_variable_3 inpt_variable_4]; %quantized forcing vaeiable
end
