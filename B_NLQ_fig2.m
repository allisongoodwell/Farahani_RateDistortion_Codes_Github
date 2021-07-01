%03032020 n_level Quantization, Temperature
% Illustration of rate-distortion theory applied to air temperature
% to quantize the data to different precisions. N=2,3,4,5 levels of
% quantization
% Find the representation points and final decision threshold
% Find Representation points calculated with the fixed binning method.
% Figure 2 in paper
clear
close all


% Read forcing data excel file and save each variable
    % xlsx_filename = 'SoyMaize_GooseCreek_Forcing2018.xlsx';
    % sheet_number =1;
    % data_range = 'K14306:k24001';
    % T = xlsread(xlsx_filename,sheet_number,data_range);
    
load T    %load Temperature
T = sort(T) ; %
N=5; % Levels of quantization N
%%
[N1,edges] = histcounts(T,'Normalization','probability'); %Determine probability mass function (PMF)of focing data

for n_threshold = 1:4 % Levels of quantization N = n_threshold+1
    [pdf_fixedbin, Coords] = compute_pdf(T,n_threshold+1); % add points showing means for N-fixed binnig method
    
    subplot(2,2,n_threshold)
    hold on
    p1 = histogram(T,'Normalization','pdf');
    p4 = plot(Coords,0,'r*','linewidth',8);
    
    set(gca,'FontSize',20)    
    xlabel({'Temperature, \circ C'}, 'FontSize',20)
    ylabel('Density', 'FontSize',20)
    grid
    
    fprintf('working on threshold case: %d \n', n_threshold)
    cnt = ones(n_threshold,1) ;
    iteration = ones(n_threshold,1)*1e10 ;
    error = ones(n_threshold,1)*inf ;
    tol =ones(n_threshold,1)* .0001;
    
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
        
        equal_range = linspace(min(T),max(T),n_threshold+1) ;
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
        diff = abs(ave_Ex_T - t).^2 ;
        
        if  diff <= error
            error = diff ;
            t_final = t ; %find final thresholds
            
            mean_final = Ex_T; %find final corresponding representation levels
            
        end
        
        cnt = cnt + 1 ;
    end
    
    % plot final thresholds and means
    Titleindex=["N=2" "N=3" "N=4" "N=5"];
    subplot(2,2,n_threshold) 
    for i = 1:n_threshold
        for j= 1:n_threshold+1
            hold on
            p2 = xline(t_final(i,:),'--','LineWidth',2);
            p3 = xline(mean_final(j,:),'r','LineWidth',2);
            title (Titleindex(n_threshold), 'FontSize',20);
            hold off
            
        end
        fprintf('after %d iterations, error = %5.5f, VPD = % 5.2fkPa\n', cnt(1,1), error(i,1),t_final(i))
        legend([p1 p2 p3],{'PMF','Final Threshold','Mean'},'Location','northeast', 'FontSize',15)

    end
    fprintf('mean = % 5.2fkPa\n', mean_final)
end
%%
