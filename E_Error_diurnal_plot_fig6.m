%11162020 
% Add Soy and Maize fraction
% MLCan diurnal error plotting
%Diurnal quantized model results for {LE},{SH}, and {Fc}.
% left panel indicates comparison between observation, full model and quantized model(Case 1 and Case 7) diurnal cycle,
% middle panel indicates the diurnal cycle error (RMSE) of full model and quantized model for Case 1 and Case 7 
% right panel indicates the difference between quantized model (Case 1) and full model RMSE of {LE},{SH}, and {Fc}.
clear
close all
%% Model rsults for each spicies
Soy_result_Model = load ('soybean_3.9.2020.mat','LE_eco_store','H_eco_store','Fc_eco_store','doy');
doy_Model = Soy_result_Model.doy;
% DOY 149-250
LE_Soy_Model = Soy_result_Model.LE_eco_store';
H_Soy_Model = Soy_result_Model.H_eco_store';
Fc_Soy_Model = Soy_result_Model.Fc_eco_store';

load Maize_result_Model.mat
% DOY 149-250
LE_Maize_Model = Maize_result_Model.LE_eco_store';
H_Maize_Model = Maize_result_Model.H_eco_store';
Fc_Maize_Model = Maize_result_Model.Fc_eco_store';

%% Quantized Model rsults for each spicies
case_folder = 'Results_addedfraction/Ta only';
filePattern = fullfile(case_folder, '*.mat');
theFiles = dir(filePattern);

Fc_MQ = load(fullfile(theFiles(1).folder, theFiles(1).name)).Fc_MQ;
H_MQ = load(fullfile(theFiles(2).folder, theFiles(2).name)).H_MQ;
LE_MQ = load(fullfile(theFiles(3).folder, theFiles(3).name)).LE_MQ;
% load(fullfile(theFiles(1).folder, theFiles(1).name))
% load(fullfile(theFiles(2).folder, theFiles(2).name))
% load(fullfile(theFiles(3).folder, theFiles(3).name))

%% Add fraction weight for results
load cornsoy_fractions.mat
Maize_fraction = data.Maize(73633:83328,:)./100;    %for DOY 150-250 year 2018
Soy_fraction = data.Soybean(73633:83328,:)./100;

% model result
LE_Model = LE_Soy_Model.*Soy_fraction + LE_Maize_Model.*Maize_fraction;
H_Model = H_Soy_Model.*Soy_fraction + H_Maize_Model.*Maize_fraction;
Fc_Model = Fc_Soy_Model.*Soy_fraction + Fc_Maize_Model.*Maize_fraction;


%% Direct measurment
load Directmeasurment.mat 
Fc = Directmeasurment(:,1);
Fc = Fc.*10^6./44009.5; %change unit from mg/(m^2 s) to Umol/(m^2 s) 
                        %1 mole co2 = 44009.5 mgram %1 mol = 10^6 umol
LE = Directmeasurment(:,2);
H = Directmeasurment(:,3);

%save original data
Fc_orig=Fc;
LE_orig = LE;
H_orig=H;
%remove outliers and make them nans
Fc(Fc<-30 | Fc>10)=nan;
LE(LE<-100 | LE>500)=nan;
H(H<-120 | H>300)=nan;

% Interpolation for missing data
Data = [Fc,LE,H];
for i = 2:length(Data)
    for k = 1:3
        if isnan(Data(i,k)) %if data is nan due to a gap
            Data_previous = Data(i-1,k); %last point
            for jj = 1:3000 %just a large number, longer than any gap
                Data_next = Data(i+jj,k); %next point: need to find a non-nan value
                if ~isnan(Data_next)
                    gap_length = jj+1;
                    break 
                end
            end            
            if gap_length <= 8 % if gap is less than 2 hrs,do linear interpolation, otherwise leave as nan 
                slope = (Data_next-Data_previous)./gap_length;
                Data(i,k)= Data_previous + slope;
            end
        end
    end
end
Fc_RealData = Data(:,1);
LE_RealData = Data(:,2);
H_RealData = Data(:,3);

%% Diurnal cycle
steps_per_day = 4*24; %15 minute data, each day
step = 1/steps_per_day; %15 minute step
decimal_hour_vect = doy_Model-floor(doy_Model);
hours_vect = 0:.25:23.75;
num_file = size(Fc_MQ,2);

for lehfc = 1:3 %loop through LE, H, and Fc output variables
    ct=1; %reset time counter
    
    for hr = 0:step:(1-step)   %loop of 15-minute segments
        indices = find(decimal_hour_vect>=hr-.001 & decimal_hour_vect<=hr+.001);

        if lehfc ==1
            RealData_ind = LE_RealData(indices);
            Model_ind = LE_Model(indices);
            MQ_ind = LE_MQ(indices,:);
        elseif lehfc ==2
            RealData_ind = H_RealData(indices);
            Model_ind = H_Model(indices);
            MQ_ind = H_MQ(indices,:);
        else
            RealData_ind = Fc_RealData(indices);
            Model_ind = Fc_Model(indices);
            MQ_ind = Fc_MQ(indices,:);
        end
        
        ind_var = find(isnan(RealData_ind)); %list of nan indexes in observed data
        
        RealData_ind(ind_var)=[]; %remove nan index values
        Model_ind(ind_var)=[];
        MQ_ind(ind_var,:)=[];
        
        RealData_hourly(ct)= mean(RealData_ind);
        Model_hourly(ct)= mean(Model_ind);
        Model_error_diurnal(ct,:) = sqrt(mean((RealData_ind - Model_ind).^2));

        for i=1:num_file %loop of quantization cases, find mean at that time of day and RMSE
            MQ_hourly(ct,i)= mean(MQ_ind(:,i));
            MQ_error_diurnal(ct,i) = sqrt(mean((RealData_ind - MQ_ind(:,i)).^2));
        end
        
        
        ct=ct+1; %increment time (hour of day)
    end
    %keep the RMSE for each time point, but average every hour
    Model_AVg_1 =  mean(reshape(Model_error_diurnal, 4, []));
    for i=1:num_file %loop of quantization cases, find mean for every 1 hour and RMSE       
        MQ_AVg_1(:,i) =  mean(reshape(MQ_error_diurnal(:,i), 4, []));
    end
    %now need to attribute errors for different variables
    if lehfc ==1
        MQ_error_diurnal_LE = MQ_error_diurnal;
        Model_error_diurnal_LE = Model_error_diurnal;       
        MQ_error_avg_LE = MQ_AVg_1;
        Model_error_avg_LE = Model_AVg_1;        
        LE_RealData_hourly = RealData_hourly;
        LE_Model_hourly = Model_hourly;
        LE_MQ_hourly = MQ_hourly;
    elseif lehfc ==2
        MQ_error_diurnal_H = MQ_error_diurnal;
        Model_error_diurnal_H = Model_error_diurnal;
        MQ_error_avg_H = MQ_AVg_1;
        Model_error_avg_H = Model_AVg_1;
        H_RealData_hourly = RealData_hourly;
        H_Model_hourly = Model_hourly;
        H_MQ_hourly = MQ_hourly;
    else
        MQ_error_diurnal_Fc = MQ_error_diurnal;
        Model_error_diurnal_Fc = Model_error_diurnal;
        MQ_error_avg_Fc = MQ_AVg_1;
        Model_error_avg_Fc = Model_AVg_1;
        Fc_RealData_hourly = RealData_hourly;
        Fc_Model_hourly = Model_hourly;
        Fc_MQ_hourly = MQ_hourly;
    end
    
end %end of loop through variables
%% Case 7: Quantized Model rsults for each spicies
case_folder_7 = 'Results_addedfraction/Ta_U_VPD';
filePattern_7 = fullfile(case_folder_7, '*.mat');
theFiles_7 = dir(filePattern_7);

Fc_MQ_7 = load(fullfile(theFiles_7(1).folder, theFiles_7(1).name)).Fc_MQ;
H_MQ_7 = load(fullfile(theFiles_7(2).folder, theFiles_7(2).name)).H_MQ;
LE_MQ_7 = load(fullfile(theFiles_7(3).folder, theFiles_7(3).name)).LE_MQ;
% Diurnal cycle
num_file_7 = size(Fc_MQ_7,2);

for lehfc = 1:3 %loop through LE, H, and Fc output variables
    ct=1; %reset time counter
    
    for hr = 0:step:(1-step)   %loop of 15-minute segments
        indices = find(decimal_hour_vect>=hr-.001 & decimal_hour_vect<=hr+.001);

        if lehfc ==1
            RealData_ind = LE_RealData(indices);
            Model_ind = LE_Model(indices);
            MQ_ind = LE_MQ_7(indices,:);
        elseif lehfc ==2
            RealData_ind = H_RealData(indices);
            Model_ind = H_Model(indices);
            MQ_ind = H_MQ_7(indices,:);
        else
            RealData_ind = Fc_RealData(indices);
            Model_ind = Fc_Model(indices);
            MQ_ind = Fc_MQ_7(indices,:);
        end
        
        ind_var = find(isnan(RealData_ind)); %list of nan indexes in observed data
        
        RealData_ind(ind_var)=[]; %remove nan index values
        Model_ind(ind_var)=[];
        MQ_ind(ind_var,:)=[];
        
        RealData_hourly(ct)= mean(RealData_ind);
        Model_hourly(ct)= mean(Model_ind);
        Model_error_diurnal(ct,:) = sqrt(mean((RealData_ind - Model_ind).^2));

        for i=1:num_file_7 %loop of quantization cases, find mean at that time of day and RMSE
            MQ_hourly(ct,i)= mean(MQ_ind(:,i));
            MQ_error_diurnal(ct,i) = sqrt(mean((RealData_ind - MQ_ind(:,i)).^2));
        end
        
        
        ct=ct+1; %increment time (hour of day)
    end
    %keep the RMSE for each time point, but average every hour
    Model_AVg_1 =  mean(reshape(Model_error_diurnal, 4, []));
    for i=1:num_file_7 %loop of quantization cases, find mean for every 1 hour and RMSE       
        MQ_AVg_7(:,i) =  mean(reshape(MQ_error_diurnal(:,i), 4, []));
    end
    %now need to attribute errors for different variables
    if lehfc ==1
        MQ_error_diurnal_LE_7 = MQ_error_diurnal;
        MQ_error_avg_LE_7 = MQ_AVg_7;
        LE_MQ_hourly_7 = MQ_hourly;
    elseif lehfc ==2
        MQ_error_diurnal_H_7 = MQ_error_diurnal;
        MQ_error_avg_H_7 = MQ_AVg_7;
        H_MQ_hourly_7 = MQ_hourly;
    else
        MQ_error_diurnal_Fc_7 = MQ_error_diurnal;
        MQ_error_avg_Fc_7 = MQ_AVg_7;
        Fc_MQ_hourly_7 = MQ_hourly;
    end
    
end %end of loop through variables

%% Some plotting
fig = figure;
for i = 1:num_file
    
    subplot(3,3,1)
    plot([1:24],mean(reshape(LE_RealData_hourly(1,:), 4, [])),'color','k','LineWidth',1.5)
    hold on
    set(gca,'FontSize',11,'XTick',(0:3:24))
    plot([1:24],mean(reshape(LE_Model_hourly(1,:), 4, [])),'color','b','LineWidth',1.5)
    plot([1:24],mean(reshape(LE_MQ_hourly(:,i), 4, [])),'color','r','LineWidth',1.5)
    plot([1:24],mean(reshape(LE_MQ_hourly_7(:,i), 4, [])),'color','g','LineWidth',1.5)
    legend({'Obs','Full Model','QM_{Case1}', 'QM_{Case7}'}, 'FontSize',8,'orientation','horizontal')
    title('Diurnal fluxes', 'FontSize',12)
    ylabel('LE (W/m^2)', 'FontSize',12)
    xlim([1 24])       
    subplot(3,3,4)
    plot([1:24],mean(reshape(H_RealData_hourly(1,:), 4, [])),'color','k','LineWidth',1.5)
    hold on
    set(gca,'FontSize',11,'XTick',(0:3:24))
    plot([1:24],mean(reshape(H_Model_hourly(1,:), 4, [])),'color','b','LineWidth',1.5)
    plot([1:24],mean(reshape(H_MQ_hourly(:,i), 4, [])),'color','r','LineWidth',1.5)
    plot([1:24],mean(reshape(H_MQ_hourly_7(:,i), 4, [])),'color','g','LineWidth',1.5)
  
    ylabel('SH (W/m^2)', 'FontSize',12)
    xlim([1 24])   
    subplot(3,3,7)
    plot([1:24],mean(reshape(Fc_RealData_hourly(1,:), 4, [])),'color','k','LineWidth',1.5)
    hold on
    set(gca,'FontSize',11,'XTick',(0:3:24))
    plot([1:24],mean(reshape(Fc_Model_hourly(1,:), 4, [])),'color','b','LineWidth',1.5)
    plot([1:24],mean(reshape(Fc_MQ_hourly(:,i), 4, [])),'color','r','LineWidth',1.5)
    plot([1:24],mean(reshape(Fc_MQ_hourly_7(:,i), 4, [])),'color','g','LineWidth',1.5)   
    xlabel('Hour', 'FontSize',12)
    ylabel('Fc (\mumol/m^2s)', 'FontSize',12)
    xlim([1 24])
    
    subplot(3,3,2)
    plot([1:24],Model_error_avg_LE,'Color','b','LineWidth',1.5)
    hold on
    plot([1:24],MQ_error_avg_LE(:,i),'color','r','LineWidth',1.5)
    plot([1:24],MQ_error_avg_LE_7(:,i),'color','g','LineWidth',1.5)
    set(gca,'FontSize',11,'XTick',(0:3:24))
    legend({'Full Model','QM_{Case1}','QM_{Case7}'}, 'FontSize',8,'orientation','horizontal')
    title('Diurnal RMSE', 'FontSize',12)
    xlim([1 24])
    
    subplot(3,3,5)
    plot([1:24],Model_error_avg_H,'Color','b','LineWidth',1.5)    
    hold on
    plot([1:24],MQ_error_avg_H(:,i),'color','r','LineWidth',1.5)
    plot([1:24],MQ_error_avg_H_7(:,i),'color','g','LineWidth',1.5)
    set(gca,'FontSize',11,'XTick',(0:3:24))
    xlim([1 24])  
    
    subplot(3,3,8)
    plot([1:24],Model_error_avg_Fc,'Color','b','LineWidth',1.5)
    hold on
    set(gca,'FontSize',11,'XTick',(0:3:24))
    plot([1:24],MQ_error_avg_Fc(:,i),'color','r','LineWidth',1.5)
    plot([1:24],MQ_error_avg_Fc_7(:,i),'color','g','LineWidth',1.5)
    xlabel('Hour', 'FontSize',12)
    xlim([1 24])
    
    subplot(3,3,3)   
    C = {'r';'g';'m'; 'c'} ;
    plot([1:24],(MQ_error_avg_LE(:,i)-Model_error_avg_LE'),'color',C{i},'LineWidth',1.5)
    hold on
    set(gca,'FontSize',11,'XTick',(0:3:24))
    legend({'N_T =2' 'N_T=3' 'N_T=4' 'N_T=5'}, 'FontSize',8,'orientation','horizontal')
    title('\Delta_{RMSE}', 'FontSize',10)
    xlim([1 24])    
    subplot(3,3,6)
    plot([1:24],(MQ_error_avg_H(:,i)-Model_error_avg_H'),'color',C{i},'LineWidth',1.5)
    hold on
    set(gca,'FontSize',11,'XTick',(0:3:24))
    xlim([1 24])  
    subplot(3,3,9)
    plot([1:24],(MQ_error_avg_Fc(:,i)-Model_error_avg_Fc'),'color',C{i},'LineWidth',1.5)
    hold on
    set(gca,'FontSize',11,'XTick',(0:3:24))
    xlabel('Hour', 'FontSize',12)
    xlim([1 24])

end

for i1 = 1:9
  subplot(3,3,i1)
  grid on
end