% 3/26/2021: Allison making some edits
% goals: look between different quantization cases
% look at error between regular model and quant cases (rather than performance compared to observed data)
% Quantized versus full model performance for all cases,
% defined as the quantized diurnal error divided by full model diurnal error

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
steps_per_day = 24; %AEG: changed from 4*24 to hourly step
step = 1/steps_per_day;
decimal_hour_vect = doy_Model-floor(doy_Model);
hours_vect = 0:1:23;



%loop through result cases
cvect = jet(7); %color for each quantization category

cvecthour = hsv(24); %color for each hour

%matrices for whole-big-matrix of all diurnal results...for all cases
LE_MQ_vs_model_error_ALL=[];
H_MQ_vs_model_error_ALL=[];
Fc_MQ_vs_model_error_ALL=[];

LE_MQ_error_diurnal_ALL =[];
H_MQ_error_diurnal_ALL =[];
Fc_MQ_error_diurnal_ALL =[];

ctplot =1;
for q = 1:7 %7 "categories of quantization
    
    if q ==1
        case_folder = 'Results_addedfraction/Ta only';
    elseif q==2
        case_folder = 'Results_addedfraction/U only';
    elseif q ==3
        case_folder = 'Results_addedfraction/VPD only';
    elseif q ==4
        case_folder = 'Results_addedfraction/Ta_U';
    elseif q == 5
        case_folder = 'Results_addedfraction/Ta_VPD';
    elseif q==6
        case_folder = 'Results_addedfraction/VPD_U';
    elseif q ==7
        case_folder = 'Results_addedfraction/Ta_U_VPD';
    end
    
    filePattern = fullfile(case_folder, '*.mat');
    theFiles = dir(filePattern);

    load(fullfile(theFiles(1).folder, theFiles(1).name))
    load(fullfile(theFiles(2).folder, theFiles(2).name))
    load(fullfile(theFiles(3).folder, theFiles(3).name))
    
    num_file = size(Fc_MQ,2);
    cvectfile = jet(num_file); %color map for different models within each case...
    
    
    %this loads Fc_MQ, H_MQ and LE_MQ, these already include
    %corn and soy models combined
    
    
    for lehfc = 1:3 %loop through LE, H, and Fc output variables
        ct=1; %reset time counter
        
        MQ_hourly=[];
        MQ_error_diurnal=[];
        MQ_vs_model_error=[];
        
        
        for hr = 0:step:(1-step)   %loop of hourly segments
            indices = find(decimal_hour_vect>=hr-.001 & decimal_hour_vect<=hr+step+.001);
            
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
                MQ_vs_model_error(ct,i) = sqrt(mean((Model_ind - MQ_ind(:,i)).^2));
            end
            
            
            ct=ct+1; %increment time (hour of day)
        end
        
        %now need to attribute errors for different variables
        if lehfc ==1
            MQ_error_diurnal_LE = MQ_error_diurnal;
            Model_error_diurnal_LE = Model_error_diurnal;
            LE_RealData_hourly = RealData_hourly;
            LE_Model_hourly = Model_hourly;
            LE_MQ_hourly = MQ_hourly;
            LE_MQ_vs_model_error = MQ_vs_model_error;
            
            %allocate big matices to keep all results...
            LE_MQ_vs_model_error_ALL = [LE_MQ_vs_model_error_ALL MQ_vs_model_error];
            LE_MQ_error_diurnal_ALL = [LE_MQ_error_diurnal_ALL MQ_error_diurnal./Model_error_diurnal];
        elseif lehfc ==2
            MQ_error_diurnal_H = MQ_error_diurnal;
            Model_error_diurnal_H = Model_error_diurnal;
            H_RealData_hourly = RealData_hourly;
            H_Model_hourly = Model_hourly;
            H_MQ_hourly = MQ_hourly;
            H_MQ_vs_model_error = MQ_vs_model_error;
            
            H_MQ_vs_model_error_ALL = [H_MQ_vs_model_error_ALL MQ_vs_model_error];
            H_MQ_error_diurnal_ALL = [H_MQ_error_diurnal_ALL MQ_error_diurnal./Model_error_diurnal];
        else
            MQ_error_diurnal_Fc = MQ_error_diurnal;
            Model_error_diurnal_Fc = Model_error_diurnal;
            Fc_RealData_hourly = RealData_hourly;
            Fc_Model_hourly = Model_hourly;
            Fc_MQ_hourly = MQ_hourly;
            Fc_MQ_vs_model_error = MQ_vs_model_error;
            
            Fc_MQ_vs_model_error_ALL = [Fc_MQ_vs_model_error_ALL MQ_vs_model_error];
            Fc_MQ_error_diurnal_ALL = [Fc_MQ_error_diurnal_ALL MQ_error_diurnal./Model_error_diurnal];
        end
        
    ctplot=ctplot+1;
    end %end of loop through variables
    
end

%%
figure(1)
%for all categories, quantized diurnal error divided by full model diurnal error
%here, greater than 1 --> quantized model has worse performance
% less than 1 --> quantized model has lower error (better performance)
subplot(3,1,1)
imagesc(LE_MQ_error_diurnal_ALL)
xticks(cumsum([0 4 4 4 16 16 16]))
caxis([.8 1.2])
set(gca,'FontSize',10);
ylabel('Hour of day', 'FontSize',10)
title('Quantized model performance (P_{Q,t}): LE', 'FontSize',10)
colorbar;
colormap(jet(256));


subplot(3,1,2)
imagesc(H_MQ_error_diurnal_ALL)
xticks(cumsum([0 4 4 4 16 16 16]))
caxis([.8 1.2])
set(gca,'FontSize',10);
ylabel('Hour of day', 'FontSize',10)
title('Quantized model performance (P_{Q,t}): SH', 'FontSize',10)
colorbar

subplot(3,1,3)
imagesc(Fc_MQ_error_diurnal_ALL)
xticks(cumsum([0 4 4 4 16 16 16]))
caxis([.8 1.2])
set(gca,'FontSize',10);
ylabel('Hour of day', 'FontSize',10)
title('Quantized model performance (P_{Q,t}): Fc', 'FontSize',10)
colorbar

