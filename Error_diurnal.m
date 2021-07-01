function [MQ_error_diurnal, Model_error_diurnal] =Error_diurnal(case_folder)
%used in Figure 9,10,11

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
% Results = 'Results_addedfraction/Ta only';
filePattern = fullfile(case_folder, '*.mat');
theFiles = dir(filePattern);

Fc_MQ = load(fullfile(theFiles(1).folder, theFiles(1).name)).Fc_MQ;
H_MQ = load(fullfile(theFiles(2).folder, theFiles(2).name)).H_MQ;
LE_MQ = load(fullfile(theFiles(3).folder, theFiles(3).name)).LE_MQ;

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

for i=1:num_file
    ct=1;
    for hr = 0:step:(1-step)      
        indices = find(decimal_hour_vect>=hr-.001 & decimal_hour_vect<=hr+.001);

        LE_RealData_ind = LE_RealData(indices);
        %remove NaN value
        ind_LE = zeros(1) ;
        cnt_LE = 1 ;
        for it_row = 1:size(LE_RealData_ind,1)           
            if ~isnan(LE_RealData_ind(it_row,1))
                ind_LE(cnt_LE,1) = it_row ;
                cnt_LE = cnt_LE + 1 ;
            end
        end
        LE_RealData_inds = LE_RealData_ind(ind_LE,:);
        
        LE_Model_ind = LE_Model(indices); 
        LE_Model_inds = LE_Model_ind(ind_LE,:); %remove same indices of NaN value in RealData from model data
        
        LE_MQ_ind(:,i) = LE_MQ(indices,i);
        LE_MQ_inds = zeros(size(ind_LE,1),1);
        LE_MQ_inds(:,i) = LE_MQ_ind(ind_LE,i); %remove same indices of NaN value in RealData from qantized data
        
        H_RealData_ind = H_RealData(indices);
        ind_H = zeros(1) ;
        cnt_H = 1 ;
        for it_row = 1:size(H_RealData_ind,1)           
            if ~isnan(H_RealData_ind(it_row,1))
                ind_H(cnt_H,1) = it_row ;
                cnt_H = cnt_H + 1 ;
            end
        end
        H_RealData_inds = H_RealData_ind(ind_H,:);
        
        H_Model_ind = H_Model(indices);
        H_Model_inds = H_Model_ind(ind_H,:);
        
        H_MQ_ind(:,i) = H_MQ(indices,i);
        H_MQ_inds = zeros(size(ind_H,1),1);
        H_MQ_inds(:,i) = H_MQ_ind(ind_H,i);
        
        Fc_RealData_ind = Fc_RealData(indices);
        ind_Fc = zeros(1) ;
        cnt_Fc = 1 ;
        for it_row = 1:size(Fc_RealData_ind,1)           
            if ~isnan(Fc_RealData_ind(it_row,1))
                ind_Fc(cnt_Fc,1) = it_row ;
                cnt_Fc = cnt_Fc + 1 ;
            end
        end
        Fc_RealData_inds = Fc_RealData_ind(ind_Fc,:);
        
        Fc_Model_ind = Fc_Model(indices);
        Fc_Model_inds = Fc_Model_ind(ind_Fc,:);
        
        Fc_MQ_ind(:,i) = Fc_MQ(indices,i);
        Fc_MQ_inds = zeros(size(ind_Fc,1),1);
        Fc_MQ_inds(:,i) = Fc_MQ_ind(ind_Fc,i);
        
        LE_RealData_hourly(ct)= mean(LE_RealData_inds);
        H_RealData_hourly(ct)= mean(H_RealData_inds);
        Fc_RealData_hourly(ct)= mean(Fc_RealData_inds);
        
        LE_Model_hourly(ct)= mean(LE_Model_inds);
        H_Model_hourly(ct)= mean(H_Model_inds);
        Fc_Model_hourly(ct)= mean(Fc_Model_inds);
        
        LE_MQ_hourly(ct,i)= mean(LE_MQ_inds(:,i));
        H_MQ_hourly(ct,i)= mean(H_MQ_inds(:,i));
        Fc_MQ_hourly(ct,i)= mean(Fc_MQ_inds(:,i)); 
        
        % Diurnal Error Calculation
        Model_error_diurnal_LE(ct,:) = sqrt(mean((LE_RealData_inds - LE_Model_inds).^2));
        Model_error_diurnal_H(ct,:) = sqrt(mean((H_RealData_inds - H_Model_inds).^2));
        Model_error_diurnal_Fc(ct,:) = sqrt(mean((Fc_RealData_inds - Fc_Model_inds).^2));
            
        MQ_error_diurnal_LE(ct,i) = sqrt(mean((LE_RealData_inds - LE_MQ_inds(:,i)).^2));
        MQ_error_diurnal_H(ct,i)  = sqrt(mean((H_RealData_inds - H_MQ_inds(:,i)).^2));
        MQ_error_diurnal_Fc(ct,i) = sqrt(mean((Fc_RealData_inds - Fc_MQ_inds(:,i)).^2));
       
    ct=ct+1;
    
    end
end
m =1;
for i =1:size(MQ_error_diurnal_LE,2)
    MQ_error_diurnal(:,:,m) = [MQ_error_diurnal_LE(:,i)  MQ_error_diurnal_H(:,i)  MQ_error_diurnal_Fc(:,i)];
    m = m+1;
end
Model_error_diurnal = [Model_error_diurnal_LE  Model_error_diurnal_H  Model_error_diurnal_Fc];
end
