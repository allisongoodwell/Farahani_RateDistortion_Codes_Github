%11162020
%Difference between quantized and full model error ($\Delta_{RMSE_t}$) of {LE}
% within five 20-day time windows through the growing season.
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

load(fullfile(theFiles(1).folder, theFiles(1).name))
load(fullfile(theFiles(2).folder, theFiles(2).name))
load(fullfile(theFiles(3).folder, theFiles(3).name))

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
n=20; %20 days, each window
steps_per_window = steps_per_day*n;
step = 1/steps_per_day; %15 minute step
% decimal_hour_vect = doy_Model-floor(doy_Model);
hours_vect = 0:.25:23.75;
num_file = size(Fc_MQ,2);
decimal_hour_vect_MQ = doy_Model(1:steps_per_window,1)-floor(doy_Model(1:steps_per_window,1));

%divide quantized values to "5" windows, each 20 days
for i=1:num_file
    for w = 1:5 %loop through 5 20-day windows
        LE_MQ_w(:,i,w) = LE_MQ((w-1)*steps_per_window+1:w*steps_per_window,i);
        H_MQ_w(:,i,w) = H_MQ((w-1)*steps_per_window+1:w*steps_per_window,i);
        Fc_MQ_w(:,i,w) = Fc_MQ((w-1)*steps_per_window+1:w*steps_per_window,i);
        
    end
end


%divide observation and model values to "5" windows, each 20 days
for w=1:5
    LE_RealData_w(:,w)=LE_RealData((w-1)*steps_per_window+1:w*steps_per_window,:);
    H_RealData_w(:,w)=H_RealData((w-1)*steps_per_window+1:w*steps_per_window,:);
    Fc_RealData_w(:,w)=Fc_RealData((w-1)*steps_per_window+1:w*steps_per_window,:);
    
    LE_Model_w(:,w)=LE_Model((w-1)*steps_per_window+1:w*steps_per_window,:);
    H_Model_w(:,w)=H_Model((w-1)*steps_per_window+1:w*steps_per_window,:);
    Fc_Model_w(:,w)=Fc_Model((w-1)*steps_per_window+1:w*steps_per_window,:);
    
end


for lehfc = 1:3 %loop through LE, H, and Fc output variables
    for i=1:num_file
        ct=1; %reset time counter
        for hr = 0:step:(1-step)   %loop of 15-minute segments
            indices = find(decimal_hour_vect_MQ>=hr-.001 & decimal_hour_vect_MQ<=hr+.001);
            for w=1:5
                if lehfc ==1
                    RealData_ind = LE_RealData_w(indices,w);
                    Model_ind = LE_Model_w(indices,w);
                    MQ_ind = LE_MQ_w(indices,i,w);
                elseif lehfc ==2
                    RealData_ind = H_RealData_w(indices,w);
                    Model_ind = H_Model_w(indices,w);
                    MQ_ind = H_MQ_w(indices,i,w);
                    
                else
                    RealData_ind = Fc_RealData_w(indices,w);
                    Model_ind = Fc_Model_w(indices,w);
                    MQ_ind = Fc_MQ_w(indices,i,w);
                    
                end
                
                ind_var = find(isnan(RealData_ind)); %list of nan indexes in observed data
                
                RealData_ind(ind_var)=[]; %remove nan index values
                Model_ind(ind_var)=[];
                MQ_ind(ind_var,:,:)=[];
                
                
                RealData_hourly(ct,w)= mean(RealData_ind);
                Model_hourly(ct,w)= mean(Model_ind);
                Model_error_diurnal(ct,w) = sqrt(mean((RealData_ind - Model_ind).^2));
                
                MQ_hourly(ct,i,w)= mean(MQ_ind);
                MQ_error_diurnal(ct,i,w) = sqrt(mean((RealData_ind - MQ_ind).^2));
                
            end
            
            ct=ct+1; %increment time (hour of day)
        end
        
        %keep the RMSE for each time point, but average every hour
        for w =1:5
            Model_AVg_1(:,w) =  mean(reshape(Model_error_diurnal(:,w), 4, []));
            MQ_AVg_1(:,i,w) =  mean(reshape(MQ_error_diurnal(:,i,w), 4, []));
        end
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
end

%end of loop through variables



%% Plot

fig = figure;
Titleindex=["W=1","W=2" "W=3" "W=4" "W=5"];
for i = 1:num_file
    for w =1:5
        subplot(1,5,w) ;
        set(gca,'FontSize',12,'XTick',(0:4:24))
        C = {'r';'g';'b';'k'} ;
        plot([1:24],MQ_error_avg_LE (:,i,w)-Model_error_avg_LE(:,w),'color',C{i},'LineWidth',1.2);
        grid on
        hold on
        
        title (Titleindex(w), 'FontSize',12);
        ylim([-15 25]);
        xlim([1 24]);
        legend({'N_T =2' 'N_T=3' 'N_T=4' 'N_T=5'}, 'FontSize',8)
        han=axes(fig,'visible','off');
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        ylabel(han, '\Delta_{RMSE_t}for LE', 'FontSize',13);
        xlabel(han,'Hour of day', 'FontSize',13);
       
    end
end