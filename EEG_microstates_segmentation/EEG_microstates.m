close all
clear all
%% INSTRUCTIONS
% mat file for each subject which contains:
%  
% 'data' ? The EEG (channels x time x trials).
% 'labels' ? Labels for each trial, denoting the stimulus type (1 = famous, 2 = unfamiliar, 3 = scrambled).
% 'time' ? Timepoints in ms. t=0 is when the stimulus (picture) is shown, and what happens before should not be included in your microstate segmentation.
% 'fs' ? sampling frequency
% 'prep_quality' ? vector with a letter for each trial denoting a subjective rating of the level of preprocessing (g = good, o = ok quality, b = bad quality).
% 'chanlocs' ? Struct with information on the positions of the electrodes. Is used in the topoplot function to plot the scalp maps of timepoints of the microstates.
%  
% Since you have downloaded EEGlab you also have the topoplot function. You can make a figure with it using:
% figure
% topoplot(data_timepoint, chanlocs);
%  
% To check if you understand how the data is structured you can see if you can plot the average scalp map for a specific time point, e.g., 250 ms after stimulus. You can compare it with the same timepoint in this video I made from averaging all trials for subject 3.
%% IMPORT FILE
load Sub01.mat

%% PLAY with the topoplot
%{
figure(1)
time_index = 1:size(data,2);
for t = time_index
topoplot(data(:,t,1),chanlocs,'style','map')
title(['Time ',num2str(time(t)),' ms'])
drawnow
end
%}
%% INITIALIZE PARAMETERS
trial =120;
V_t = data(:,:,trial);
N_s = size(data,1);
b = 50; 
N_T = size(data,2);

%% FIND the best number of microstates
N_mu_max = 8;
N_mu_min = 2;
N_mu_vector = N_mu_min:N_mu_max;
N_simulations = 10*3;
for N_sim = 1:N_simulations %simulations for averaging fit
    for i = 1:length(N_mu_vector)
        N_mu = N_mu_vector(i);
        [~,~,~,~,sigma2_mu(N_sim,i) ] = modified_Kmean(V_t,N_mu);
        sigma2_mcv(i) = sigma2_mu(N_sim,i)*((N_s-1)^-1 * (N_s-1-N_mu))^-2;
    end
    sigma2_mcv_mat(N_sim,:) = sigma2_mcv;
    [~,bestK_i] = min(sigma2_mcv);
    bestK(N_sim) = N_mu_vector(bestK_i);
end


%% SHOW the best K
figure(2)
histogram(bestK)
title(['\mu_k: ',num2str(mean(bestK)),', \sigma_k: ',num2str(std(bestK))]);

%% CLUSTERING
N_mu_best = round(mean(bestK));
[L_t,R2,R2_s,Gamma_k,sigma2_mu ] = modified_Kmean(V_t,N_mu_best);

%% FEATURE EXTRACTION
sampling_time = 1e-3*abs((time(2) - time(1)));
GFP_t = GFP(V_t) ;
[ average_lifespan,frequency,coverage,amplitude,transition_matrix,GEV ] = u_state_statistics( L_t,sampling_time,GFP_t );

%% EVOLVING USTATE IN TIME
N_mu_colors = rand(length(1:N_mu_best),3);
figure(3)
plot_ustate_gfp(L_t,GFP_t,time,N_mu_colors)
title('U-state evolution')



