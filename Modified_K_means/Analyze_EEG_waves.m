%% Read some data
clear all %#ok<CLALL>
close all
[hdr, V_t] = edfread('S001R01.edf') ;
V_t = V_t(1:end-1,:);

%% Visualize EEG waves
figure(1)

hold on
max_potential_value = max(max(V_t(1:end,:)));
min_potential_value = min(min(V_t(1:end,:)));
peak2peak_max_amplitude = max_potential_value - min_potential_value;
for i = 1:size(V_t,1)-1
    plot(i*peak2peak_max_amplitude + V_t(i,:))
end
title('EEG waveforms')
hold off

N_T = size(V_t,2); %number of samples - number of time instants
N_s = size(V_t,1); %number of electrodes

% Optimize number of microstate
disp('Optimize number of microstate...')
K = 3:5;
for i=1:length(K)
    N_mu = K(i);
    [~,~,~,~,sigma2_mu ] = modified_Kmean(V_t,N_mu); 
    sigma2_mcv(i) = sigma2_mu*((N_s-1)^-1 * (N_s-1-N_mu))^-2;
end
[~,bestK] = min(sigma2_mcv);
bestK = K(bestK);
figure
plot(K,sigma2_mcv,'-o');
title('Predictive residual variance')
ylabel('\sigma^2_{mcv}');
xlabel('N_{\mu}')

% recompute for optimal number of microstate
[L_t,R2,R2_s,Gamma_k,sigma2_mu ] = modified_Kmean(V_t,bestK);
N_mu = bestK;

%% show us the microstates
figure
for k = 1:N_mu
    subplot(1,N_mu,k)
    microstate =  topographic_brain_map_64(Gamma_k{k}');
    imagesc(microstate, [0,max(max(microstate))]);
    title(['Microstate ',num2str(k)])
end
%% check results in real time
figure
for t = 1:N_T
   V_t_tmp = V_t(1:end,t)'; % current vector   
   imagesc(topographic_brain_map_64(V_t_tmp),[0,256])
   colormap('gray')
   title(['EEG in real time, t = ',num2str(1/N_T*t*10), ' s'])
   xlabel(['Current microstate: ',num2str(L_t(t))])
    drawnow                
end