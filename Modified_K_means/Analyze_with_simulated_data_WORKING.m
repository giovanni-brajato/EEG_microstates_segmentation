%% Simulate some data
clear all %#ok<CLALL>
close all

N_mu_true = 5;
N_s = 64;
N_T = 1e3*N_mu_true;
a_kt_true = 2*rand(N_s,N_mu_true) - 1; % uniform in [-1,1]
V_t = zeros(N_s,N_T);
L_t_true = [];
for k = 1:N_mu_true
  L_t_true = [L_t_true, k*ones(1,N_T/N_mu_true)];
  idx = find(L_t_true == k);
  Gamma_k_true(k,:) = randn(1,length(idx));
  V_t(:,idx) = a_kt_true(:,k)*Gamma_k_true(k,:) + 0.1*randn(N_s,length(idx));
end
% sigma2_noise_true = 0.01;
% 
% 
% H = ((eye(N_s) - ones(N_s,N_s)/N_s));
% Gamma_k_true_r = 2*rand(N_mu_true,N_s) -1;
% Gamma_k_true = sigma2_noise_true.*(H*Gamma_k_true_r.').';
% 
% for k = 1:N_mu_true %normalization
%     Gamma_k_true(k,:) = Gamma_k_true(k,:)./norm(Gamma_k_true(k,:));
% end
% 
% C =2*rand(N_T,N_s) -1;
% beta = sigma2_noise_true;
% W = eye(N_s,N_s);
% 
% E_t = beta*H*W*C.';
% 
% 
% 
%  V_t = a_kt_true(L_t_true,:)*Gamma_k_true(L_t_true,:) +E_t.';
%  V_t = V_t.';

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
K = 3:8;
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
figure(3)
for k = 1:N_mu
    subplot(1,N_mu,k)
    microstate =  topographic_brain_map_64(Gamma_k{k}');
    imagesc(microstate, [0,max(max(microstate))]);
    title(['Microstate ',num2str(k)])
end
figure(4)
for k = 1:N_mu_true
    subplot(1,N_mu_true,k)
    microstate =  topographic_brain_map_64(Gamma_k_true(k,:));
    imagesc(microstate, [0,max(max(microstate))]);
    title(['True Microstate ',num2str(k)])
end
%% check results in real time
%{
figure
for t = 1:N_T
   V_t_tmp = V_t(1:end,t)'; % current vector   
   imagesc(topographic_brain_map_64(V_t_tmp),[0,256])
   colormap('gray')
   title(['EEG in real time, t = ',num2str(1/N_T*t*10), ' s'])
   xlabel(['Current microstate: ',num2str(L_t(t))])
    drawnow                
end
 %}
%% microstate sequence
figure(6)
plot(L_t_true)
hold on
plot(L_t)
hold off
legend('True labels','estimated labels')

figure(7)
[UU,VV,SS] = svd(V_t);
plot3(SS(1,:),SS(2,:),SS(3,:),'.','MarkerSize',0.1);
grid on