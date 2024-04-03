%% Simulate some data
clear all %#ok<CLALL>
close all

%% loop for hyperparameters beta and lambda
b_vector = 0:9;
lambda_vector = 10;
sigma2_mcv_curve = cell(length(b_vector)*length(lambda_vector),2);
for b = b_vector
    for lambda = lambda_vector
    
     sigma2_mcv_curve{b +1,2} = ['b: ' ,num2str(b),' ,','\lambda: ',num2str(lambda)]; % the second cell contains the (beta,lambda) pair
    disp(['Testing the pair ',sigma2_mcv_curve{b +1,2}])

%% simulations
b_true = 5; % after how many steps we change ustate?
N_simulations = 10; % number of averaging simulations
K_vector = 2:15; % How many N_mu's to try
sigma2_mcv_matrix = zeros(N_simulations,length(K_vector));
N_mu_true = 5; % true number of microstates
N_s = 64; % number of channels
N_T = 1e2*b_true;  % number of time samples
for N_exp = 1:N_simulations
disp(['Experiment number ',num2str(N_exp)])


% random transition after 5 time steps
a_kt_true = zeros(N_mu_true,N_T); % initialized to zero
L_t_true = []; % true labels vectors
time_index = 1;
prev_label = ceil(rand*N_mu_true);
while time_index <N_T
    next_label = ceil(rand*N_mu_true);
    while next_label == prev_label
    next_label = ceil(rand*N_mu_true); % random transition to one of the microstates (it can be the same - markov matrix is in the form 1'1/N_u)
    end
    L_t_true = [L_t_true, next_label*ones(1,b_true)]; 
    a_kt_true(next_label,time_index:time_index +b_true -1 ) =2*rand(1,b_true) - 1; % only one is nonzero -> non overlapping property
    time_index = time_index + b_true;
    prev_label = next_label;
end

sigma2_noise_true = 0.01; % noise variance for the linear model


H = ((eye(N_s) - ones(N_s,N_s)/N_s)); % linear average refence transformation matrix
Gamma_k_true_r = 2*rand(N_mu_true,N_s) -1;
Gamma_k_true = sigma2_noise_true.*(H*Gamma_k_true_r.').'; % build the microstates vectors

for k = 1:N_mu_true %normalization
    Gamma_k_true(k,:) = Gamma_k_true(k,:)./norm(Gamma_k_true(k,:));
end

C =2*rand(N_T,N_s) -1;
beta = sigma2_noise_true;


W = eye(N_s,N_s); %unitary matrix for uncorrelated noise
E_t = beta*H*W*C.';

% special matrix for spatial correlation
%{
for i = 1:N_s
    for j = 1:N_s
        if abs(i-j) <= 1 | abs(i-j) == (N_s-1)
        W(i,j) = 1/3;
        else
            W(i,j) = 0;
        end
    end
end



% for spatial correlation - in the paperyou can read that now the
% covariance structure of the noise is different
E_t =  beta*H*W*W.'*H*W*C.'; % not sure if it is correct
%}


 V_t = a_kt_true(L_t_true,:).'*Gamma_k_true(L_t_true,:) +E_t.';
 V_t = V_t.';

%% Visualize EEG waves
%{
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
%}
% Optimize number of microstate
disp('Optimize number of microstate...')
for i=1:length(K_vector)
    disp(['Trying to fit ',num2str(i),' microstates'])
    N_mu = K_vector(i);
    [~,~,~,~,sigma2_mu ] = modified_Kmean(V_t,N_mu,b,lambda); 
    sigma2_mcv(i) = sigma2_mu*((N_s-1)^-1 * (N_s-1-N_mu))^-2;
end
[~,bestK] = min(sigma2_mcv);
bestK = K_vector(bestK);
sigma2_mcv_matrix(N_exp,:) = sigma2_mcv;
%{
figure
plot(K_vector,sigma2_mcv,'-o');
title('Predictive residual variance')
ylabel('\sigma^2_{mcv}');
xlabel('N_{\mu}')
%}
% recompute for optimal number of microstate
[L_t,R2,R2_s,Gamma_k,sigma2_mu ] = modified_Kmean(V_t,bestK,b,lambda);
N_mu = bestK;

%% show us the microstates
%{
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
%}
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
%{
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
pause
%}
end

sigma2_mcv_curve{b + 1,1} = sigma2_mcv_matrix;
% plot error bars
%{
 figure(1)
errorbar(mean(sigma2_mcv_matrix,1),std(sigma2_mcv_matrix,1))
title('Predictive residual variance')
ylabel('\sigma^2_{mcv}');
xlabel('N_{\mu}')
%}
    end
end

%% plot some results
figure(9)
hold on
for i = 1:size(sigma2_mcv_curve,1)
    if size(sigma2_mcv_curve{i,1},1)==1
        plot(sigma2_mcv_curve{i,1})
    else
        errorbar(mean(sigma2_mcv_curve{i,1},1),std(sigma2_mcv_curve{i,1},1))
%         plot(mean(sigma2_mcv_curve{i,1},1),'Color',[i*25.6/256,0,0])
    end
end
hold off
title('Predictive residual variance')
ylabel('\sigma^2_{mcv}');
xlabel('N_{\mu}')
 legend(sigma2_mcv_curve{:,2})