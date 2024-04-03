clear all
close all

% b loop
N_s = 64;
b =6;
N_T = 1e2*b;
sigma2_noise_true = 0.1;
noise_type = 'uncorrelated';
transition_matrix_true = [0,1,0,0;
                                         0,0,1,0;
                                         0,0,0,1;
                                         1,0,0,0];
%for b = 1:10
%% model generation
    disp(['Generating model for ',num2str(b),' ...'])
     [ V_t,a_kt_true,L_t_true  ] = generate_eeg_linear_model(transition_matrix_true,N_s,N_T,sigma2_noise_true);
     % Optimize number of microstate
disp('Optimize number of microstate...')
K = 1:6;
N_simulation = 10*3;
for N_sim = 1:N_simulation
    for i=1:length(K)
        N_mu = K(i);
        [~,~,~,~,sigma2_mu(N_sim,i) ] = modified_Kmean(V_t,N_mu);
        sigma2_mcv(i) = sigma2_mu(N_sim,i)*((N_s-1)^-1 * (N_s-1-N_mu))^-2;
    end
    sigma2_mcv_mat(N_sim,:) = sigma2_mcv;
    [~,bestK_i] = min(sigma2_mcv);
    bestK(N_sim) = K(bestK_i);
end;
%% results of optimization
figure(1)
errorbar(K,mean(sigma2_mcv_mat,1),std(sigma2_mcv_mat,1),'-o')
title('Predictive residual variance')
ylabel('\sigma^2_{mcv}');
xlabel('N_{\mu}')
figure(2)
histogram(bestK,'Normalization','PDF')


bestK_val = mode(bestK);
misclass_err = 1-sum(bestK==bestK_val)/N_simulation;
%% recompute for optimal number of microstate
[L_t,R2,R2_s,Gamma_k,sigma2_mu ] = modified_Kmean(V_t,bestK_val);


%% statistics
sampling_time = 1;
GFP_t = GFP(V_t) ;
%true statistic
[ average_lifespan_true,frequency_true,coverage_true,amplitude_true,transition_matrix_true,GEV_true ] = u_state_statistics( L_t_true,sampling_time,GFP_t);
%calculated statistic
[ average_lifespan,frequency,coverage,amplitude,transition_matrix,GEV ] = u_state_statistics( L_t,sampling_time,GFP_t );
     

%% check the model

figure(7)
[UU,VV,SS] = svd(V_t);
plot3(SS(1,:),SS(2,:),SS(3,:),'.');
grid on
%end

%% VISUALIZE USTATE
