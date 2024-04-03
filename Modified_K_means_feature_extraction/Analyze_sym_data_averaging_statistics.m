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
            
            
            
            
            
            %% Optimize number of microstate
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
            
            % recompute for optimal number of microstate
            [L_t,R2,R2_s,Gamma_k,sigma2_mu ] = modified_Kmean(V_t,bestK,b,lambda);
            N_mu = bestK;
            
            %% extract statistics
            disp('Extracting statistics...')
            % Average duration aka lifespan - average length of time when a mu-state
            % remains stable whenever it appears
            N_apperance_mu_states = zeros(N_mu,1);
            lifespan_mu_states = zeros(N_mu,1);
            prev_mu_state = L_t(1);
            N_apperance_mu_states(prev_mu_state) = 1;
            lifespan_mu_states(prev_mu_state) = 1;
            for t = 1:N_T-1
                next_mu_state = L_t(t+1);
                if next_mu_state == prev_mu_state
                    lifespan_mu_states(next_mu_state) = lifespan_mu_states(next_mu_state) +1;
                else
                    N_apperance_mu_states(next_mu_state) = N_apperance_mu_states(next_mu_state) +1;
                end
                prev_mu_state = next_mu_state;
            end
            average_lifespan_mu_states = lifespan_mu_states./N_apperance_mu_states; % this value is expressed in samples
            
            % frequency - average number of time per second that a mu_state becomes
            % dominant - we need the proper time scale
            frequencies_mu_states = N_apperance_mu_states;
            
            %coverage - fraction of the total recording when a mu_state is dominant
            coverage = lifespan_mu_states;
            
            %topographical shape - representation of a topographical map
            topographic_mu_states = cell(N_mu,1);
            for i = 1:N_mu
                topographic_mu_states{i} = topographic_brain_map_64(Gamma_k{i}');
            end
            
            %amplitude - average GFP during mu-state dominance
            % calculate GFP
            GFP_t = zeros(1,N_T);
            amplitude_mu_states = zeros(N_mu,1);
            for t = 1:N_T
                GFP_t(t) = sqrt(sum((V_t(:,t) - mean(V_t(:,t))).^2));
                amplitude_mu_states(L_t(t)) = amplitude_mu_states(L_t(t)) +  GFP_t(t);
            end
            amplitude_mu_states = amplitude_mu_states./lifespan_mu_states;
            
            % GEV - global explained variance - % of the total variance explained by a
            % u-state -> SEE the references, I have no idea how to calcuate it
            
            %Transition probabilities - we model a markov matrix, transition matrix,
            %and we normalize the rows. The element (i,j) represent the transition
            %probability from u-state i to u-state j. It would be nice to also define a
            %transtion matrix also when building the model and compare both
            Transition_matrix = zeros(N_mu,N_mu);
            for t = 1:N_T-1
                Transition_matrix(L_t(t),L_t(t+1)) = Transition_matrix(L_t(t),L_t(t+1)) +1;
            end
            for k = 1:N_mu % normalization
                Transition_matrix(k,:) = Transition_matrix(k,:)./sum(Transition_matrix(k,:));
            end
            
            
            
            
            
            
            
            
        end
        
        
        sigma2_mcv_curve{b + 1,1} = sigma2_mcv_matrix;
        
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