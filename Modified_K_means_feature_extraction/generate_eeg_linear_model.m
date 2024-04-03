function [ V_t,a_kt_true,Gamma_k_true,E_t  ] = generate_eeg_linear_model(b_true, transition_matrix_true,N_s,N_T,sigma2_noise_true,noise_type )
%GENERATE_EEG_LINEAR_MODEL: generate N_s x N_T EEG signals according to a
%linear model  - ref. Roberto D. Pascual-Marqui, Christoph M. Michael,
%Dietrich Lehmann "Segmentation of Brain Electrical Activity into Microstates: Model estimation and Validation"
%	INPUT
%       - b_true: number of samples which we assure a microstate is
%       resting, i.e not changing.
%       - transition_matrix_true: matrix of ustates transition, its rows
%       have to sum up to one
%       - N_s: number of channels
%       - N_T: number of time samples
%       - sigma2_noise_true: variance of the noise
%       - noise_type: 'correlated' or 'uncorrelated'
%   OUTPUT
%       - V_t: EEG signals
%       - a_kt_true: weights of the ustates, according to the
%       non-overlapping property
%       - Gamma_k_true: generated ustates
%       - E_t: vector of noise
%% initialization
N_mu_true = size(transition_matrix_true,1);
a_kt_true = zeros(N_mu_true,N_T); % initialized to zero
L_t_true = zeros(1,N_T); % true labels vectors

%% labels and weights
L_t_true(1) = ceil(rand*N_mu_true); %random starting state
t = 1;
while t <=N_T
    % resting state
    L_t_true(t:t+b_true-1) = repmat(L_t_true(t),1,b_true);
    a_kt_true(L_t_true(t:t+b_true-1),t:t+b_true-1) =2*rand(1,b_true) - 1; % only one is nonzero -> non overlapping property
    t = t + b_true;
    %changing state
    state_select = cumsum(transition_matrix_true(L_t_true(t-1),:));
    L_t_true(t) = min(find(state_select > rand));
end


%% ustates
H = ((eye(N_s) - ones(N_s,N_s)/N_s)); % linear average refence transformation matrix
Gamma_k_true_r = 2*rand(N_mu_true,N_s) -1;
Gamma_k_true = sigma2_noise_true.*(H*Gamma_k_true_r.').'; % build the microstates vectors
for k = 1:N_mu_true %normalization
    Gamma_k_true(k,:) = Gamma_k_true(k,:)./norm(Gamma_k_true(k,:));
end

%% noise
C =2*rand(N_T,N_s) -1;
beta = sigma2_noise_true;
W = zeros(N_s,N_s);
if strcmp(noise_type,'uncorrelated')
    W = eye(N_s,N_s); %unitary matrix for uncorrelated noise
    E_t = beta*H*W*C.';
elseif strcmp(noise_type,'correlated')
    for i = 1:N_s
        for j = 1:N_s
            if abs(i-j) <= 1 || abs(i-j) == (N_s-1)
                W(i,j) = 1/3;
            else
                W(i,j) = 0;
            end
        end
    end
    E_t =  beta*H*W*W.'*H*W*C.';
end
%% EEG signals
V_t = (a_kt_true(L_t_true,:).'*Gamma_k_true(L_t_true,:) +E_t.').';
end

