function [ L_t,R2,R2_s,Gamma_k,sigma2_mu ] = modified_Kmean( V_t, K, b_in, lambda_in )
% This function apply modified K-mean algorithm in two steps for microstate
% clustering. The rwo spets consisnt in basic clustering and smoothing.
% INPUT:  V_t: matrix of raw data from EEG.
%              The matrix must have C x T dimension, where C is the number
%              of channels and T is the temporal dimension.
%         K:   number of cluster the algorithm is supposed to look for.
%         b:   Integer number of time instant over which the smoothing algorithm is
%              looking at; [default = 3]
%         lambda: weight decay, penalty; [default = 5]
% OUTPUT: L_t: Label vector. Microstate assigned to each time instant. It
%              will have dimension 1 x T.
%         R2:  goodness of fit before smoothing.
%         R2_s:goodness of fit after smoothing.

switch nargin
    case 3
        b = b_in;
        lambda = 5;
    case 4
        b = b_in;
        lambda = lambda_in;
    otherwise
        b = 3;
        lambda = 5;
end;
if b>0 && mod(b,1)~=0
    error('b must be an natual number');
end;

N_T = size(V_t,2); %number of samples - number of time instants
N_s = size(V_t,1); %number of electrodes


%% BASIC N-TH microstate algortihm
disp('Basic N-microstate algorithm...')
N_mu = K; % change here the model - number of mirostates
eps = 1e-6; %determines when to stop the algorithm
sigma2_0 = 0;
sigma2_mu_vec = sigma2_0;
% initialization 2b: initialize the labels in random uniform fashion
% taken from Table I
L_t = ceil(N_mu*rand(1,N_T));

% initialize
a_kt = zeros(N_mu,N_T);
S_k = cell(N_mu,1);
Gamma_k = S_k;

iter = 0;
loop = true;
come_back = false;
while loop
    iter = iter+1;
    while come_back
        for t = 1:N_T
            possible_labels = zeros(N_mu,1);
            for k = 1:N_mu
                possible_labels(k) = (V_t(1:end,t)'*Gamma_k{k})^2;
            end
            [~,L_t(t)] = max(possible_labels);
        end
        come_back = false;
    end
    for k = 1:N_mu
        S_k{k} =  zeros(N_s,N_s);
        for t = 1:N_T
            % calculate the matrices
            S_k{k} = S_k{k} + (L_t(t) == k)*V_t(1:end,t)*V_t(1:end,t)'; %upate only if they have the label
        end
        % diagonalize S_k - hopefully is not singluar
        [V, D] = eig(S_k{k});
        % find the largest eigenvalue
        [~,max_v_index] = max(max(D));
        % find the correspondent eigenvector
        Gamma_k{k} = V(:,max_v_index);
    end
    % compute the noise estimation variance (nonpredictive variance)
    sigma2_mu = 0;
    for t = 1:N_T
        sigma2_mu = sigma2_mu + (V_t(1:end,t)'*V_t(1:end,t) - (Gamma_k{L_t(t)}'*V_t(1:end,t))^2);
    end
    sigma2_mu = (1/N_T)*(1/(N_s-1))*sigma2_mu;
    sigma2_mu_vec = [sigma2_mu_vec,sigma2_mu]; %#ok<AGROW>
    % stoppig criterion check
    if abs(sigma2_mu - sigma2_0) < eps*sigma2_mu
        loop = false;
    else
        sigma2_0 = sigma2_mu;
        come_back = true;
    end
    %plot(1:iter,sigma2_mu_vec(2:end));title('\sigma^2_{\mu} convergence'); xlabel('iteration');drawnow
end
sigma2_D = 0;
for t = 1:N_T
    kappa = L_t(t); %current label
    for k = 1:N_mu
        if k ~= kappa
            a_kt(k,t) = 0;
        end
    end
    a_kt(kappa,t) = V_t(1:end,t)'*Gamma_k{kappa};
    sigma2_D = sigma2_D + V_t(1:end,t)'*V_t(1:end,t);
end
sigma2_D = (1/N_T)*(1/(N_s-1))*sigma2_D;
% goodness of fit
R2 = 1 - sigma2_mu/sigma2_D;
%L_t_save = L_t;



%% SEGMENTATION SMOOTHING algortihm
disp('Segmentation smoothing algorithm...');
%initialize
sigma2_0 = 0;
%b = 3;          % window size hyper-parameter
%lambda = 5;     % smoothn. penalty hyper-parameter
LAMBDA_t = L_t; % copy for substituition, alg. will update automatically L_t
e = sigma2_mu;
N_b = zeros(N_mu,length(1+b:N_T-b));

loop = true;
come_back = false;
while loop
    while come_back
        LAMBDA_t = L_t;
        come_back = false;
    end
    for t = 1+b:N_T-b
        for k= 1:N_mu
            N_b(k,t) = sum(LAMBDA_t(t-b:t+b) == k);
            possible_labels(k) = ( V_t(:,t)'*V_t(:,t)-(Gamma_k{k}'*V_t(:,t))^2 ) / (2*e*(N_s-1))-lambda*N_b(k,t);
        end;
        [~, LAMBDA_t(t)] = min( possible_labels );
        L_t(t) = LAMBDA_t(t);        
    end;
    sigma2_mu = 0;
    for t = 1:N_T
        sigma2_mu = sigma2_mu + (V_t(1:end,t)'*V_t(1:end,t) - (Gamma_k{L_t(t)}'*V_t(1:end,t))^2);
    end
    sigma2_mu = (1/N_T)*(1/(N_s-1))*sigma2_mu;
    % stoppig criterion check
    if abs(sigma2_mu - sigma2_0) < eps*sigma2_mu
        loop = false;
    else
        sigma2_0 = sigma2_mu;
        come_back = true;
    end    
end;
sigma2_D_seg = 0;
for t = 1:N_T
    kappa = L_t(t); %current label
    for k = 1:N_mu
        if k ~= kappa
            a_kt(k,t) = 0;
        end
    end
    a_kt(kappa,t) = V_t(1:end,t)'*Gamma_k{kappa};
    sigma2_D_seg = sigma2_D_seg + V_t(1:end,t)'*V_t(1:end,t);
end
sigma2_D_seg = (1/N_T)*(1/(N_s-1))*sigma2_D_seg;
% goodness of fit
R2_s = 1 - sigma2_mu/sigma2_D_seg;


end

