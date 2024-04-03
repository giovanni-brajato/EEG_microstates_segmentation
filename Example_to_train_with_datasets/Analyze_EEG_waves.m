%% Read some data
close all
[hdr, EEG_waveforms] = edfread('S001R01.edf') ;
EEG_waveforms = EEG_waveforms(1:end-1,:);

%% Visualize EEG waves
figure(1)

hold on
max_potential_value = max(max(EEG_waveforms(1:end,:)));
min_potential_value = min(min(EEG_waveforms(1:end,:)));
peak2peak_max_amplitude = max_potential_value - min_potential_value;
for i = 1:size(EEG_waveforms,1)-1
    plot(i*peak2peak_max_amplitude + EEG_waveforms(i,:))
end
title('EEG waveforms')
hold off

N_T = size(EEG_waveforms,2); %number of samples
N_s = size(EEG_waveforms,1); %number of electrodes


%% BASIC N-TH microstate algortihm

N_microstates = 4; % change here the model
convergence_criterion = 1e-6; %determines when stop the algorithm
sigma_0 = 0;
% initialization 2b: initialize the labels in random uniform fashion
Labels = ceil(N_microstates*rand(1,N_T));

a_kt = zeros(N_microstates,N_T);

 S_k = cell(N_microstates,1);
 Gamma_k = S_k;
 
 loop = true;
 come_back = false;
 while loop
     while come_back
     for t = 1:N_T
         possible_labels = zeros(N_microstates,1);
         for k = 1:N_microstates
             possible_labels(k) = (EEG_waveforms(1:end,t)'*Gamma_k{k})^2;
         end
         [~,index] = max(possible_labels);
         Labels(t) = index;
         
     end
     come_back = false;
     end
 for k = 1:N_microstates
    S_k{k} =  zeros(N_s,N_s);
    for t = 1:N_T
        % calculate the matrices
        S_k{k} = S_k{k} + (Labels(t) == k)*EEG_waveforms(1:end,t)*EEG_waveforms(1:end,t)'; %upate only if they have the label
    end
    % diagonalize S_k - hopefully is not singluar
    [V, D] = eig(S_k{k});
    % find the largest eigenvalue
    [~,max_v_index] = max(max(D));
    % find the correspondent eigenvector
    Gamma_k{k} = V(:,max_v_index) ;
    % compute the noise estimation variance
    

end
sigma_mi = 0;
for t = 1:N_T
    sigma_mi = sigma_mi + (EEG_waveforms(1:end,t)'*EEG_waveforms(1:end,t) - (Gamma_k{Labels(t)}'*EEG_waveforms(1:end,t))^2);
end
sigma_mi = (1/N_T)*(1/(N_s-1))*sigma_mi;

if abs(sigma_mi - sigma_0) < convergence_criterion*sigma_mi
    loop = false;
else
    sigma_0 = sigma_mi;
    come_back = true;
end
 end
 sigma_D = 0;
 for t = 1:N_T
     curr_label = Labels(t);
     for k = 1:N_microstates
         if k ~= curr_label
             a_kt(k,t) = 0;
         end
     end
     a_kt(curr_label,t) = EEG_waveforms(1:end,t)'*Gamma_k{curr_label};
     sigma_D = sigma_D + EEG_waveforms(1:end,t)'*EEG_waveforms(1:end,t);
 end
 sigma_D = (1/N_T)*(1/(N_s-1))*sigma_D;
 
R_squared = 1 - sigma_mi/sigma_D;

%% show us the microstates
figure(3)
for k = 1:N_microstates
    subplot(1,N_microstates,k)
    microstate =  topographic_brain_map_64(Gamma_k{k}');
    imagesc(microstate, [0,max(max(microstate))]);
    title(['Microstate ',num2str(k)])
    
end

%% check results in real time
figure(2)
for t = 1:N_T
    V_t = EEG_waveforms(1:end,t)'; %  current vector   
   imagesc(topographic_brain_map_64(V_t),[0,256])     
   colormap('gray')
   title(['EEG in real time, t = ',num2str(1/N_T*t), 'ms'])
   xlabel(['Current microstate: ',num2str(Labels(t))])
    drawnow                
end