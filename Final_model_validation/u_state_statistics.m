function [ average_lifespan,frequency,coverage,amplitude,transition_matrix,GEV ] = u_state_statistics( u_state_sequence,sampling_time,gfp_t )
% FUNCTION: extract statistics from u_states
% INPUT
% - u_state_sequence: sequence of u_states, definied by integer numbers
% from 1 to N_mu
% - sampling_time: sampling tim, needed to return correct values
% - GFP: global field power
% OUTPUT
% - average_lifespan: average length of time when a mu-state
% remains stable whenever it appears
% - frequency:  average number of time per second that a mu_state becomes
% dominant
% - coverage:  fraction of the total recording when a mu_state is dominant
% - amplitude: average GFP during mu-state dominance
% - transition_matrix: contains the transition probabilities from one u_state to another. We model a markov matrix, transition matrix,
%and we normalize the rows. The element (i,j) represent the transition
%probability from u-state i to u-state j. 
% - GEV: % of the total variance explained by a
% u-state
%% INITIALIZATION
N_T = length(u_state_sequence); %total samples
N_u = max(u_state_sequence); 
N_apperance_u_states = zeros(N_u,1); %number of times a ustate occur
lifespan_u_states = zeros(N_u,1);
amplitude = zeros(N_u,1);
transition_matrix = zeros(N_u,N_u);
%% LIFESPAN
prev_u_state = u_state_sequence(1);
N_apperance_u_states(prev_u_state) = 1;
lifespan_u_states(prev_u_state) = 1;
amplitude(prev_u_state) = amplitude(prev_u_state) +  gfp_t(1);
for t = 1:N_T-1
   next_u_state = u_state_sequence(t+1);
   amplitude(next_u_state) = amplitude(next_u_state) + gfp_t(t+1);
   transition_matrix(prev_u_state,next_u_state) = transition_matrix(prev_u_state,next_u_state) +1;
   if next_u_state == prev_u_state
       lifespan_u_states(next_u_state) = lifespan_u_states(next_u_state) +1;
   else
      N_apperance_u_states(next_u_state) = N_apperance_u_states(next_u_state) +1;
   end
   prev_u_state = next_u_state;
end
average_lifespan = lifespan_u_states./N_apperance_u_states;
%% FREQUENCY
samples_per_second = 1/sampling_time;
frequency = N_apperance_u_states/samples_per_second;
%% COVERAGE
coverage = lifespan_u_states;
%% AMPLITUDE
amplitude = amplitude./lifespan_u_states;
%% TRANSITION MATRIX
for k = 1:N_u % normalization
    transition_matrix(k,:) = transition_matrix(k,:)./sum(transition_matrix(k,:));
end
%% GEV - to be implemented soon...
GEV = NaN;

end

