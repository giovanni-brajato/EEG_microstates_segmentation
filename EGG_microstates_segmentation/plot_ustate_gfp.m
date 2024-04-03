function [ ustate_plot ] = plot_ustate_gfp(L_t,GFP_t ,time,N_mu_colors )
%PLOT_USTATE_GFP plot the microstate sequence under the GFP signal
% INPUT
%		-	L_t: label sequence
%   - GFP_t: Global Field Power (standard deviation of the signals)
%		- time: time sequence in ms
% 	- N_mu_colors: colors associated for each ustate
% OUTPUT
%		plot the ustate sequence
%
%
%% how many ustates?
N_mu = unique(L_t);

%% plot the results
legend_string = cell(1,length(N_mu));
hold on
for k = 1:length(N_mu)
    area(time(k:k+1),[0,0],0,'FaceColor', N_mu_colors(k,:),'EdgeAlpha',0);
    legend_string{k} = ['\mu-state ', num2str(k)];
end
legend(legend_string)
t = 1;
while t < length(L_t)
     curr_u_state = L_t(t);
     curr_occourrence = find(curr_u_state == [L_t(t:end),0,curr_u_state]);
     curr_duration = find((diff(curr_occourrence) ~= 1),1) ;
     area(time(t:t + curr_duration -1),GFP_t(t:t + curr_duration -1),0,'FaceColor', N_mu_colors(curr_u_state,:),'EdgeAlpha',0);
     t = t + curr_duration;
end
plot(time,GFP_t,'Color',[0,0,0],'LineWidth',1);
hold off
xlabel('Time [ms]')
ylabel('GFP')


