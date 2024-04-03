function [ gfp ] = GFP( V_t )
%GFP calculate global field power
N_T = size(V_t,2);
gfp = zeros(1,N_T);
for t = 1:N_T
   gfp(t) = sqrt(sum((V_t(:,t) - mean(V_t(:,t))).^2));
end

end

