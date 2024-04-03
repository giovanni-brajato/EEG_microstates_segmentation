pca_coeff = pca(V_t','Centered','off','NumComponents',2);
V_t_pca = pca_coeff'*V_t;

Gamma_k_mat = [Gamma_k{1},Gamma_k{2},Gamma_k{3},Gamma_k{4},Gamma_k{5}];
Gamma_k_mat_pca = pca_coeff'*Gamma_k_mat;

plot(V_t_pca(1,:),V_t_pca(2,:),'.','MarkerSize',0.1);
figure
plot(Gamma_k_mat_pca(1,:),Gamma_k_mat_pca(2,:),'.','MarkerSize',10);
rectangle('Position',[-1 -1 2 2],'Curvature',[1 1],'LineStyle','--');
axis equal
xlim([-1.2 1.2])
ylim([-1.2 1.2])
grid on