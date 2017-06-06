function [E_hf,w_up,w_down,nu,nd,convg]=HFAF(N_size,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r)
convg=1;
itr_num=500;
norm_tol=10^(-3);
norm_diff=1;
N_row=N_size(1);
N_col=N_size(2);
N_sites=N_row*N_col;

nu=ones(N_row,N_col)*N_up/N_sites;
nd=ones(N_row,N_col)*N_down/N_sites;

m=0.2;
for j=1:N_row
for k=1:N_col
nu(j,k)=nu(j,k)+(-1)^(j+k)*m;
nd(j,k)=nd(j,k)-(-1)^(j+k)*m;
end
end
nu=nu*N_up/(sum(nu(:)));
nd=nd*N_down/(sum(nd(:)));

for itr_id=1:itr_num
 itr_id
[H_up,H_down]=Hbd_obc_HF(N_size,nu,nd,t_up,t_down,epsilon_up,epsilon_down,U,V0,r);
[w_up,E_up]=eig(full(H_up));
[w_down,E_down]=eig(full(H_down));
nu2=zeros(N_row,N_col);
for j=1:N_up
    density=partdens2d(N_size,1,0,w_up(:,j),1);
    nu2=nu2+density;
end
nd2=zeros(N_row,N_col);
for j=1:N_down
    density=partdens2d(N_size,0,1,w_down(:,j),-1);
    nd2=nd2+density;
end
norm_diff_p=norm_diff;
norm_diff=(norm(nu2(:)-nu(:))+norm(nd2(:)-nd(:)))
%break when desired accuracy reached
if norm_diff<norm_tol
    break
end

 alpha=0.2*rand()+0.1;
 
nu=(1-alpha)*nu+alpha*nu2;
nd=(1-alpha)*nd+alpha*nd2;
end
itr_id
norm_diff
if norm_diff>norm_tol
    fprintf('did not converge\n');
    convg=0;
end
E_hf=sum(diag(E_up(1:N_up,1:N_up)))+sum(diag(E_down(1:N_down,1:N_down)));
end