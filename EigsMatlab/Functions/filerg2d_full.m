function [E0,E_dif,Nu,Nd,N_fill]=filerg2d_full(N_size,t_up,t_down,epsilon_up,epsilon_down,U,V0,r)
N_sites=N_size(1)*N_size(2);
for k=0:2*N_sites
N_fill(k+1)=k;
[E0(k+1),Nu(k+1),Nd(k+1)]=filgs2d_full(N_size,N_fill(k+1),t_up,t_down,epsilon_up,epsilon_down,U,V0,r);
end
E_dif(1)=E0(2);
for k=2:2*N_sites
E_dif(k)=E0(k+1)-E0(k);
end
end
%scatter(N_fill,E0)
%title('Ground state energy vs electron number','interpreter','latex')
%xlabel('electron number','interpreter','latex') % x-axis label
%ylabel('Energy','interpreter','latex') % y-axis label
%ax = gca;
%ax.XTick = (0:1:12);
%gtext('$N_{sites}=6,t= 1, U=10$','interpreter','latex')
