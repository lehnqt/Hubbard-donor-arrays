function density = occupation2d(N_size,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r)
H=extHbd2d(N_size,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r);
[w,~]=eigs(H,1,'sa');
density_up=partdens2d(N_size,N_up,N_down,w,1);
density_down=partdens2d(N_size,N_up,N_down,w,-1);
density=density_up+density_down;
end