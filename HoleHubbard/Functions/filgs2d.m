function [Egs,N_up,N_down] = filgs2d(N_size,N_fill,t_up,t_down,epsilon_up,epsilon_down,U,V0,r)
N_up=ceil(N_fill/2);
N_down=N_fill-N_up;
H=extHbd2d(N_size,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r);
Egs=eigs(H,1,'sa');
end     