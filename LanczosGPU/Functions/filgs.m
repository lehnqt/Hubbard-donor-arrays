function [Egs,N_up,N_down] = filgs(N_size,N_fill,t_up,t_down,epsilon_up,epsilon_down,U,V0,r)
N_up=ceil(N_fill/2);
N_down=N_fill-N_up;
[H_up,H_down,n_up,n_down,ns_up,ns_down]= Hbd_OBC(N_size,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r);
 Egs=lanczos(H_up,H_down,n_up,n_down,ns_up,ns_down);
end     