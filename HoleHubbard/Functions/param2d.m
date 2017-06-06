t_up=cell(1,2);

t_up{1}=-hopp_amph*ones(N_row,N_col-1);
t_up{2}=-hopp_ampv*ones(N_row-1, N_col);
t_down=t_up;

U=EC*ones(N_row,N_col);

epsilon_up=-EB*ones(N_row,N_col);
epsilon_down=epsilon_up;