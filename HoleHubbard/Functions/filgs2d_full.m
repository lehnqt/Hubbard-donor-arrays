function [Egs,N_up,N_down] = filgs2d_full(N_size,N_fill,t_up,t_down,epsilon_up,epsilon_down,U,V0,r)
    N_sites=N_size(1)*N_size(2);
    lower_n=ceil(N_fill/2);
    upper_n=min(N_sites,N_fill);
    for j=lower_n:upper_n
        N_up=j;
        N_down=N_fill-N_up;
      H0=extHbd2d(N_size,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r);
      E0=eigs(H0,1,'sa');
        erg(j-lower_n+1)=E0;
    end
[Egs,id]=min(erg(:));
N_up=id-1+lower_n;
N_down=N_fill-N_up;
end     