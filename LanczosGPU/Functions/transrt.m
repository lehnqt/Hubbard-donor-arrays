function [E,Ep,Q] = transrt(N_size,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r,spin)  
N_row=N_size(1);
N_col=N_size(2);
[H_up,H_down,n_up,n_down,ns_up,ns_down]= Hbd_OBC(N_size,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r);
     [E,wn]=lanczos(H_up,H_down,n_up,n_down,ns_up,ns_down);
    if spin==1
    if N_up==0
        Q=0;
    else
    [H_up,H_down,n_up,n_down,ns_up,ns_down]= Hbd_OBC(N_size,N_up-1,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r);
     [Ep,wnp]=lanczos(H_up,H_down,n_up,n_down,ns_up,ns_down);
            Ql=0;
            Qr=0;      
            for p=1:N_row
                Ql=Ql+(wnp{1}'*cj2d(N_size,N_up,N_down,wn{1},p,1,1))^2;
                Qr=Qr+(wnp{1}'*cj2d(N_size,N_up,N_down,wn{1},p,N_col,1))^2;
            end
            Q=1/(1/Ql+1/Qr);
    end
    end
    if spin==-1
    if N_up==0
        Q=0;
    else
    [H_up,H_down,n_up,n_down,ns_up,ns_down]= Hbd_OBC(N_size,N_up,N_down-1,t_up,t_down,epsilon_up,epsilon_down,U,V0,r);
    [Ep,wnp]=lanczos(H_up,H_down,n_up,n_down,ns_up,ns_down);
            Ql=0;
            Qr=0;      
            for p=1:N_row
                Ql=Ql+(wnp{1}'*cj2d(N_size,N_up,N_down,wn{1},p,1,1))^2;
                Qr=Qr+(wnp{1}'*cj2d(N_size,N_up,N_down,wn{1},p,N_col,1))^2;
            end
            Q=1/(1/Ql+1/Qr);
    end
    end
end