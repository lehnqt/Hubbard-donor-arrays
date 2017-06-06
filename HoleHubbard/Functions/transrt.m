function [E,Ep,Q] = transrt(N_size,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r,spin)  
N_row=N_size(1);
N_col=N_size(2);
H = extHbd2d(N_size,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r);
     [wn,E]=eigs(H,1,'sa');
     if spin==1
    if N_up==0
        Q=0;
    else
    H = extHbd2d(N_size,N_up-1,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r);
     [wnp,Ep]=eigs(H,1,'sa');
            Ql=0;
            Qr=0;      
            for p=1:N_row
                Ql=Ql+(wnp'*cj2d(N_size,N_up,N_down,wn,p,1,1))^2;
                Qr=Qr+(wnp'*cj2d(N_size,N_up,N_down,wn,p,N_col,1))^2;
            end
            Q=1/(1/Ql+1/Qr);
        
    end
     end
      if spin==-1
    if N_down==0
        Q=0;
    else
    H = extHbd2d(N_size,N_up,N_down-1,t_up,t_down,epsilon_up,epsilon_down,U,V0,r);
     [wnp,Ep]=eigs(H,1,'sa');
            Ql=0;
            Qr=0;      
            for p=1:N_row
                Ql=Ql+(wnp'*cj2d(N_size,N_up,N_down,wn,p,1,-1))^2;
                Qr=Qr+(wnp'*cj2d(N_size,N_up,N_down,wn,p,N_col,-1))^2;
            end
            Q=1/(1/Ql+1/Qr);
        
    end
     end
end