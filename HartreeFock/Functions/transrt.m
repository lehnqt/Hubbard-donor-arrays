function [E,Q,convg] = transrt(N_size,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r,spin)  
N_row=N_size(1);
N_col=N_size(2);

[E_hf,w_up,w_down,~,~,convg]=HFAF(N_size,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r);
 E=E_hf;
if spin==1
    if N_up==0
        Q=0;
    else
         Ql=0;
         Qr=0;
          for p=1:N_row
       Ql=Ql+(cj2d(N_size,1,0,w_up(:,N_up),p,1,1))^2;
       Qr=Qr+(cj2d(N_size,1,0,w_up(:,N_up),p,N_col,1))^2;
          end
       Q=1/(1/Ql+1/Qr);
    end
end

if spin==-1
    if N_down==0
        Q=0;
    else
         Ql=0;
         Qr=0;
          for p=1:N_row
       Ql=Ql+(cj2d(N_size,1,0,w_down(:,N_down),p,1,-1))^2;
       Qr=Qr+(cj2d(N_size,1,0,w_down(:,N_down),p,N_col,-1))^2;
          end
       Q=1/(1/Ql+1/Qr);
    end
end
end