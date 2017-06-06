function [En,Q_up,Q_down,Dn] = transrt_L2d(N_size,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r,Dmax)  
N_row=N_size(1);
N_col=N_size(2);
Hn=extHbd2d(N_size,N_up,N_down, t_up,t_down, epsilon_up,epsilon_down,U,V0,r);
    Dn=min(length(Hn),Dmax);
    [wn,En]=eigs(Hn,Dn,'sa');
    En=diag(En);
    if N_up==0
        Q_up=zeros(1,Dn);
    else
    Hnp=extHbd2d(N_size,N_up-1,N_down, t_up,t_down, epsilon_up,epsilon_down,U,V0,r);
    Dnp=min(length(Hnp),Dmax);
    [wnp,~]=eigs(Hnp,Dnp,'sa');
    for j=1:Dnp
        for k=1:Dn
            Ql=0;
            Qr=0;      
            for p=1:N_row
                Ql=Ql+((wnp(:,j))'*cj2d(N_size,N_up,N_down,wn(:,k),p,1,1))^2;
                Qr=Qr+((wnp(:,j))'*cj2d(N_size,N_up,N_down,wn(:,k),p,N_col,1))^2;
            end
            Q_up(j,k)=1/((1/Ql+1/Qr));
        end
    end
    end
    if N_down==0
        Q_down=zeros(1,Dn);
    else
    Hnp=extHbd2d(N_size,N_up,N_down-1, t_up,t_down, epsilon_up,epsilon_down,U,V0,r);
    Dnp=min(length(Hnp),Dmax);
    [wnp,~]=eigs(Hnp,Dnp,'sa');
    for j=1:Dnp
        for k=1:Dn
            Ql=0;
            Qr=0; 
            for p=1:N_row
                Ql=Ql+((wnp(:,j))'*cj2d(N_size,N_up,N_down,wn(:,k),p,1,-1))^2;
                Qr=Qr+((wnp(:,j))'*cj2d(N_size,N_up,N_down,wn(:,k),p,N_col,-1))^2;
            end
               Q_down(j,k)=1/((1/Ql+1/Qr));
        end
    end
    end
end