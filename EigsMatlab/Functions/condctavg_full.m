function g_req = condctavg_full(rc,N_size,hopp,epsilon_up,epsilon_down,U,V0,r,kT,d_lim,ex_lim) 
Mref=[1,0,0;0,1,0;0,0,-1];
 Mrot=[0,-1,0;1,0,0;0,0,1];
 M=Mref*Mrot;
N_row=N_size(1);
N_col=N_size(2);
 N_sites=N_row*N_col;
t_up=cell(1,2);
t_up{1}=zeros(N_row,N_col-1);
t_up{2}=zeros(N_row-1,N_col);
disid=randsample(length(rc),N_row*N_col,true);
disid=reshape(disid,[N_row,N_col]);
 for j=1:N_row
     for k=1:N_col-1
        id=(disid(j,k)-1)*length(rc)+disid(j,k+1);
         t_up{1}(j,k)=hopp(id)*(1+0.01*(rand()-1/2));
     end
 end
 
 for j=1:N_row-1
     for k=1:N_col
         r1=rc(disid(j,k),:)*M';
         r2=rc(disid(j+1,k),:)*M';
        [~,id1]=ismember(r1,rc,'rows');
        [~,id2]=ismember(r2,rc,'rows');
        if id1*id2==0
            fprintf('displacement not found');
            break
        end
        id=(id1-1)*length(rc)+id2;
         t_up{2}(j,k)=hopp(id)*(1+0.01*(rand()-1/2));
     end
 end   

t_down=t_up;

Dmax=maxdim2d(N_size,t_up,t_down,epsilon_up,epsilon_down,U,V0,r,kT,d_lim,ex_lim);
[~,E_dif,~]=filerg2d_full(N_size,t_up,t_down,epsilon_up,epsilon_down,U,V0,r);

fprintf('calculating transition rate\n');
for N_up=0:N_sites
    for N_down=0:N_sites
        N_mat(N_up+1,N_down+1)=N_up+N_down;
        [E,Q1,Q2,D] = transrt_L2d(N_size,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r,Dmax);
        se=length(E);
        su=size(Q1);
        sd=size(Q2);
        En(N_up+1,N_down+1,1:se)=E;
        Q_up(N_up+1,N_down+1,1:su(1),1:su(2))=Q1;
        Q_down(N_up+1,N_down+1,1:sd(1),1:sd(2))=Q2;
        Dn(N_up+1,N_down+1)=D;
    end
end

chem=E_dif;
g_req=zeros(length(chem),1);
for k=1:length(chem)
    g_req(k)=reqcondct2d(En,Q_up,Q_down,Dn,N_mat,chem(k),kT,ex_lim);
end
end