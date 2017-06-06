function Q = transrtdist(rc,N_size,N_up,N_down,hopp,epsilon_up,epsilon_down,U,V0,r,spin) 
Mref=[1,0,0;0,1,0;0,0,-1];
 Mrot=[0,-1,0;1,0,0;0,0,1];
 M=Mref*Mrot;
N_row=N_size(1);
N_col=N_size(2);
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

fprintf('calculating transition rate\n');
[~,~,Q] = transrt(N_size,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r,spin); 
end