 d_avg=4.6;
 hopp_amph=7.55;
 hopp_ampv=7.55;
 Wnr=22.21;
 EC=43.86;
 V0=123;
dr=0;
EB=45.59;
  T=0.1;
kT=1.38*(10^(-23))*T*10^3/(1.6*10^(-19));
d_lim=1;
ex_lim=10; 
  
N_size=[3,3];
N_row=N_size(1);
N_col=N_size(2);
N_sites=N_row*N_col;
N_up=ceil(N_sites/2);
N_down=N_sites-N_up;
  
r=cell(1,2);
   r{1}=zeros(N_row,N_col);
   r{2}=zeros(N_row,N_col);
   for j=1:N_row
       for k=1:N_col
           r{1}(j,k)=(k-1)*d_avg+dr*(rand(1)-1/2);
           r{2}(j,k)=(j-1)*d_avg+dr*(rand(1)-1/2);
       end
   end
   
   t_up=cell(1,2);

t_up{1}=-hopp_amph*ones(N_row,N_col-1);
t_up{2}=-hopp_ampv*ones(N_row-1, N_col);
t_down=t_up;

U=EC*ones(N_row,N_col);

epsilon_up=-EB*ones(N_row,N_col);
for j1=1:N_row
    for k1=1:N_col
    for j2=1:N_row
        for k2=1:N_col
        if (j1-j2)^2+(k1-k2)^2>0
            epsilon_up(j1,k1)=epsilon_up(j1,k1)-V0/norm([r{1}(j1,k1);r{2}(j1,k1)]-[r{1}(j2,k2);r{2}(j2,k2)]);
        end
        end
    end
    end
end
epsilon_down=epsilon_up;