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
  
N_size=[10,10];
N_row=N_size(1);
N_col=N_size(2);
N_sites=N_row*N_col;
N_up=ceil(N_sites/4);
N_down=floor(N_sites/4);
  
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
epsilon_down=epsilon_up;