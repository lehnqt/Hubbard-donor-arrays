d_avg=4.60751;
EB=45.59;
EC=43.86;
V0=123;
N_size=[2,4];
N_row=N_size(1);
N_col=N_size(2);
N_sites=N_row*N_col;
 
T=0.1;
kT=1.38*(10^(-23))*T*10^3/(1.6*10^(-19));
d_lim=1;
ex_lim=10;
N_sampl=2000; 

   r=cell(1,2);
   r{1}=zeros(N_row,N_col);
   r{2}=zeros(N_row,N_col);
   for j=1:N_row
       for k=1:N_col
           r{1}(j,k)=(k-1)*d_avg;
           r{2}(j,k)=(j-1)*d_avg;
       end
   end
   
U=EC*ones(N_row,N_col);

epsilon_up=-EB*ones(N_row,N_col);
epsilon_down=epsilon_up;

hopp=load('hopping.txt');
rc=load('displacement.txt');
hopp=-hopp';

G=zeros(2*N_sites,N_sampl);

parfor sampl_id=1:N_sampl
    sampl_id
G(:,sampl_id)= condctavg_full(rc,N_size,hopp,epsilon_up,epsilon_down,U,V0,r,kT,d_lim,ex_lim); 
end
filename=strcat('condct_dist_',num2str(N_row),'by',num2str(N_col),'_e-1K.mat');
save(filename)
