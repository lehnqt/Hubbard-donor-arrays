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

hopp=load('hopping.txt');
rc=load('displacement.txt');
hopp=-hopp';

G=zeros(2*N_sites,Nsampl);

parfor sampl_id=1:N_sampl
    sampl_id
G(:,sampl_id)= condctavg_full(rc,N_size,hopp,epsilon_up,epsilon_down,U,V0,r,kT,d_lim,ex_lim); 
end
filename=strcat('condct_dist_',num2str(N_row),'by',num2str(N_col),'_e-1K.mat');
save(filename)