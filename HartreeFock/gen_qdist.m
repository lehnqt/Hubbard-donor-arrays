d=4.6;
EB=45.59;
EC=43.86;
V0=123;

% N_size=[4,4];

N_row=N_size(1);
N_col=N_size(2);
N_sites=N_row*N_col;

ff=1/2;
N_up=ceil(N_sites*ff);
N_down=floor(N_sites*ff);

N_sampl=1000; 

   r=cell(1,2);
   r{1}=zeros(N_row,N_col);
   r{2}=zeros(N_row,N_col);
   for j=1:N_row
       for k=1:N_col
           r{1}(j,k)=(k-1)*d;
           r{2}(j,k)=(j-1)*d;
       end
   end
   
U=EC*ones(N_row,N_col);

epsilon_up=-EB*ones(N_row,N_col);
epsilon_down=epsilon_up;

hopp=load('hopping.txt');
rc=load('displacement.txt');
hopp=-hopp';

parfor idx=1:56
    rng('shuffle');
end

parfor sampl_id=1:N_sampl
sampl_id
[Q(sampl_id),convg(sampl_id)]= transrtdist(rc,N_size,N_up,N_down,hopp,epsilon_up,epsilon_down,U,V0,r,1); 
end
fn=strcat('dat_HF_qdist',num2str(N_row),'by',num2str(N_col),'fill',num2str(N_up),'u',num2str(N_down),'d',num2str(N_sampl/1000),'k');
save(fn)