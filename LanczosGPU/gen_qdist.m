d=4.6;
EB=45.59;
EC=43.86;
V0=123;
N_size=[1,4];
N_row=N_size(1);
N_col=N_size(2);
N_sites=N_row*N_col;

ff=1/2;
N_up=ceil(N_sites*ff);
N_down=floor(N_sites*ff);
 
 Mref=[1,0,0;0,1,0;0,0,-1];
 Mrot=[0,-1,0;1,0,0;0,0,1];
 M=Mref*Mrot;

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
s=load('separation.txt');
rc=load('displacement.txt');
hopp=-hopp';

spin=1;

spmd
gd=gpuDevice;
idx=gd.Index;
disp(['Using GPU ',num2str(idx)]);
end

parfor idx=1:10
    rng('shuffle');
end

parfor sampl_id=1:N_sampl
sampl_id
Q(sampl_id)= transrtdist(rc,N_size,N_up,N_down,hopp,epsilon_up,epsilon_down,U,V0,r,spin); 
end
fn=strcat('dat_lczgpu_qdist',num2str(N_row),'by',num2str(N_col),'fill',num2str(N_up),'u',num2str(N_down),'d');
save(fn)