d=4.6;
hopp_amph=7.55;
hopp_ampv=7.55;
EC=43.86;
% V0=0;
EB=45.59;
  
% N_size=[1,6];
N_row=N_size(1);
N_col=N_size(2);
N_sites=N_row*N_col;

% ff=1/2;
N_up=ceil(N_sites*ff);
N_down=floor(N_sites*ff);

% ff=1/2;
% N_up=ceil(N_sites*ff);
% N_down=floor(N_sites*ff);
  
r=cell(1,2);
   r{1}=zeros(N_row,N_col);
   r{2}=zeros(N_row,N_col);
   for j=1:N_row
       for k=1:N_col
           r{1}(j,k)=(k-1)*d;
           r{2}(j,k)=(j-1)*d;
       end
   end
   
t_up=cell(1,2);
t_up{1}=-hopp_amph*(1+0.01*(rand(N_row,N_col-1)-1/2));
t_up{2}=-hopp_ampv*(1+0.01*(rand(N_row-1, N_col)-1/2));
t_down=t_up;

U=EC*ones(N_row,N_col);
epsilon_up=-EB*ones(N_row,N_col);
epsilon_down=epsilon_up;

spin=1; 
[~,~,Q] = transrt(N_size,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r,spin); 
% fn=strcat('dat_eigsml_qorder',num2str(N_row),'by',num2str(N_col));
% save(fn)