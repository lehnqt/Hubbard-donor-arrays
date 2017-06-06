function density = partdens2d(N_size,N_up,N_down,w,sigma)
N_row=N_size(1);
N_col=N_size(2);
density=zeros(N_row,N_col);
for j=1:N_row
    for k=1:N_col
    density(j,k)=spdm2d(N_size,N_up,N_down,w,j,k,j,k,sigma);
    end
end
end
