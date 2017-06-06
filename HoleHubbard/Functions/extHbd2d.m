function H_tot= extHbd2d(N_size,N_up,N_down,t_up,t_down,epsilon_up,epsilon_down,U,V0,r)

N_row=N_size(1);
N_col=N_size(2);
N_sites=N_row*N_col;

T=epsilon_up';
epsilon_up=T(:);
T=epsilon_down';
epsilon_down=T(:);
T=U';
U=T(:);
th_up=t_up{1};
tv_up=t_up{2};
th_down=t_down{1};
tv_down=t_down{2};

basis_up=uperm(de2bi(2^N_sites - 2^(N_sites-N_up)));
if basis_up==0
    basis_up=zeros(1,N_sites);
end
dec_up=bi2de(basis_up,'left-msb');
row_index_up=[];
column_index_up=[];
hopp_amplitude_up=[];

if N_up>0
    if N_col>1
for j=1:(N_row)
    for k=1:(N_col-1)
        l=(j-1)*N_col+k;
    [~,L1,L2]=intersect((dec_up+2^(N_sites-l)-2^(N_sites-l-1)).*(1-basis_up(:,l)).* basis_up(:,l+1),dec_up);
row_index_up=[row_index_up;L1];
column_index_up=[column_index_up;L2];
hopp_amplitude_up= [hopp_amplitude_up; ones(length(L1),1)*th_up(j,k)];
    end
end
    end
    if N_row>1
for k=1:(N_col)
    for j=1:(N_row-1)
        l1=(j-1)*N_col+k;
        l2=j*N_col+k;
    [~,L1,L2]=intersect((dec_up+2^(N_sites-l1)- 2^(N_sites-l2)).*(1-basis_up(:,l1)*(1-eq(l1,l2))).* basis_up(:,l2),dec_up);
    exponent = (-1) .^(sum(basis_up(L2,l1:l2),2)-1);
row_index_up=[row_index_up;L1];
column_index_up=[column_index_up;L2];
hopp_amplitude_up= [hopp_amplitude_up; ones(length(L1),1)*tv_up(j,k).* exponent];
    end
end
    end
end

%similarly for spin down
basis_down=uperm(de2bi(2^N_sites - 2^(N_sites-N_down)));
if basis_down==0
    basis_down=zeros(1,N_sites);
end

dec_down=bi2de(basis_down,'left-msb');
row_index_down=[];
column_index_down=[];
hopp_amplitude_down=[];

if N_down>0
    if N_col>1
for j=1:(N_row)
    for k=1:(N_col-1)
        l=(j-1)*N_col+k;
    [~,L1,L2]=intersect((dec_down+2^(N_sites-l)-2^(N_sites-l-1)).*(1-basis_down(:,l)).* basis_down(:,l+1),dec_down);
row_index_down=[row_index_down;L1];
column_index_down=[column_index_down;L2];
hopp_amplitude_down= [hopp_amplitude_down; ones(length(L1),1)*th_down(j,k)];
    end
end
    end
    if N_row>1
for k=1:(N_col)
    for j=1:(N_row-1)
        l1=(j-1)*N_col+k;
        l2=j*N_col+k;
    [~,L1,L2]=intersect((dec_down+2^(N_sites-l1)- 2^(N_sites-l2)).*(1-basis_down(:,l1)*(1-eq(l1,l2))).* basis_down(:,l2),dec_down);
    exponent = (-1) .^(sum(basis_down(L2,l1:l2),2)-1);
row_index_down=[row_index_down;L1];
column_index_down=[column_index_down;L2];
hopp_amplitude_down= [hopp_amplitude_down; ones(length(L1),1)*tv_down(j,k).* exponent];
    end
end
    end
end

%Generate the kinetic part of the Hamiltonian
size_up=nchoosek(N_sites,N_up);
size_down=nchoosek(N_sites,N_down);
H_up=sparse([row_index_up; column_index_up], [column_index_up; row_index_up],[hopp_amplitude_up; hopp_amplitude_up],size_up,size_up);
H_down=sparse([row_index_down; column_index_down], [column_index_down; row_index_down],[hopp_amplitude_down; hopp_amplitude_down],size_down,size_down);
I_up=speye(size_up);
I_down=speye(size_down);
T1=KronProd({I_down,H_up},[1 2],[size_down size_up],1);
T2=KronProd({H_down,I_up},[1 2],[size_down size_up],1);
H_tot=sparse(T1)+sparse(T2);
nj_up=cell(1,N_sites);
nj_down=cell(1,N_sites);
for j=1:N_sites
    nj_up{j}=spdiags(basis_up(:,j),0,size_up,size_up);
    nj_up_tensor=KronProd({I_down,nj_up{j}},[1 2],[size_down size_up],epsilon_up(j));
    nj_down{j}=spdiags(basis_down(:,j),0,size_down,size_down);
    nj_down_tensor=KronProd({nj_down{j},I_up},[1 2],[size_down size_up],epsilon_down(j));
    Uj=KronProd({nj_down{j},nj_up{j}},[1 2],[size_down size_up],U(j));
    H_tot=H_tot+sparse(nj_up_tensor)+sparse(nj_down_tensor)+sparse(Uj);
end
if V0~=0
 for k=1:N_sites-1
     kx=mod(k-1,N_col)+1;
     ky=floor((k-1)/N_col)+1;
     for j=k+1:N_sites
         jx=mod(j-1,N_col)+1;
         jy=floor((j-1)/N_col)+1;
                W=V0/norm([r{1}(ky,kx);r{2}(ky,kx)]-[r{1}(jy,jx);r{2}(jy,jx)]);
                W1=KronProd({nj_down{k},nj_up{j}},[1 2],[size_down size_up],W);
                 W2=KronProd({nj_down{j},nj_up{k}},[1 2],[size_down size_up],W);
                  W3=KronProd({I_down,nj_up{j}*nj_up{k}},[1 2],[size_down size_up],W);
                   W4=KronProd({nj_down{j}*nj_down{k},I_up},[1 2],[size_down size_up],W);
    H_tot=H_tot+sparse(W1)+sparse(W2)+sparse(W3)+sparse(W4);
     end
 end
end
end