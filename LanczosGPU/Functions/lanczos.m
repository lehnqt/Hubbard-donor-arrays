function [Eg,wg]=lanczos(H_up,H_down,n_up,n_down,ns_up,ns_down)
ne=1;
m=15;

N_sites=length(n_up);

m=2*ne+m;
wg=cell(1,ne);

D_up=length(H_up);
D_down=length(H_down);
%fix CUDA bug
if D_up==1
    H_up=full(H_up);
    for j=1:N_sites
        n_up{j}=full(n_up{j});
        ns_up{j}=full(ns_up{j});
    end
end

if D_down==1
    H_down=full(H_down);
for j=1:N_sites
        n_down{j}=full(n_down{j});
        ns_down{j}=full(ns_down{j});
end
end

m=min(m,D_up*D_down);

alpha=zeros(1,m);
beta=zeros(1,m);
[wT_up,~]=eigs(gather(H_up),1,'sa');
[wT_down,~]=eigs(gather(H_down),1,'sa');
v=wT_down*wT_up';
v1=gpuArray(v/norm(v(:)));
v0=gpuArray(0);
beta(1)=0;
j=1;

%energy%
while j<m
w=full(H_down*v1)+full(fstimes(v1,H_up));
for l=1:N_sites
w=w+full(fstimes(n_down{l}*v1,ns_up{l}));
w=w+full(fstimes(ns_down{l}*v1,n_up{l}));
end
alpha(j)=gather((w(:))'*v1(:));
w=w-alpha(j)*v1-beta(j)*v0;
beta(j+1)=norm(w(:));
v0=v1;
v1=w/beta(j+1);
j=j+1;
end

w=full(H_down*v1)+full(fstimes(v1,H_up));

for l=1:N_sites
w=w+full(fstimes(n_down{l}*v1,ns_up{l}));
w=w+full(fstimes(ns_down{l}*v1,n_up{l}));
end

alpha(m)=gather((w(:))'*v1(:));
T=diag(alpha,0)+diag(beta(2:m),1)+diag(beta(2:m),-1);
[ug,Eg]=eigs(T,ne,'sa');
Eg=diag(Eg);

if nargout>1 

%wavefunction%
v1=gpuArray(v/norm(v(:)));
for idx=1:ne
wg{idx}=ug(1,idx)*v1;
end


j=1;
while j<m
w=full(H_down*v1)+full(fstimes(v1,H_up));
for l=1:N_sites
w=w+full(fstimes(n_down{l}*v1,ns_up{l}));
w=w+full(fstimes(ns_down{l}*v1,n_up{l}));
end

alpha(j)=gather((w(:))'*v1(:));
w=w-alpha(j)*v1-beta(j)*v0;
beta(j+1)=norm(w(:));
v0=v1;
v1=w/beta(j+1);
j=j+1;

for idx=1:ne
wg{idx}=wg{idx}+ug(j,idx)*v1;
end
end

for idx=1:ne
N=norm(wg{idx}(:));
wg{idx}=gather(wg{idx}(:)/N);
end
end
end