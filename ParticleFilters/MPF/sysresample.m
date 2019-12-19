function i=sysresample(q)
qc=cumsum(q);    M=length(q);
u=([0:M-1]+rand(1))/M;
i=zeros(1,M);    k=1;
for j=1:M
  while (qc(k)<u(j))
    k=k+1;
  end
  i(j)=k;
end
end