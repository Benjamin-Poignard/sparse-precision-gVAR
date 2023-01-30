function P=Pmatrix(S,gamma)

N=length(S(1,:));
for i=1:N
    for j=1:i
        if j<i
            if abs(S(i,j))>10^(-9)
                P(i,j)=abs(S(i,j))^(-gamma);
            else P(i,j)=10^9;
            end;
        else P(i,j)=0;
        end;
        P(j,i)=P(i,j);
    end;
end;