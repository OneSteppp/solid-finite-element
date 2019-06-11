function [RA,RB,n,X]=liezhu(A,b)
B=[A b]; n=length(b); RA=rank(A); 
RB=rank(B);zhica=RB-RA;
if zhica>0
disp('请注意：因为RA~=RB，所以此方程组无解.')
return
end
if RA==RB
   if RA==n
disp('请注意：因为RA=RB=n，所以此方程组有唯一解.') 
 X=zeros(n,1); C=zeros(1,n+1);
      for p= 1:n-1
[Y,j]=max(abs(B(p:n,p))); C=B(p,:);
B(p,:)= B(j+p-1,:); B(j+p-1,:)=C;
for k=p+1:n
             m= B(k,p)/ B(p,p);
 B(k,p:n+1)= B(k,p:n+1)-m* B(p,p:n+1);
end
end
         b=B(1:n,n+1);A=B(1:n,1:n); X(n)=b(n)/A(n,n);
      for q=n-1:-1:1
         X(q)=(b(q)-sum(A(q,q+1:n)*X(q+1:n)))/A(q,q);
      end
else
         disp('请注意：因为RA=RB<n，所以此方程组有无穷多解.')
   end
end
