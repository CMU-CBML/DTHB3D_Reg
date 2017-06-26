function ders = DerBasisFun(p,m,U,i,u,n)
%--------------------------------------------------------------
%function Nip = OneBasisFun(p,m,U,i,u)
% NURBS-Book modified (algorithm A2.4)
% evalute basis function at a particular index
%INPUT:
% i          : current knotspan(if inout is from FindSpan function, then keep it as i. Else it is matlab index-1)
% u          : evaluation point
% p          : degree of the basis functions
% U          : knot vector
% m          : (Highest index of the knot span used-1)
% n          : the degree of the derivative 
%OUTPUT:
% ders       : matrix (n ,1 ), the value returned is the basis function
% value evaluated at that index.
%              
%--------------------------------------------------------------

ders = zeros(n+1,1);
ND = zeros(n+2,1);
NN = zeros(p+1,p+1);
saved = 0;
if(u<U(1,i+1) || u>=U(1,i+p+2)),
    for k = 0:n,
        ders(k+1,1) = 0;
    end
    return;
end

for j = 0:p,
    if(u>=U(1,i+j+1) && u<U(1,i+j+2)),
        NN(j+1,1) = 1;
    else
        NN(j+1,1) = 0;
    end
end

for k = 1:p,
    if(NN(1,k)==0)
        saved = 0;
    else
        saved = (u - U(1,i+1))*NN(1,k)/(U(1,i+k+1)-U(1,i+1));
    end
    
    for j = 0:p-k,
        Uleft = U(1,i+j+2);
        Uright = U(1,i+j+k+2);
        if(NN(j+2,k) == 0),
            NN(j+1,k+1) = saved;
            saved = 0;
        else
            temp = NN(j+2,k)/(Uright-Uleft);
            NN(j+1,k+1) = saved+(Uright-u)*temp;
            saved = (u-Uleft)*temp;
        end
    end
end

ders(1,1) = NN(1,p+1);
a00 = 1;
if((U(i+p+1)-U(i+1))~=0),
    
a10 = a00/(U(i+p+1)-U(i+1));
else
    a10 = 0;
end
if((U(i+p+2)-U(i+2))~=0)
    
a11 = -a00/(U(i+p+2)-U(i+2));
else
    a11=0;
end

ders(2,1) = factorial(p)/factorial(p-1)*(a10*NN(1,p)+a11*NN(2,p));

end
        
            
            
            
            

        
        
     
        