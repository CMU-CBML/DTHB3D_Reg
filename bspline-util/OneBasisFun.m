function Nip = OneBasisFun(p,m,U,i,u)
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
%OUTPUT:
% Nip       : matrix (1 ,1 ), the value returned is the basis function
% value evaluated at that index.
%              
%--------------------------------------------------------------
NN = zeros(1,p+1);
if((i==0&&u==U(1,1)) || ((i==m-p-1)&&u==U(1,m+1))),
    Nip = 1;
    return;
end

if(u<U(1,i+1)|| u>U(1,i+p+2)),
    Nip = 0;
    return;
end

for j = 0:p,
    if(u>= U(1,i+j+1)&& u<U(i+j+2)),
        NN(1,j+1) = 1;
    else
        NN(1,j+1) = 0;
    end
end   
    for k = 1:p,
        if(NN(1,1)==0),
            saved = 0;
        else
            saved = (u-U(1,i+1))*NN(1,1)/(U(1,i+k+1)-U(1,i+1));
        end
        
        for j=0:p-k,
            Uleft = U(1,i+j+2);
            Uright = U(1,i+j+k+2);
            if(NN(1,j+2)==0),
                NN(1,j+1) = saved;
                saved=0;
            else
                temp = NN(1,j+2)/(Uright-Uleft);
                NN(1,j+1) = saved+(Uright-u)*temp;
                saved = (u-Uleft)*temp;
            end
        end
    end
    Nip = NN(1,1);
end

             
            
    
