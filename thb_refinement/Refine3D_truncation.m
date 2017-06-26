function [Ifinal,Efinal, Pfinal] = Refine3D_truncation(index,level,Dem,Elem,Pm,knotvectorU,knotvectorV,knotvectorW,pU,pV,pW)

%This function performs refinement of a B-spline 
%INPUT:
%index: B-psline global index at a particular refinement level
%level: Refinement level of the B-spline
%Dem: B-spline array structure
%Elem: Element array structure
%Pm: control point array structure
%knotvectorU: knot vector in u direction
%knotvectorV: knot vector in v direction
%knotvectorW: knot vector in w direction
%pU: degree of splines in u direction
%pV: degree of splines in u direction
%pW: degree of splines in u direction

%OUTPUT:
%Ifinal: B-spline array structure after refinement
%Efinal: Element array structure after refinement
%Pfinal: Control point array structure after refinement

%control points for current and next level refinement
P1 = Pm{level,1};
P2 = Pm{level+1,1};

%knot vectors for current and next level refinement
knotcu = knotvectorU{level,1};
knotfu = knotvectorU{level+1,1};

knotcv = knotvectorV{level,1};
knotfv = knotvectorV{level+1,1};

knotcw = knotvectorW{level,1};
knotfw = knotvectorW{level+1,1};

%the global index of a bspline
Bindex = index;

%supp_b is the support cell indices of bspline with global index Bindex
supp_b = Dem{level,5}(Bindex,:);

b_ct = 0;
Btemp = zeros(1,1);

%number of control points at level l and l+1
nobcu = size(knotcu,2) - pU - 1;
nobfu = size(knotfu,2) - pU - 1;

nobcv = size(knotcv,2) - pV - 1;
nobfv = size(knotfv,2) - pV - 1;

nobcw = size(knotcw,2) - pW - 1;
nobfw = size(knotfw,2) - pW - 1;

%============================================================================
refineBS = [];
ring_nbr = 1; %only one-ring neighbourhood of the splines are truncated for degree = 2

% Element refinement
for i=1:size(supp_b,2)%Loop over the support cells of the bspline
    if(supp_b(1,i)~=0 && Elem{level,2}(supp_b(1,i),1)==1 ) %zeros values appending, that's why only consider non-zero
        %&& element in the support is active
        Elem{level,2}(supp_b(1,i),1)=0; % Make the element inactive
        refineBS = [refineBS, Elem{level,3}(supp_b(1,i),:)]; %add the B-splines which are non-zero over the element to the to-be-refined array
        CIcell = Elem{level,4}(supp_b(1,i),:); %children elements of the element
        noci = size(CIcell,2); %number of children elements
        for k = 1:noci  %Loop over their children cells
            ind = CIcell(1,k);
            Elem{level+1,2}(ind,1) = 1;  %Make the children elements active
        end
    end
end

refineBS = unique(refineBS);

%% Bspline refinement
if(isempty(refineBS)==0) %we have non-zero size of the array refineBS
    for i =1:size(refineBS,2) %loop over to-be-refined splines
        supp_b = Dem{level,5}(refineBS(1,i),:); %loop over the elements in the support of the spline
        supp_ct = 0; %support cell count
        inact = 0; %inactive cell count
        for j = 1:size(supp_b,2)
            if(supp_b(1,j)~=0)
                supp_ct = supp_ct+1;
                if(Elem{level,2}(supp_b(1,j),1)==0)
                    inact = inact+1;
                end
            end
        end
        
        if(supp_ct==inact) %if support cell count == inactive cell count, refine the B-spline
            if(Dem{level,2}(refineBS(1,i),1)==1) %only refine the B-splines that are initially active
                Dem{level,2}(refineBS(1,i),1) = 0; %set B-spline as passive
                Dem{level,6}(refineBS(1,i),1) = 1; %flag B-spline as passive and not inactve
                nob_child = size(Dem{level,3},2); %number of children B-splines of the spline
                Ci = Dem{level,3}(refineBS(1,i),:); %children B-spline array of the splines
                for j = 1:nob_child %loop over the children splines
                    if(Ci(1,j)~=0)
                        Dem{level+1,2}(Ci(1,j),1) = 1; %set children B-splines as active
                    end
                end
                
                if(Dem{level,7}(refineBS(1,i),1)==0) %store the inactive B-splines in temp array
                    b_ct = b_ct + 1;
                    Btemp(b_ct,1) = refineBS(1,i);
                end
                
                %computing the global indices of the splines in one-ring
                %neighbourhood for truncation
                idu = Dem{level,1}(refineBS(1,i),1);
                idv = Dem{level,1}(refineBS(1,i),2);
                idw = Dem{level,1}(refineBS(1,i),3);
                
                ofx1 = ring_nbr;
                ofx2 = ring_nbr;
                ofy1 = ring_nbr;
                ofy2 = ring_nbr;
                ofz1 = ring_nbr;
                ofz2 = ring_nbr;
                for j = 1:(pU-1)
                    if(idu == j)
                        ofx1 = j-1;
                    end
                end
                
                for j = nobcu:-1:nobcu-(pU-2)
                    if(idu == j)
                        ofx2 = nobcu-j;
                    end
                end
                
                for j = 1:(pV-1)
                    if(idv == j)
                        ofy1 = j-1;
                    end
                end
                
                for j = nobcv:-1:nobcv-(pV-2)
                    if(idv == j)
                        ofy2 = nobcv-j;
                    end
                end
                
                for j = 1:(pW-1)
                    if(idw == j)
                        ofz1 = j-1;
                    end
                end
                
                for j = nobcw:-1:nobcw-(pW-2)
                    if(idw == j)
                        ofz2 = nobcw-j;
                    end
                end
                
                for t3 = (idw-ofz1):1:(idw+ofz2)
                    for t2 = (idv-ofy1):1:(idv+ofy2)
                        for t1 = (idu-ofx1):1:(idu+ofx2)
                            gg = nobcv*nobcu*(t3-1)+nobcu*(t2-1)+t1;
                            %if truncated spline is active and never set as
                            %truncated
                            if(Dem{level,2}(gg,1)==1 && Dem{level,7}(gg,1)==0)
                                %flag the spline as truncated
                                Dem{level,7}(gg,1) = 1;
                                b_ct = b_ct+1;
                                %add the truncated spline in temp array
                                Btemp(b_ct,1) = gg;
                            end
                        end
                    end
                end
            end
        end
    end
end

%Compute the control point coordinates for the next level
for i = 1:size(Btemp,1)
    if(Btemp(1,1)~=0)
        
        itdu= Dem{level,1}(Btemp(i,1),1);
        itdv= Dem{level,1}(Btemp(i,1),2);
        itdw= Dem{level,1}(Btemp(i,1),3);
        
        newkcu = knotcu(itdu:itdu+pU+1);
        newkcv = knotcv(itdv:itdv+pV+1);
        newkcw = knotcw(itdw:itdw+pW+1);
        
        %find the starting knot index for insertion, u direction
        if(itdu < pU+1)
            usu = itdu;
        else
            usu = FindSpan_mex(nobfu,pU,knotcu(1,itdu),knotfu)+1;
        end
        
        %find the ending point for knot insertion, u direction
        if((nobcu-itdu)<pU)
            ueu= size(knotvectorU{level+1,1},2)-(nobcu-itdu);
        else
            ueu = FindSpan_mex(nobfu,pU,knotcu(1,itdu+pU+1),knotfu)+1;
        end
        
        %find the starting knot index for insertion, v direction
        if(itdv < pV+1)
            usv = itdv;
        else
            usv = FindSpan_mex(nobfv,pV,knotcv(1,itdv),knotfv)+1;
        end
        
        %find the ending point for knot insertion, v direction
        if((nobcv-itdv)<pV)
            uev = size(knotvectorV{level+1,1},2)-(nobcv-itdv);
        else
            uev = FindSpan_mex(nobfv,pV,knotcv(1,itdv+pV+1),knotfv)+1;
        end
        
        %find the starting knot index for insertion, w direction
        if(itdw < pW+1)
            usw = itdw;
        else
            usw = FindSpan_mex(nobfw,pW,knotcw(1,itdw),knotfw)+1;
        end
        
        %find the ending point for knot insertion, w direction
        if((nobcw-itdw)<pW)
            uew = size(knotvectorW{level+1,1},2)-(nobcw-itdw);
        else
            uew = FindSpan_mex(nobfw,pW,knotcw(1,itdw+pW+1),knotfw)+1;
        end
        
        %extract the knot vector containing the finer splines
        newkfu = knotfu(usu:ueu);
        newkfv = knotfv(usv:uev);
        newkfw = knotfw(usw:uew);
        
        unb1 = size(newkcu,2)-pU-1;
        unb2 = size(newkfu,2)-pU-1;
        
        vnb1 = size(newkcv,2)-pV-1;
        vnb2 = size(newkfv,2)-pV-1;
        
        wnb1 = size(newkcw,2)-pW-1;
        wnb2 = size(newkfw,2)-pW-1;
        
        %compute the refinement coefficients according to knot insertion
        %algorithm
        TmatU =  Tmatrix(newkcu,newkfu,pU);
        TmatV =  Tmatrix(newkcv,newkfv,pV);
        TmatW =  Tmatrix(newkcw,newkfw,pW);
        
        TmatU = TmatU(1:unb2,1:unb1);
        TmatV = TmatV(1:vnb2,1:vnb1);
        TmatW = TmatW(1:wnb2,1:wnb1);
        
        %compute the next level control points
        Pnew1 = zeros(unb2,vnb2,wnb2);
        Pnew2 = Pnew1;
        Pnew3 = Pnew1;
        for kk3=1:wnb2
            for kk2 =1:vnb2
                for kk1 = 1:unb2
                    Pnew1(kk1,kk2,kk3) = TmatU(kk1,1)*TmatV(kk2,1)*TmatW(kk3,1)*P1(Btemp(i,1),1);
                    Pnew2(kk1,kk2,kk3) = TmatU(kk1,1)*TmatV(kk2,1)*TmatW(kk3,1)*P1(Btemp(i,1),2);
                    Pnew3(kk1,kk2,kk3) = TmatU(kk1,1)*TmatV(kk2,1)*TmatW(kk3,1)*P1(Btemp(i,1),3);
                end
            end
        end
        
        for kk1=1:unb2
            for kk2 =1:vnb2
                for kk3 = 1:wnb2
                    ing = nobfv*nobfu*(usw+(kk3-1)-1) + nobfu*(usv+(kk2-1)-1)+(usu+(kk1-1));
                    
                    P2(ing,1) = P2(ing,1)+Pnew1(kk1,kk2,kk3);
                    P2(ing,2) = P2(ing,2)+Pnew2(kk1,kk2,kk3);
                    P2(ing,3) = P2(ing,3)+Pnew3(kk1,kk2,kk3);
                    
                end
            end
        end
    end
end

%update the element and B-spline hierarchy with the updated values
Pm{level+1,1} = P2;
Pfinal = Pm;
Ifinal = Dem;
Efinal = Elem;
end





