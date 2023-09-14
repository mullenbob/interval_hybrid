function disp=trusssolve(angle,iii,signx)
method=2;
nskip=20;
% sensitivity code for hybrid solutions   RLM 2023-09-14
global ndof;
global delta;
global nel;
global tx;
global ty;
global force;
global nnd
global nel
global ndof
global telement
global tE;
global tA;
global fxA;
global fyA;
global resxA;
global resyA;
global alphaA;
global beta;
global nintvload;
global nintvgroup;
global groupid;
global gamma;
global beta;
global nodeid;
persistent K;
persistent D1
persistent A1;
persistent R0;
persistent CC;
persistent L;

if (iii == 1)  % build stiffness once
    L=zeros(ndof,nel+nintvgroup+nintvload);
K=zeros(ndof);
D1=zeros(nel,1);
A1=zeros(nel,ndof);
  %assenble stiffness matrix
  for e=1:nel
     i=telement(e,1);
     j=telement(e,2);
     txx=tx(j)-tx(i);
     tyy=ty(j)-ty(i);
     l=sqrt(txx^2+tyy^2);
     c=txx/l;
     s=tyy/l;
     B1= [c s] /l;
     D1(e) = (tE(e)*tA(e)*l);
     k11=(tE(e)*tA(e)/l)*[(c)^2 s*c;s*c s^2];
    K((2*i-1):2*i,(2*i-1):2*i)=K((2*i-1):2*i,(2*i-1):2*i)+k11;
    K((2*j-1):2*j,(2*j-1):2*j)=K((2*j-1):2*j,(2*j-1):2*j)+k11;
    K((2*i-1):2*i,(2*j-1):2*j)=K((2*i-1):2*i,(2*j-1):2*j)-k11;
    K((2*j-1):2*j,(2*i-1):2*i)=K((2*j-1):2*j,(2*i-1):2*i)-k11;
    A1(e,(2*i-1):2*i)=A1(e,(2*i-1):2*i)+B1;
    A1(e,(2*j-1):2*j)=A1(e,(2*j-1):2*j)-B1;
  end
%%
% force fixed dof
for e=1:nnd
    
if (resxA(e) == 0)
        ii=e*2-1;
        A1(:,ii)=0;
        for i=1:nnd
            i1=i*2-1;
            K(ii,i1)=0.;
            K(i1,ii)=0.;
            i1=i1+1;
            K(ii,i1)=0.;
            K(i1,ii)=0.;
        end
        K(ii,ii)=1.0;
end
if (resyA(e) == 0)
        ii=e*2;
        A1(:,ii)=0;
        for i=1:nnd
            i1=i*2-1;
            K(ii,i1)=0.;
            K(i1,ii)=0.;
            i1=i1+1;
            K(ii,i1)=0.;
            K(i1,ii)=0.;
            end
        K(ii,ii)=1.0;
end
end
% force fixed dof
for e=1:nnd
    
if (resxA(e) == 0)
        ii=e*2-1;
        A1(:,ii)=0;
        for i=1:nnd
            K(ii,i)=0.;
            K(i,ii)=0.;
        end
        K(ii,ii)=1.0;
end
if (resyA(e) == 0)
        ii=e*2;
        A1(:,ii)=0;
        for i=1:nnd
            K(ii,i)=0.;
            K(i,ii)=0.;
        end
        K(ii,ii)=1.0;
end
end
CC=inv(K);
end % build stiffness once
% build force vector
ii=0;
for i=1:nnd
ii=ii+1;
force(ii)=fxA(i)*cos(angle)+fyA(i)*sin(angle);
ii=ii+1;
force(ii)=fyA(i)*cos(angle)-fxA(i)*sin(angle);
end


R0=CC*(force');
%temp=R0(2*nnd-1);
% debug step  check A1 and D matrix   calculate applied forces and
% reactions   THe resuts check out   zero terms are 10^(-13)  
% need to add B.C. to A1
% Dint=diag(D1);
% F3=transpose(A1)*Dint*A1*R0;
% KK2=transpose(A1)*Dint*A1;
% for e=1:nnd
% 
% if (resxA(e) == 0)
%         ii=e*2-1;
%         A1(:,ii)=0;
%         for i=1:nnd
%             KK2(ii,i)=0.;
%             KK2(i,ii)=0.;
%         end
%         KK2(ii,ii)=1.0;
% end
% if (resyA(e) == 0)
%         ii=e*2;
%         A1(:,ii)=0;
%         for i=1:nnd
%             KK2(ii,i)=0.;
%             KK2(i,ii)=0.;
%         end
%         KK2(ii,ii)=1.0;
% end
% end
% CC2=inv(KK2);
% R2=CC2*transpose(force);
%calculated sensitivity for each element value
% skip sensitivity direction assuming similar from term to term
if (mod(iii,20) == 1)
v1SAVE=A1*R0;
v0=A1*CC*transpose(force);
u0=CC*transpose(force);   %work with u version of fixed point
for npoint=1:nel
 
%v1=A1*R0;   %THIS CAN BE DONW ONLY ONCE FOR EACH PERTUBATIONS
v1=v1SAVE;   %change to increating only alpha
Dint=D1;
Dint(npoint)=Dint(npoint)*(1+alphaA(npoint));

% look at only doing iteraction in U 
if (method == 2)
    temp1=A1.*(D1-Dint);
   
    temp2=CC*(transpose(temp1)*A1);
    temp3=temp2*R0;
R1=u0+temp3;
else
D0=D1;
DD=(D0-Dint).*v1;

%pull out force term   %this can be done only once for all pertabations
%v0=A1*CC*transpose(force);
for i1=1:1
    
v1=(v0+A1*CC*transpose(A1)*DD);
DD=(D0-Dint).*v1;
end


% split multipication into two parts for tighter bounds
R1=CC*transpose(force)+CC*transpose(A1)*DD;
end
L(:,npoint)=sign(R1-R0);


end
%add in group sensitivity
if (nintvgroup > 0)
    for ng=1:nintvgroup
 
v1=v1SAVE;  % change to saved value
%change to increating only alpha

if (method == 2)
Dint=D1;
for i1=1:nel
    if groupid(i1) == ng
Dint(i1)=Dint(i1)*(1+gamma(ng));
    end
end

    temp1=A1.*(D1-Dint);
   
    temp2=CC*(transpose(temp1)*A1);
    temp3=temp2*R0;
R1=u0+temp3;
else
    Dint=diag(D1);
D0=diag(D1);
for i1=1:nel
    if groupid(i1) == ng
Dint(i1,i1)=Dint(i1,i1)*(1+gamma(ng));
    end
end
DD=(D0-Dint)*v1;


%pull out force term
v0=A1*CC*transpose(force);
for i1=1:1
    
v1=(v0+A1*CC*transpose(A1)*DD);
DD=(D0-Dint)*v1;
end


% split multipication into two parts for tighter bounds
R1=CC*transpose(force)+CC*transpose(A1)*DD;
end


L(:,nel+ng)=sign(R1-R0);
    end
end
end %  supress sensitivity L calculations on each iteration
%add in load sensitivity
 forcetemp=transpose(force);
if (nintvload > 0)
   
    for ng=1:nintvload
       
        i1=nodeid(ng);
        forcetemp(i1)=force(i1)*(1+beta(ng));
  %note added only for force problem for paper
        forcetemp(i1+1)=force(i1+1)*(1+beta(ng));

% split multipication into two parts for tighter bounds
R1=CC*forcetemp;

L(:,nel+nintvgroup+ng)=sign(R1-R0);
    end
end
%end  % end of initial sensistive calculation
%% find sensitivity solution for top right node R0(2*nnd-1);
if (signx == 1)
    Dint=D1;
if nintvgroup == 0
for npoint=1:nel
Dint(npoint)=Dint(npoint)*(1+alphaA(npoint)*L(2*nnd-1,npoint));
end
else
    for npoint=1:nel
        ig=groupid(npoint);
        if (ig == 0 || gamma(ig) == 0)
            Dint(npoint)=Dint(npoint)*(1+alphaA(npoint)*L(2*nnd-1,npoint));
        else
Dint(npoint)=Dint(npoint)*(1+alphaA(npoint)*L(2*nnd-1,npoint))* ...
    (1+gamma(ig)*L(2*nnd-1,nel+ig));
 
    end
    end
end
KK2=transpose(A1.*Dint)*A1;

for e=1:nnd
    
if (resxA(e) == 0)
        ii=e*2-1;
        A1(:,ii)=0;
         for i=1:nnd
            i1=i*2-1;
            KK2(ii,i1)=0.;
            KK2(i1,ii)=0.;
            i1=i1+1;
            KK2(ii,i1)=0.;
            KK2(i1,ii)=0.;
        end
        KK2(ii,ii)=1.0;
        
end
if (resyA(e) == 0)
        ii=e*2;
        A1(:,ii)=0;
        for i=1:nnd
            i1=i*2-1;
            KK2(ii,i1)=0.;
            KK2(i1,ii)=0.;
            i1=i1+1;
            KK2(ii,i1)=0.;
            KK2(i1,ii)=0.;
        end
        KK2(ii,ii)=1.0;
end
end
forcetemp=transpose(force);
for i1=1:nintvload
    i2=nodeid(i1);
    forcetemp(i2)=forcetemp(i2)*(1+beta(i1)*L(2*nnd-1,nel+nintvgroup+i1));
  % added for angle problem
    forcetemp(i2+1)=forcetemp(i2+1)*(1+beta(i1)*L(2*nnd-1,nel+nintvgroup+i1));
end

R2=inv(KK2)*forcetemp;
R3=R2;
else %lower bound
Dint=D1;
if nintvgroup == 0
for npoint=1:nel
Dint(npoint)=Dint(npoint)*(1-alphaA(npoint)*L(2*nnd-1,npoint));
end
else
    for npoint=1:nel
        ig=groupid(npoint);
        if (ig == 0 || gamma(ig) == 0)
            Dint(npoint)=Dint(npoint)*(1-alphaA(npoint)*L(2*nnd-1,npoint));
        else
Dint(npoint)=Dint(npoint)*(1-alphaA(npoint)*L(2*nnd-1,npoint))* ...
    (1-gamma(ig)*L(2*nnd-1,nel+ig));
 
    end
    end
end

KK2=transpose(A1.*Dint)*A1;
for e=1:nnd
    
if (resxA(e) == 0)
        ii=e*2-1;
        A1(:,ii)=0;
         for i=1:nnd
            i1=i*2-1;
            KK2(ii,i1)=0.;
            KK2(i1,ii)=0.;
            i1=i1+1;
            KK2(ii,i1)=0.;
            KK2(i1,ii)=0.;
        end
        KK2(ii,ii)=1.0;
        
end
if (resyA(e) == 0)
        ii=e*2;
        A1(:,ii)=0;
        for i=1:nnd
            i1=i*2-1;
            KK2(ii,i1)=0.;
            KK2(i1,ii)=0.;
            i1=i1+1;
            KK2(ii,i1)=0.;
            KK2(i1,ii)=0.;
        end
        KK2(ii,ii)=1.0;
end
end
forcetemp=transpose(force);
for i1=1:nintvload
    i2=nodeid(i1);
   
    forcetemp(i2)=forcetemp(i2)*(1-beta(i1)*L(2*nnd-1,nel+nintvgroup+i1));
    forcetemp(i2+1)=forcetemp(i2+1)*(1-beta(i1)*L(2*nnd-1,nel+nintvgroup+i1));
   %  added for angle problem
end
R3=inv(KK2)*forcetemp;
R2=R3;
end

disp=[R2(2*nnd-1),R3(2*nnd-1)];
return  
end