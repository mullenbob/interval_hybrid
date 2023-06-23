%LINEAR INTERVAL FINTE ELEMENT ANALYSIS OF PLANE TRUSS
%SECONDARY VARIABLES WITH THE SAME ACCURACY OF THE PRIMARY ONES
%PROGRAM BY Prof. RAFI MUHANNA AND Prof.  ROBERT MULLEN
% addpath 'G:\My Drive\Documents\MATLAB\Intlab_V8\Intlab_V8'
%  startintlab();
 intvalinit('DisplayInfsup');
clear
clc
format long
global strainflag;  %this is set to 1 if strain values are needed.


tic
intvalinit('DisplayInfsup');
name="popova20group";
inp= fopen(name+'.inp','r');
out =fopen(name+'.out','w');
strainflag=1;   % if set to one, strains will also be calculated
fid=out;
nintvload=0
nintvgroup=0
[nn,nnA,ne,xA,yA,ndof,ndofA, resxA,resyA,fxA,fyA,BetaxA,BetayA, ...
    A,E,elementA,alfaA,ng,Gamma,betaxA,betayA,Alfa,gammao,gn,gamma1]=readtruss2(name,inp,out);

%     nvars=ne+nintvload+nintvgroup;
% lb=zeros(ne+nintvgroup+nintvload,1);
% ub=lb;
% x0=lb;
angle = 0.;
temp=ifep(out,8,angle,nn,nnA,ne,xA,yA,ndof,ndofA, resxA,resyA,fxA,fyA,BetaxA,BetayA, ...
    A,E,elementA,alfaA,ng,Gamma,betaxA,betayA,Alfa,gammao,gn,gamma1);
inf(temp)
sup(temp)
%  %set initial values for modulus values
% for i1=1:ne
%     x0(i1)=E(i1);
%     lb(i1)=x0(i1)*(1.-alfaA(i1)); 
%     ub(i1)=x0(i1)*(1.+ alfaA(i1));
% end
% if (nintvgroup > 0)
%     for i1=1:nintvgroup
%     i2=i1+nel;
%     x0(i2)=1.;
%     lb(i2)=(1-gamma(i1));
%     ub(i2)=(1.+gamma(i1));
%     end
% end
% if (nintvload > 0)
% for i1=1:nintvload
%     i2=i1+nel+nintvgroup;
%     x0(i2)=1.;
%     lb(i2)=(1-beta(i1));
%     ub(i2)=(1.+beta(i1));
% end
% end
% %old matrix version = ebe x naa x lamba x strain x strain
% %new version is ebe (twodof per element)*naa*lamba*strain  where lamba and
% %strain are ne terms.
% options = optimoptions('particleswarm','FunctionTolerance',1.e-10);
% [x,fval,exitflag,output] =particleswarm(@(x)ifep(out,8,x,nn,nnA,ne,xA,yA,ndof,ndofA, resxA,resyA,fxA,fyA,BetaxA,BetayA, ...
%     A,E,elementA,alfaA,ng,Gamma,betaxA,betayA,Alfa,gammao,gn,gamma1), ...
%     nvars,lb,ub,options) 



function [result]=ifep(out,ivalue,optvar,nn,nnA,ne,xA,yA,ndof,ndofA, resxA,resyA,fxA,fyA,BetaxA,BetayA, ...
    A,E,elementA,alfaA,ng,Gamma,betaxA,betayA,Alfa,gammao,gn,gamma1)
global strainflag;
%--------------------------------------------------------------------------
Del_alpha = infsup(-alfaA,alfaA);

fid=out;
nforce=0;
ForceA=infsup(zeros(2*nnA,1),zeros(2*nnA,1)); %allocate space for ForceA  RLM 9/8/22
for i=1:nnA
  if(abs(fxA(i))>0.0)
      nforce = nforce+1;
  end
  if(abs(fyA(i))>0.0)
      nforce = nforce+1;
  end
  ForceA(2*i-1) = BetaxA(i)*fxA(i);
  ForceA(2*i)   = BetayA(i)*fyA(i);
end
nE = 0; %number of elements with uncertain E
for i=1:ne
    if(alfaA(i)>0.0)
        nE = nE+1;
    end
end
% for e=1:ne  no longer used
%     element(e,1)=2*e-1;
%     element(e,2)=2*e;
% end
% fx = zeros(nn,1);
% fy = zeros(nn,1);
% loadflag = zeros(nnA,1);
% for i=1:nnA
%     loadflag(i) = 0;
%     if(fxA(i)~=0||fyA(i)~=0)
%         loadflag(i) =1;
%     end
% end
% % for e=1:ne
%    i1 = elementA(e,1); %nodes for assembled model
%    i2 = elementA(e,2); 
%    j1 = 2*e-1;%nodes for EBE model
%    j2 = 2*e;
%    % %construct map between ebe nodes and assembled model nodes
%    % ebenode(j1)=i1;
%    % ebenode(j2)=i2;  no longer used
%    % 
%    %restraints for EBE model
%    resxEBE(2*e-1) = resxA(i1);
%    resyEBE(2*e-1) = resyA(i1);
%    resxEBE(2*e)   = resxA(i2);
%    resyEBE(2*e)   = resyA(i2);
%    %forces for EBE model
%    % need this code to not replace loading 
%    % if(loadflag(i1)==1)
%        fx(j1)=fxA(i1);
%        fy(j1)=fyA(i1);
%        loadflag(i1)=0;
% end
%    if(loadflag(i2)==1)
%       fx(j2)=fxA(i2);
%       fy(j2)=fyA(i2); 
%       loadflag(i2)=0;
%    end

% betaxEBE = zeros(nn,1);  no longer used
% betayEBE = zeros(nn,1);
% %load uncertainty for EBE model
% for i=1:nn
%   j=ebenode(i);
%   betaxEBE(i) = betaxA(j);
%   betayEBE(i) = betayA(j);
% end

%nlamda = nf*i2-nzeros;  %Complex form   there are always two constraints for each element when it is a truss
nlamda = ne*2 ;  
ndof=2*ne;
nf=nnA;
ndof1=2*nnA;
%nodfnew is the size of the combined stiffness matrix
ndofnew = ndof+ndof1+nlamda;
if (strainflag == 1) 
    ndofnew=ndofnew+2*ne;
end
fprintf(1,'ndof = %d ndof1=%d nlamda= %d ndofnew=%d, strainflag = %d\n',ndof,ndof1,nlamda,ndofnew,strainflag);
%a = zeros(ne,ndofnew);
a = zeros (ne,2*ne);
aA = zeros(ne,ndofA);
Le=zeros(ne,1);
cosine=zeros(ne,1);
sine=zeros(ne,1);
Doo=zeros(ne,ne);
%allocated space for D arrays
Do=Doo*infsup(1.,1.);
D=Do;
DoD=Do;
Dom=Doo;
 for e=1:ne
  connA=elementA(e,:);
  x1 = xA(connA(1));
  y1 = yA(connA(1));
  x2 = xA(connA(2));
  y2 = yA(connA(2));
  le=sqrt((x2-x1)^2+(y2-y1)^2);
  Le(e) = le;
  c=(x2-x1)/le;
  s=(y2-y1)/le;
  cosine(e) = c;
  sine(e) = s;
  %dof=[2*conn(1)-1 2*conn(1) 2*conn(2)-1 2*conn(2)] ;
  dof = [2*e-1 2*e];
 
  dofA=[2*connA(1)-1 2*connA(1) 2*connA(2)-1 2*connA(2)];
   %____________global transformation matrix_________________
 % % R=[ c -s  0  0;
 %      s  c  0  0;
 %      0  0  c -s;
 %      0  0  s  c];
  %RLM CHANGE TO TWO DOF FOR EACH ELEMENT NOT FOUR
  a(e,dof)=[1 -1];
 % a(e,dof)= eye(4)*[1 0 -1  0]'; % RLM I do not see why we need the identy matrix
  aA(e,dofA)= [c s -c -s]' ; %R*[1 0 -1  0]'; % again   why not just vector [c s -c -c]'
  %____________diagonal member matrix_______________________

 Doo(e,e)=A(e)*E(e)/le; 
 Do(e,e)=A(e)*E(e)/le*infsup(1.0,1.0);
 D(e,e)=(A(e)*E(e)/le*Alfa(e));
 DoD(e,e)=Do(e,e)-D(e,e) ; 
 Dom(e,e) = (E(e)*A(e)/le)*(1.+alfaA(e)*gammao(e));
end
aa=zeros(ne,4*ne+2*nf+strainflag*2*ne);
aa(:,1:2*ne)=a;
%----------------formulation of the new approach-------------------------

% % change to just a single force matrix with values starting after the EBE elements
% force=zeros(ndofnew,1)*infsup(0.,0.);
% j=2*ne;
% for i=1:nf
%     j=j+1;
%     force(j)  = infsup(1-betaxA(i),1+betaxA(i))*fx(i);
%      j=j+1;
%      force(j)    = infsup(1-betayA(i),1+betayA(i))*fy(i);
% 
% end
% one constraint term for each EBE node
i2=1;
cons1=zeros(2*ne+ndof1,2*ne);
for i1=1:ne
    cons1(i2,i2)=1;
    i3=elementA(i1,1)+ne;
    i4=2*i3-1;
    cons1(i4,i2)=-cosine(i1);
    i4=i4+1;
    cons1(i4,i2)=-sine(i1);
    i3=elementA(i1,2)+ne;
    i2=i2+1;
    cons1(i2,i2)=1.;
    i4=2*i3-1;
    cons1(i4,i2)=-cosine(i1);
    i4=i4+1;
    cons1(i4,i2)=-sine(i1);
    i2=i2+1;
end
if (strainflag== 1)
    i2=0;
    cons2=zeros(2*ne,ne);
    for i1=1:ne
    i2=i2+1;
    cons2(i2,i1)=-1/Le(i1);
    i2=i2+1;
    cons2(i2,i1)=1./Le(i1);
    end
istrain=2*ne+2*nf+2*ne+2*ne;
end
inostrain=2*ne+2*nf+2*ne;
%two centered matriese K no group and Km with 1+gamma*delta term
if (strainflag == 1)
    Km=zeros(istrain,istrain);
    Km(1:2*ne,inostrain+1:inostrain+ne)=cons2;
    Km(inostrain+1:inostrain+ne,1:2*ne)=cons2';
    Km(inostrain+1:inostrain+ne,inostrain+ne+1:inostrain+2*ne)=-eye(ne,ne);
    Km(inostrain+ne+1:inostrain+2*ne,inostrain+1:inostrain+ne)=-eye(ne,ne);
    K=Km;
else
    Km=zeros(inostrain,inostrain);
    K=Km;
end    
    Km(1:2*ne,1:2*ne)=a'*Dom*a;
    Km(1:2*ne+2*nf,2*ne+2*nf+1:inostrain)=cons1;
    Km(2*ne+2*nf+1:inostrain,1:2*ne+2*nf)=cons1';
    K(1:2*ne,1:2*ne)=a'*Doo*a;
    K(1:2*ne+2*nf,2*ne+2*nf+1:inostrain)=cons1;
    K(2*ne+2*nf+1:inostrain,1:2*ne+2*nf)=cons1';
bn2 = 0;
for i=1:nf
     if(resxA(i)==0)
           bn2 = bn2+1;
           ifix3(bn2) = 2*ne+ 2*i-1;
     end
     if(resyA(i)==0)
          bn2 = bn2+1;
          ifix3(bn2) = 2*ne+ 2*i;
    end
end
Km(ifix3,:)=zeros(bn2,4*ne+nf*2+strainflag*2*ne);
Km(:,ifix3)=zeros(4*ne+2*nf+strainflag*2*ne,bn2);
Km(ifix3,ifix3)=eye(length(ifix3));
K(ifix3,:)=zeros(bn2,4*ne+nf*2+strainflag*2*ne);
K(:,ifix3)=zeros(4*ne+2*nf+strainflag*2*ne,bn2);
K(ifix3,ifix3)=eye(length(ifix3));
%check for centered solution
Cm=inv(Km);
C=inv(K);
%  create force matrix for assembled nodes
ffm=zeros(4*ne+2*nf+strainflag*2*ne,1);
ff=ffm*infsup(0.,0.);
ff(2*ne+1:2*ne+2*nf,1)=ForceA;
ffm(2*ne+1:2*ne+2*nf,1)=mid(ForceA);
% uu is midpoint solution
uu=C*ffm;

%*********************   Initial Enclosure  *******************
w=ones(ne,1);
w1=w-mag(DoD)*mag(aa*C*aa')*w;
w2=mag(DoD)*mag(aa*C*ff);
% allocate alpha for efficency
alpha=zeros(ne,1);
for i=1:ne
    alpha(i)=(w2(i)/w1(i));
end
alphamax=max(alpha(:,:));           
dd=infsup(-alphamax*w,alphamax*w);  
u1 = C*ff + C*aa'*dd;
%************************************************************
 % for i=1:bn2    
 %   aa(:,ifix2(i))=0; 
 % end 
v=aa*u1;
d=(DoD)*v;
v(:,1)=v;
d(:,1)=d; 
for i=1:10
    v(:,i+1)=intersect(((aa*C)*ff+(aa*C*aa')*d(:,i)),(v(:,i)));
    d(:,i+1)=intersect((DoD*v(:,i+1)),(d(:,i)));
end
u1=(C*ff)+(C*aa')*d(:,11); 
%u1 calculated using first approach-Element uncertainty only
% the above is just to verify non group solutions
% Group Uncertainty lambda matrix 
if (ng > 0)
LambdaG=zeros(ne,ng);

    for i = 1:ne
    j = gn(i);
    if j~=0
        LambdaG(i,j)=1;
    end
end
%LambdaG

Del_gamma = gamma1*(infsup(-1,1));  % centered interval multiplier for groups

else
    %no groups just make it all zeros
LambdaG=zeros(ne,1);
Del_gamma=infsup(0.,0);
end
LambdaE = eye(ne,ne);  %may not need if E is EBE but keep for reading equations
%--------------------------------------------------------------------------
vs = (aa*Cm)*ff;
vg = (aa*C)*ff;

for i = 1:10
    vg = hull(((aa*C)*ff - (((aa*C*aa'*Do)*diag(vg))*LambdaE)*Del_alpha),vg);
    vs = hull(((aa*Cm)*ff-((aa*Cm*aa'*Do)*diag(vs)*LambdaE*Del_alpha)+ ...
        ((aa*Cm*aa'*Do)*diag(vs)*LambdaG*Del_gamma)),vs);
end
us = (Cm*ff-((Cm*aa'*Do)*diag(vs)*LambdaE*Del_alpha)+((Cm*aa'*Do)*diag(vs)*LambdaG*Del_gamma));

ugs = C*ff + C*aa'*Do*diag(vg)*LambdaE*Del_alpha; 
IFEAtime = toc


 [option,uc,Cforce,Cstrain ]=groupcombo (fid,ne,ng,nnA,nforce,E,A,Alfa,Gamma,gamma1,ForceA,xA,yA,elementA,gn,resxA,resyA);
  result=u1(ivalue);



%print results
% uc,Cforce,Cstrain are the end-point combitation results
% uu is midpoint solution
% u1 is no group old intersection resutls
% us is group solution hull
% ugs is no group solution hull
outifep(fid,strainflag,option,nnA,ndof,ndofA,ndofnew,ng,gamma1,gn,ne,elementA,Alfa,E,A,Gamma,uu,u1,us,ugs,uc,Cforce,Cstrain);
return

end