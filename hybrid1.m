% hybrid 1   particle swarm -fix-point methods using
%LINEAR INTERVAL FINTE ELEMENT ANALYSIS OF PLANE TRUSS
%SECONDARY VARIABLES WITH THE SAME ACCURACY OF THE PRIMARY ONES
%PROGRAM BY Prof. RAFI MUHANNA AND Prof.  ROBERT MULLEN   only new solution
% for speed
% addpath 'G:\My Drive\Documents\MATLAB\Intlab_V8\Intlab_V8'
%  startintlab();
 intvalinit('DisplayInfsup');
clear
clc
format long
global strainflag;  %this is set to 1 if strain values are needed.
global out;
global nn;
global nnA;
global ne;
global xA;
global yA;
global ndof;
global ndof1;
global ndofA;
global resxA;
global resyA;
global fxA;
global fyA;
global BetaxA;
global BetayA;
global   A;
global E;
global elementA;
global alfaA;
global ng;
global Gamma;
global betaxA;
global betayA;
global Alfa;
global gammao;
global gn;
global gamma1;
global iii;
tic
intvalinit('DisplayInfsup');
name="popova20group";
inp= fopen(name+'.inp','r');
out =fopen(name+'.out','w');
strainflag=0;   % if set to one, strains will also be calculated
fid=out;

[nn,nnA,ne,xA,yA,ndof,ndofA, resxA,resyA,fxA,fyA,BetaxA,BetayA, ...
    A,E,elementA,alfaA,ng,Gamma,betaxA,betayA,Alfa,gammao,gn,gamma1]=readtruss2(name,inp,out);
ivalue=8;
angle=0.;
ndof1=nnA*2;
ndofnew = 2*ne+2*nnA+2*ne+2*ne;
[us]=ifep2(out,ivalue,angle,nn,nnA,ne,xA,yA,ndof,ndof1,ndofA, resxA,resyA,fxA,fyA,BetaxA,BetayA, ...
    A,E,elementA,alfaA,ng,Gamma,betaxA,betayA,Alfa,gammao,gn,gamma1);



nforce=0;  % not used anymore
% [option,uc,Cforce,Cstrain ]=groupcombo (fid,ne,ng,nnA,nforce,E,A,Alfa,Gamma,gamma1,ForceA,xA,yA,elementA,gn,resxA,resyA);
 option= 2;  % no combinations run
%outifep2(fid,strainflag,option,nnA,ndof,ndofA,ndofnew,ng,gamma1,gn,ne,elementA,Alfa,E,A,Gamma,uu,u1,us,ugs,uc,Cforce,Cstrain);
%allocate space for optimizaation program with an added angle value
nvars=1;
lb=zeros(nvars,1);
ub=lb;
x0=lb;

%set angle range
arange=pi/4.; 
% set to +/ 45 degrees
x0(nvars)=0;
lb(nvars)=-arange;
up(nvars)=arange;
fprintf(out,'hybrid1 particle swarm - fixed point one angle with range %f\n',arange);
% start analysis clock

options = optimoptions('particleswarm','FunctionTolerance',1.e-6,'Display','iter');
[x,fval,exitflag,output] =particleswarm(@funx,nvars,lb,ub,options) ;%defines a set of lower and upper bounds on the design variables, x, so that a solution is found in the range lbxub. (Set Aeq=[] and beq=[] if no linear equalities exist.)
fprintf(out,' particle swarm solution min %15.10e max %15.10e  fvalue min %15.10e  time %s function evals %d\n', min(x(1)),max(x(1)),fval/10000,toc,iii);
iii=0;
tic;

[x,fval,exitflag,output] =particleswarm(@funy,nvars,lb,ub,options) ;%defines a set of lower and upper bounds on the design variables, x, so that a solution is found in the range lbxub. (Set Aeq=[] and beq=[] if no linear equalities exist.)
fprintf(out,' particle swarm solution min %15.10e max %15.10e  fvalue min %15.10e  time %s function evals %d\n', min(x(1)),max(x(1)),-fval/10000.,toc,iii);


function disp=funx(optvar)
global iii;
global out;
global nn;
global nnA;
global ne;
global xA;
global yA;
global ndof;
global ndof1;
global ndofA;
global resxA;
global resyA;
global fxA;
global fyA;
global BetaxA;
global BetayA;
global   A;
global E;
global elementA;
global alfaA;
global ng;
global Gamma;
global betaxA;
global betayA;
global Alfa;
global gammao;
global gn;
global gamma1;
iii=iii+1;
ivalue=2*ne+2*nnA-1; % last node x direciton
[us]=ifep2(out,ivalue,optvar,nn,nnA,ne,xA,yA,ndof,ndof1,ndofA, resxA,resyA,fxA,fyA,BetaxA,BetayA, ...
    A,E,elementA,alfaA,ng,Gamma,betaxA,betayA,Alfa,gammao,gn,gamma1);
disp=inf(us(ivalue));
return
end
function disp=funy(optvar)
global iii;
global out;
global nn;
global nnA;
global ne;
global xA;
global yA;
global ndof;
global ndof1;
global ndofA;
global resxA;
global resyA;
global fxA;
global fyA;
global BetaxA;
global BetayA;
global   A;
global E;
global elementA;
global alfaA;
global ng;
global Gamma;
global betaxA;
global betayA;
global Alfa;
global gammao;
global gn;
global gamma1;
iii=iii+1;
ivalue=2*ne+2*nnA-1; % last node x direciton
[us]=ifep2(out,ivalue,optvar,nn,nnA,ne,xA,yA,ndof,ndof1,ndofA, resxA,resyA,fxA,fyA,BetaxA,BetayA, ...
    A,E,elementA,alfaA,ng,Gamma,betaxA,betayA,Alfa,gammao,gn,gamma1);
disp=-sup(us(ivalue));
return
end
function [us]=ifep2(out,ivalue,optvar,nn,nnA,ne,xA,yA,ndof,ndof1,ndofA, resxA,resyA,fxA,fyA,BetaxA,BetayA, ...
    A,E,elementA,alfaA,ng,Gamma,betaxA,betayA,Alfa,gammao,gn,gamma1)
global strainflag;
%--------------------------------------------------------------------------
Del_alpha = infsup(-alfaA,alfaA);

fid=out;
nforce=0;
ForceA=infsup(zeros(ndof1,1),zeros(ndof1,1)); %allocate space for ForceA  RLM 9/8/22
for i=1:nnA
  if(abs(fxA(i))>0.0)
      nforce = nforce+1;
  end
  if(abs(fyA(i))>0.0)
      nforce = nforce+1;
  end
  ForceA(2*i-1) = BetaxA(i)*(cos(optvar)*fxA(i)+sin(optvar)*fyA(i));
  ForceA(2*i)   = BetayA(i)*(fyA(i)*cos(optvar)-sin(optvar)*fxA(i));
end
nE = 0; %number of elements with uncertain E
for i=1:ne
    if(alfaA(i)>0.0)
        nE = nE+1;
    end
end


ndof1=2*nnA;
%nodfnew is the size of the combined stiffness matrix
ndofnew = ndof+ndof1+ndof;
if (strainflag == 1) 
    ndofnew=ndofnew+ndof;
end
%fprintf(1,'ndof = %d ndof1=%d nlamda= %d ndofnew=%d, strainflag = %d\n',ndof,ndof1,nnA,ndofnew,strainflag);
%a = zeros(ne,ndofnew);
a = zeros (ne,ndof);
aA = zeros(ne,ndofA);
Le=zeros(ne,1);
cosine=zeros(ne,1);
sine=zeros(ne,1);
Doo=zeros(ne,1); % changeto vector from diaganol 
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
  
  a(e,dof)=[1 -1];
 
  aA(e,dofA)= [c s -c -s]' ; %R*[1 0 -1  0]'; % again   why not just vector [c s -c -c]'
  %____________diagonal member matrix_______________________

 Doo(e)=A(e)*E(e)/le; 
 Do(e)=A(e)*E(e)/le*infsup(1.0,1.0);
 D(e)=(A(e)*E(e)/le*Alfa(e));
 DoD(e)=Do(e)-D(e) ; 
 Dom(e) = (E(e)*A(e)/le)*(1.+alfaA(e)*gammao(e));
 Dom2(e,e) = (E(e)*A(e)/le)*(1.+alfaA(e)*gammao(e));
end
aa=zeros(ne,4*ne+ndof1+strainflag*ndof);
aa(:,1:ndof)=a;
%----------------formulation of the new approach-------------------------


% one constraint term for each EBE node
i2=1;
cons1=zeros(ndof+ndof1,ndof);
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
    cons2=zeros(ndof,ne);
    for i1=1:ne
    i2=i2+1;
    cons2(i2,i1)=-1/Le(i1);
    i2=i2+1;
    cons2(i2,i1)=1./Le(i1);
    end
istrain=ndof+ndof1+ndof+ndof;
end
inostrain=ndof+ndof1+ndof;
%two centered matriese K no group and Km with 1+gamma*delta term
if (strainflag == 1)
    Km=zeros(istrain,istrain);
    Km(1:ndof,inostrain+1:inostrain+ne)=cons2;
    Km(inostrain+1:inostrain+ne,1:ndof)=cons2';
    Km(inostrain+1:inostrain+ne,inostrain+ne+1:inostrain+ndof)=-eye(ne,ne);
    Km(inostrain+ne+1:inostrain+ndof,inostrain+1:inostrain+ne)=-eye(ne,ne);
    K=Km;
else
    Km=zeros(inostrain,inostrain);
    K=Km;
end 

    Km(1:ndof,1:ndof)=a'*diag(Dom)*a;
       
    Km(1:ndof+ndof1,ndof+ndof1+1:inostrain)=cons1;
    Km(ndof+ndof1+1:inostrain,1:ndof+ndof1)=cons1';
    K(1:ndof,1:ndof)=a'*diag(Doo)*a;
    K(1:ndof+ndof1,ndof+ndof1+1:inostrain)=cons1;
    K(ndof+ndof1+1:inostrain,1:ndof+ndof1)=cons1';
bn2 = 0;
for i=1:nnA
     if(resxA(i)==0)
           bn2 = bn2+1;
           ifix3(bn2) = ndof+ 2*i-1;
     end
     if(resyA(i)==0)
          bn2 = bn2+1;
          ifix3(bn2) = ndof+ 2*i;
    end
end
Km(ifix3,:)=zeros(bn2,4*ne+ndof1+strainflag*ndof);
Km(:,ifix3)=zeros(4*ne+ndof1+strainflag*ndof,bn2);
Km(ifix3,ifix3)=eye(length(ifix3));
K(ifix3,:)=zeros(bn2,4*ne+ndof1+strainflag*ndof);
K(:,ifix3)=zeros(4*ne+ndof1+strainflag*ndof,bn2);
K(ifix3,ifix3)=eye(length(ifix3));
%check for centered solution
Cm=inv(Km);
C=inv(K);
%  create force matrix for assembled nodes
ffm=zeros(4*ne+ndof1+strainflag*ndof,1);
ff=ffm*infsup(0.,0.);
ff(ndof+1:ndof+ndof1,1)=ForceA;
%ffm(ndof+1:ndof+ndof1,1)=mid(ForceA);
% uu is midpoint solution
%uu=C*ffm;
uu=0;

% %*********************   Initial Enclosure  *******************
% w=ones(ne,1);
% w1=w-mag(diag(DoD))*mag(aa*C*aa')*w;
% w2=mag(diag(DoD))*mag(aa*C*ff);
% % allocate alpha for efficency
% alpha=zeros(ne,1);
% for i=1:ne
%     alpha(i)=(w2(i)/w1(i));
% end
% alphamax=max(alpha(:,:));           
% dd=infsup(-alphamax*w,alphamax*w);  
% u1 = C*ff + C*aa'*dd;
%************************************************************
 % for i=1:bn2    
 %   aa(:,ifix2(i))=0; 
 % end 
% v=aa*u1;
% d=diag(DoD)*v;
% v(:,1)=v;
% d(:,1)=d; 
% for i=1:10
%     v(:,i+1)=intersect(((aa*C)*ff+(aa*C*aa')*d(:,i)),(v(:,i)));
%     d(:,i+1)=intersect(diag(DoD)*v(:,i+1),(d(:,i)));
% end
% u1=(C*ff)+(C*aa')*d(:,11); 
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
%vg = (aa*C)*ff;

for i = 1:5
%    vg = hull(((aa*C)*ff - (((aa*C*aa'*diag(Do))*diag(vg))*LambdaE)*Del_alpha),vg);
    vs = hull(((aa*Cm)*ff-((aa*Cm*aa'*diag(Do))*diag(vs)*LambdaE*Del_alpha)+ ...
        ((aa*Cm*aa'*diag(Do))*diag(vs)*LambdaG*Del_gamma)),vs);
end
us = (Cm*ff-((Cm*aa'*diag(Do))*diag(vs)*LambdaE*Del_alpha)+((Cm*aa'*diag(Do))*diag(vs)*LambdaG*Del_gamma));

%ugs = C*ff + C*aa'*diag(Do)*diag(vg)*LambdaE*Del_alpha; 
%IFEAtime = toc





%print results
% uc,Cforce,Cstrain are the end-point combitation results
% uu is midpoint solution
% u1 is no group old intersection resutls
% us is group solution hull
% ugs is no group solution hull
%outifep(fid,strainflag,option,nnA,ndof,ndofA,ndofnew,ng,gamma1,gn,ne,elementA,Alfa,E,A,Gamma,uu,u1,us,ugs,uc,Cforce,Cstrain);
return

end
