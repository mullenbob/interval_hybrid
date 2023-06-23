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
out =fopen(name+'speed.out','w');
strainflag=0;   % if set to one, strains will also be calculated
fid=out;

[nn,nnA,ne,xA,yA,ndof,ndofA, resxA,resyA,fxA,fyA,BetaxA,BetayA, ...
    A,E,elementA,alfaA,ng,Gamma,betaxA,betayA,Alfa,gammao,gn,gamma1]=readtruss2(name,inp,out);
ivalue=8;
angle=0.;
ndof1=nnA*2;
ndofnew = 2*ne+2*nnA+2*ne+2*ne;
%[us]=ifep2(out,ivalue,angle,nn,nnA,ne,xA,yA,ndof,ndof1,ndofA, resxA,resyA,fxA,fyA,BetaxA,BetayA, ...
%    A,E,elementA,alfaA,ng,Gamma,betaxA,betayA,Alfa,gammao,gn,gamma1);



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
tic
iii=0;
options = optimoptions('particleswarm','FunctionTolerance',1.e-6,'Display','iter');
[x,fval] =particleswarm(@funx,nvars,lb,ub,options) ;%defines a set of lower and upper bounds on the design variables, x, so that a solution is found in the range lbxub. (Set Aeq=[] and beq=[] if no linear equalities exist.)
fprintf(out,' particle swarm solution min %15.10e max %15.10e  fvalue min %15.10e  time %s function evals %d\n', min(x(1)),max(x(1)),fval/10000,toc,iii);
% iii=0;   make second value the total
% tic;     make second value the total

[x,fval] =particleswarm(@funy,nvars,lb,ub,options) ;%defines a set of lower and upper bounds on the design variables, x, so that a solution is found in the range lbxub. (Set Aeq=[] and beq=[] if no linear equalities exist.)
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
%add matrices global to save from needing to resolve
global Cm
global aa
global Do
global iii

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
if (iii <2) 
    
    [aa,Do,Cm]=createIM (ne,nnA,ndof,ndof1,strainflag,ndofA,xA,yA,elementA,A,E,resxA,resyA,alfaA,gammao);
end
%  create force matrix for assembled nodes
ffm=zeros(4*ne+ndof1+strainflag*ndof,1);
ff=ffm*infsup(0.,0.);
ff(ndof+1:ndof+ndof1,1)=ForceA;
%ffm(ndof+1:ndof+ndof1,1)=mid(ForceA);
% uu is midpoint solution
%uu=C*ffm;
% uu=0;

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
%LambdaE = eye(ne,ne);  %may not need if E is EBE but keep for reading equations
%--------------------------------------------------------------------------
vs = aa*Cm*ff;


for i = 1:5

    vs = hull(((aa*Cm)*ff-((aa*Cm*aa'*diag(Do.*vs))*Del_alpha)+ ...
        ((aa*Cm*aa'*diag(Do.*vs))*LambdaG*Del_gamma)),vs);
    
end
us = (Cm*ff-((Cm*aa'*diag(Do.*vs))*Del_alpha)+((Cm*aa'*diag(Do.*vs))*LambdaG*Del_gamma));

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
