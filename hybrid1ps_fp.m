% hybrid 1   particle swarm -fix-point methods using
%LINEAR INTERVAL FINTE ELEMENT ANALYSIS OF PLANE TRUSS
%SECONDARY VARIABLES WITH THE SAME ACCURACY OF THE PRIMARY ONES
%PROGRAM BY Prof. RAFI MUHANNA AND Prof.  ROBERT MULLEN   only new solution
% for speed
% addpath 'G:\My Drive\Documents\MATLAB\Intlab_V8\Intlab_V8'
%  startintlab();
 intvalinit('DisplayInfsup');
clear
clearvars -global
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
starttime=cputime;
intvalinit('DisplayInfsup');
name="popova100group";
inp= fopen(name+'.inp','r');
out =fopen(name+'ps_fp.out','w');
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
ub(nvars)=arange;
fprintf(out,'hybrid1 particle swarm - fixed point one angle with range %f\n',arange);
% start analysis clock
tic
starttime=cputime;
iii=0;
options = optimoptions('particleswarm','FunctionTolerance',1.e-6);
[x,fval] =particleswarm(@funx,nvars,lb,ub,options) ;%defines a set of lower and upper bounds on the design variables, x, so that a solution is found in the range lbxub. (Set Aeq=[] and beq=[] if no linear equalities exist.)
fprintf(out,' particle swarm solution min %15.10e max %15.10e  fvalue min %15.10e  time %s cpu time %f function evals %d\n', min(x(1)),max(x(1)),fval,toc,cputime-starttime,iii);
% iii=0;   make second value the total
% tic;     make second value the total

[x,fval] =particleswarm(@funy,nvars,lb,ub,options) ;%defines a set of lower and upper bounds on the design variables, x, so that a solution is found in the range lbxub. (Set Aeq=[] and beq=[] if no linear equalities exist.)
fprintf(out,' particle swarm solution min %15.10e max %15.10e  fvalue max %15.10e  time %s cputime %f function evals %d\n', min(x(1)),max(x(1)),-fval,toc,cputime-starttime,iii);


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
  % note this is beta x for angle problem
   ForceA(2*i)   = BetaxA(i)*(fyA(i)*cos(optvar)-sin(optvar)*fxA(i));
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
ff=zeros(4*ne+ndof1+strainflag*ndof,1);
ff=ff*infsup(0.,0.);
ff(ndof+1:ndof+ndof1,1)=ForceA;
 
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


for i = 1:4
    vs = hull(((aa*Cm)*ff-((aa*Cm*aa'*diag(Do.*vs))*Del_alpha)+ ...
        ((aa*Cm*aa'*diag(Do.*vs))*LambdaG*Del_gamma)),vs);
    
end
us = (Cm*ff-((Cm*aa'*diag(Do.*vs))*Del_alpha)+((Cm*aa'*diag(Do.*vs))*LambdaG*Del_gamma));







%print results
% uc,Cforce,Cstrain are the end-point combitation results
% uu is midpoint solution
% u1 is no group old intersection resutls
% us is group solution hull
% ugs is no group solution hull
%outifep(fid,strainflag,option,nnA,ndof,ndofA,ndofnew,ng,gamma1,gn,ne,elementA,Alfa,E,A,Gamma,uu,u1,us,ugs,uc,Cforce,Cstrain);
return

end
