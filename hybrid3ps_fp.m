% hybrid 3  particle swarm - sensitivity method
% written 2023-06-10 specificaly for braced truss  allows each load angle
% to be independent
%  RLM  2023-06-10
% addpath 'G:\My Drive\Documents\MATLAB\Intlab_V8\Intlab_V8'
%  startintlab();
 intvalinit('DisplayInfsup');
clear
clearvars -global

format long
clc
global strainflag;  %this is set to 1 if strain values are needed.
global inp;
global out;
global iii;
global nnA;
global ne;
global ndof1;
tic
starttime=cputime;
intvalinit('DisplayInfsup');
name="popova5group";
inp= fopen(name+'.inp','r');
out =fopen(name+'_ps_fp_ind_angles.out','w');
if name == 'popova5group'
    nangle=5;
else
    if name == 'popova10group'
        nangle=10;
    else
    if name == 'popova20group'
        nangle=20;
    else
        if name == 'popova100group'
            nangle=10;
        else
            'not the right program for this problem'
        end
    end
    end
end

strainflag=1;   % if set to one, strains will also be calculated
fid=out;
readtrussx(name,inp,out);
%[nn,nnA,ne,xA,yA,ndof,ndofA, resxA,resyA,fxA,fyA,BetaxA,BetayA, ...
%    A,E,elementA,alfaA,ng,Gamma,betaxA,betayA,Alfa,gammao,gn,gamma1]=readtruss2(name,inp,out);
angle=0.;
ndof1=nnA*2;
option= 2;  % no combinations run
%outifep2(fid,strainflag,option,nnA,ndof,ndofA,ndofnew,ng,gamma1,gn,ne,elementA,Alfa,E,A,Gamma,uu,u1,us,ugs,uc,Cforce,Cstrain);
%allocate space for optimizaation program with an added angle value
nvars=nangle
lb=zeros(nvars,1);
ub=lb;
%set angle range
arange=pi/4.; 
% set to +/ 45 degrees
for i1=1:nvars
lb(i1)=-arange;
ub(i1)=arange;
end

fprintf(out,'hybrid3 particle swarm - sensitivity with angles independent angle with range %f\n',arange);
% start analysis clock
tic
starttime=cputime;
iii=0;
options = optimoptions('particleswarm','FunctionTolerance',1.e-6);
[x,fval] =particleswarm(@funx,nvars,lb,ub,options) ;%defines a set of lower and upper bounds on the design variables, x, so that a solution is found in the range lbxub. (Set Aeq=[] and beq=[] if no linear equalities exist.)
lb(nvars)=-arange;
ub(nvars)=arange;
fprintf(out,' funx particle swarm solution angle min %15.10e max % 15.10e  fvalue min %15.10e  time %s  cputime %f function evals %d\n',min(x),max(x),fval,toc,cputime-starttime,iii);
% supress the function evals split   iii=0;
% supress time split so second time is the total time   tic

[x,fval] =particleswarm(@funy,nvars,lb,ub,options) ;%defines a set of lower and upper bounds on the design variables, x, so that a solution is found in the range lbxub. (Set Aeq=[] and beq=[] if no linear equalities exist.)
fprintf(out,' funy particle swarm solution angle min %15.10e max %15.10e  fvalue min %15.10e  time %s CPU time %f function evals %d\n', min(x),max(x),-fval,toc,cputime-starttime,iii);


function disp=funx(optvar)
global iii;
global nnA;
global ne;
iii=iii+1;
ivalue=2*ne+2*nnA-1; % last node x direciton
[us]=ifep2(out,ivalue,optvar,nn,nnA,ne,xA,yA,ndof,ndof1,ndofA, resxA,resyA,fxA,fyA,BetaxA,BetayA, ...
    A,E,elementA,alfaA,ng,Gamma,betaxA,betayA,Alfa,gammao,gn,gamma1);
disp=inf(us(ivalue));
%ivalue=2*ne+2*nnA-1; % last node x direciton
%[temp]=trusssolve(optvar);
%disp=temp(2);
return
end
function disp=funy(optvar)
global iii;
global out;
global nn;
global nnA;
global ne;

iii=iii+1;
%ivalue=2*ne+2*nnA-1; % last node x direciton
%[us]=ifep2(out,ivalue,optvar,nn,nnA,ne,xA,yA,ndof,ndof1,ndofA, resxA,resyA,fxA,fyA,BetaxA,BetayA, ...
 %
 % A,E,elementA,alfaA,ng,Gamma,betaxA,betayA,Alfa,gammao,gn,gamma1);
%temp=trusssolve(optvar);
%disp=-temp(1);
ivalue=2*ne+2*nnA-1; % last node x direciton
[us]=ifep2(out,ivalue,optvar,nn,nnA,ne,xA,yA,ndof,ndof1,ndofA, resxA,resyA,fxA,fyA,BetaxA,BetayA, ...
    A,E,elementA,alfaA,ng,Gamma,betaxA,betayA,Alfa,gammao,gn,gamma1);
disp=sup(us(ivalue));
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
nangle=0;
for i=1:nnA
  if(abs(fxA(i))>0.0)
      nforce = nforce+1;
  end
  if(abs(fyA(i))>0.0)
      nforce = nforce+1;
  end
  if (fxA(i) ~= 0.)
  anglex=angle(nangle);
ii=ii+1;
ForceA(ii)=BetaA(i)*(fxA(i)*cos(anglex)+fyA(i)*sin(anglex));
ii=ii+1;
ForceA(ii)=BetaA(i)*(fyA(i)*cos(anglex)-fxA(i)*sin(anglex));
nangle=nangle+1;
    else
        ii=ii+1;
        ForceA(ii)=0.;
        ii=ii+1;
        ForceA(ii)=0.;
    end
end
 % ForceA(2*i-1) = BetaxA(i)*(cos(optvar)*fxA(i)+sin(optvar)*fyA(i));
  %% note this is beta x for angle problem
   %ForceA(2*i)   = BetaxA(i)*(fyA(i)*cos(optvar)-sin(optvar)*fxA(i));
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

function readtrussx(name,inp,out)
   
% set up global variables to exchange inifor between functions

global tx;
global ty;

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
global ngroups;
global gamma;
global groupid;
global alphaA;
global betax;
global betay;
global nintvload;
global nintvgroup;
global nodeid;
global beta;



% read in data in same forame as 
                                                                            
[mat1] = fscanf(inp,'%d',2);
nnd= mat1(1); %number of nodes - assembled model
nel  = mat1(2); %number of elements
ndof=2*nnd;

%x and y coodinates of nodes for assembled model
% tx = zeros(nnd,1);
% ty = zeros(nnd,1);
% resxA=zeros(nnd,1);
% resyA=zeros(nnd,1);
% fxA=zeros(nnd,1);
% fyA=zeros(nnd,1);

%reading nodal information for the assembled model
[mat1] = fscanf(inp,'%d %f %f %d %d %e %e  %f %f',[9,nnd]);
mat1 = mat1';
tx = mat1(:,2);
ty = mat1(:,3);
%nodal information for assembled model
resxA = mat1(:,4); %x restraint for assembled model
resyA = mat1(:,5); %x restraint for assembled model
%x and y components of forces for assembled model
fxA = mat1(:,6);
fyA = mat1(:,7);
betax = mat1(:,8);
betay = mat1(:,9);
% count number of interval loads
nintvload=0;
for i1=1:nnd
    if (betax(i1) > 0)
        nintvload=nintvload+1;
        nodeid(nintvload)=i1*2-1;
        beta(nintvload)=betax(i1);
    end
    if (betay(i1) > 0) 
        nintvload=nintvload+1;
        nintvload=nintvload+1;
        nodeid(nintvload)=i1*2;
        beta(nintvload)=betay(i1);
    end
end
%reading element information for the assembled model
%tA       = zeros(nel,1);
%tE       = zeros(nel,1);
%telement = zeros(nel,2);
[mat2] = fscanf(inp,'%d %d %d %f %e %f %d',[7,nel]);
mat2 = mat2';
telement = mat2(:,2:3); %nodal connectivity matrix for assembled model
tA       = mat2(:,4);
tE       = mat2(:,5);
alphaA    = mat2(:,6); %uncertainty of E in assembled mode
groupid = mat2(:,7); 
% read in group information
ngroups=fscanf(inp,'%d',1);

nintvgroup=0;
if (ngroups > 0)
  [mat4] = fscanf(inp,'%d %f',[2,ngroups]);

mat4=mat4';
gamma= mat4(:,2);

    for i1=1:ngroups
    if (gamma(i1) > 0)
        nintvgroup=nintvgroup+1;
    end  
    end
end
% Print out input data now 
fprintf(out,'Program hybrid3 interval  particle swarm - sensitivity with independent angles \n input file name %s \n',name);
fprintf(out,'Nodal information - Assembled  model\n');
fprintf(out,'Node   X      Y     Restraints       Fx          Fy           beta-x    beta-y\n');
for i=1:nnd
     fprintf(out,'%2d   %4.1f   %4.1f    %d    %d     %9.1f    %9.1f           %5.3f      %5.3f\n',i,tx(i),ty(i),resxA(i),resyA(i),fxA(i),fyA(i),betax(i),betay(i)); 
end
fprintf(out,'Element information- Assembled model\n');
fprintf(out,'Element    Nodes    Area            E        alfa      group no.    group gamma \n');
for i=1:nel
    temp = 0;
    if groupid(i) > 0 
        temp=gamma(groupid(i));
    end
    fprintf(out,'%2d        %2d  %2d    %f   %.3e   %7.5f         %2d           %7.5f\n',i,telement(i,1),telement(i,2),tA(i),tE(i),alphaA(i),groupid(i), temp);
end
fprintf(out,'Groups Interval Solution\n');
fprintf(out,'Number of Groups\n');
fprintf(out,'%2d\n',ngroups);
fprintf(out,'Group Number        Group Uncertainty   \n');
for i=1:ngroups
    fprintf(out,'%2d                       %7.5f    \n',i,gamma(i)); 
end
return
end
