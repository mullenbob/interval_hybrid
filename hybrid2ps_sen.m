% hybrid 2   particle swarm - sensitivity method
% written 2023-06-08
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
global inp;
global out;
global iii;
global nnA;
global ne;
global ndof1;
tic
starttime=cputime;
intvalinit('DisplayInfsup');
name="popova100group";
inp= fopen(name+'.inp','r');
out =fopen(name+'_ps_sen.out','w');
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
nvars=1;
lb=zeros(nvars,1);
ub=lb;
%set angle range
arange=pi/4.; 
% set to +/ 45 degrees
lb(nvars)=-arange;
ub(nvars)=arange;

fprintf(out,'hybrid2 particle swarm - sensitivity one angle with range %f\n',arange);
% start analysis clock
tic
starttime=cputime;
iii=0;
options = optimoptions('particleswarm','FunctionTolerance',1.e-6);
[x,fval] =particleswarm(@funx,nvars,lb,ub,options) ;%defines a set of lower and upper bounds on the design variables, x, so that a solution is found in the range lbxub. (Set Aeq=[] and beq=[] if no linear equalities exist.)
lb(nvars)=-arange;
up(nvars)=arange;
fprintf(out,' particle swarm solution angle  %15.10e   fvalue min %15.10e  time %s cputime %f function evals %d\n', x(1),fval,toc,cputime-starttime,iii);
% supress the function evals split   iii=0;
% supress time split so second time is the total time   tic

[x,fval] =particleswarm(@funy,nvars,lb,ub,options) ;%defines a set of lower and upper bounds on the design variables, x, so that a solution is found in the range lbxub. (Set Aeq=[] and beq=[] if no linear equalities exist.)
fprintf(out,' particle swarm solution angle  %15.10e   fvalue min %15.10e  time %s cputime %f  function evals %d\n', x(1),-fval,toc,cputime-starttime,iii);


function disp=funx(optvar)
global iii;

iii=iii+1;

[temp]=trusssolve(optvar);
disp=temp(2);
return
end
function disp=funy(optvar)
global iii;

iii=iii+1;

temp=trusssolve(optvar);
disp=-temp(1);
return
end
function disp=trusssolve(angle)
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
% build force vector
ii=0;
for i=1:nnd
ii=ii+1;
force(ii)=fxA(i)*cos(angle)+fyA(i)*sin(angle);
ii=ii+1;
force(ii)=fyA(i)*cos(angle)-fxA(i)*sin(angle);
end

CC=inv(K);
R0=CC*(force');
temp=R0(2*nnd-1);
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

for npoint=1:nel
 
v1=A1*R0;
%change to increating only alpha
Dint=diag(D1);
Dint(npoint,npoint)=Dint(npoint,npoint)*(1+alphaA(npoint));
D0=diag(D1);
DD=(D0-Dint)*v1;

%pull out force term
v0=A1*CC*transpose(force);
for i1=1:1
    
v1=(v0+A1*CC*transpose(A1)*DD);
DD=(D0-Dint)*v1;
end


% split multipication into two parts for tighter bounds
R1=CC*transpose(force)+CC*transpose(A1)*DD;

L(:,npoint)=sign(R1-R0);


end
%add in group sensitivity
if (nintvgroup > 0)
    for ng=1:nintvgroup
 
v1=A1*R0;
%change to increating only alpha
Dint=diag(D1);
for i1=1:nel
    if groupid(i1) == ng
Dint(i1,i1)=Dint(i1,i1)*(1+gamma(ng));
    end
end
D0=diag(D1);
DD=(D0-Dint)*v1;

%pull out force term
v0=A1*CC*transpose(force);
for i1=1:1
    
v1=(v0+A1*CC*transpose(A1)*DD);
DD=(D0-Dint)*v1;
end


% split multipication into two parts for tighter bounds
R1=CC*transpose(force)+CC*transpose(A1)*DD;



L(:,nel+ng)=sign(R1-R0);
    end
end
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
%% find sensitivity solution for top right node R0(2*nnd-1);
Dint=diag(D1);
if nintvgroup == 0
for npoint=1:nel
Dint(npoint,npoint)=Dint(npoint,npoint)*(1+alphaA(npoint)*L(2*nnd-1,npoint));
end
else
    for npoint=1:nel
        ig=groupid(npoint);
        if (ig == 0 || gamma(ig) == 0)
            Dint(npoint,npoint)=Dint(npoint,npoint)*(1+alphaA(npoint)*L(2*nnd-1,npoint));
        else
Dint(npoint,npoint)=Dint(npoint,npoint)*(1+alphaA(npoint)*L(2*nnd-1,npoint))* ...
    (1+gamma(ig)*L(2*nnd-1,nel+ig));
 
    end
    end
end
KK2=transpose(A1)*Dint*A1;

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
%lower bound
Dint=diag(D1);
if nintvgroup == 0
for npoint=1:nel
Dint(npoint,npoint)=Dint(npoint,npoint)*(1-alphaA(npoint)*L(2*nnd-1,npoint));
end
else
    for npoint=1:nel
        ig=groupid(npoint);
        if (ig == 0 || gamma(ig) == 0)
            Dint(npoint,npoint)=Dint(npoint,npoint)*(1-alphaA(npoint)*L(2*nnd-1,npoint));
        else
Dint(npoint,npoint)=Dint(npoint,npoint)*(1-alphaA(npoint)*L(2*nnd-1,npoint))* ...
    (1-gamma(ig)*L(2*nnd-1,nel+ig));
 
    end
    end
end

KK2=transpose(A1)*Dint*A1;
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

disp=[R2(2*nnd-1),R3(2*nnd-1)];
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
fprintf(out,'Program hybrid interval  particle swarm - sensitivity \n input file name %s \n',name);
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
