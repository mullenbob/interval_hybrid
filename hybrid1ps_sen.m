%function [ ul,uu,timex,feval]=hybrid1ps_sen(name)
% hybrid 2   particle swarm - sensitivity method
% written 2023-06-08
%LINEAR INTERVAL FINTE ELEMENT ANALYSIS OF PLANE TRUSS
%SECONDARY VARIABLES WITH THE SAME ACCURACY OF THE PRIMARY ONES
%PROGRAM BY Prof. RAFI MUHANNA AND Prof.  ROBERT MULLEN   only new solution
% for speed
% addpath 'G:\My Drive\Documents\MATLAB\Intlab_V8\Intlab_V8'
%  startintlab();
clear
clearvars -global
clc
format long
global strainflag;  %this is set to 1 if strain values are needed.
global inp;
global out;
global iii;
global nnA;
global ndof1;
tic
starttime=cputime;
%intvalinit('DisplayInfsup');
name="popova20group";
inp= fopen(name+'.inp','r');
out =fopen(name+'_ps_sen.out','w');
strainflag=1;   % if set to one, strains will also be calculated
fid=out;
%set angle range
arange=pi/4.; 
fprintf(out,'hybrid1 particle swarm - sensitivity one angle with range %f\n   Run on : %s\n',arange,datetime('now'));
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

% set to +/ 45 degrees
lb(nvars)=-arange;
ub(nvars)=arange;

fprintf(out,'hybrid1 particle swarm - sensitivity one angle with range %f\n   Run on : %s\n',arange,datetime('now'));
% start analysis clock
tic
starttime=cputime;
iii=0;
options = optimoptions('particleswarm','FunctionTolerance',1.e-6);
[x,fval] =particleswarm(@funx,nvars,lb,ub,options) ;%defines a set of lower and upper bounds on the design variables, x, so that a solution is found in the range lbxub. (Set Aeq=[] and beq=[] if no linear equalities exist.)
lb(nvars)=-arange;
ub(nvars)=arange;
fprintf(out,' particle swarm solution angle  %15.10e   fvalue min %15.10e  time %s cputime %f function evals %d\n', x(1),fval,toc,cputime-starttime,iii);
% supress the function evals split   iii=0;
% supress time split so second time is the total time   tic

[x,fval2] =particleswarm(@funy,nvars,lb,ub,options) ;%defines a set of lower and upper bounds on the design variables, x, so that a solution is found in the range lbxub. (Set Aeq=[] and beq=[] if no linear equalities exist.)
fprintf(out,' particle swarm solution angle  %15.10e   fvalue min %15.10e  time %s cputime %f  function evals %d\n', x(1),-fval2,toc,cputime-starttime,iii);
% ul=fval;
% uu=-fval2;
% timex=cputime-starttime;
% feval=iii;
% return
% end
function disp=funx(optvar)
global iii;

iii=iii+1;

[temp]=trusssolve(optvar,iii,2);
disp=temp(2);
return
end
function disp=funy(optvar)
global iii;

iii=iii+1;

temp=trusssolve(optvar,iii,1);
disp=-temp(1);
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
