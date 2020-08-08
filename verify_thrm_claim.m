% verifies thereom 2.1 (checkunconditionallystabilityfull, where S is S_{ii}) and claim 3.1 (getconsistvars) in the paper 

% (3.6)
load('si2nd.mat')
lam=3/872;
[S,flag]=checkunconditionallystabilityfull(gamma,theta,lam);
[A,B,C,D,E,F,G,H,I]=getconsistvars(gamma,theta,1);

% (3.7)
load('si3rd.mat')
lam=18/28567;
[S,flag]=checkunconditionallystabilityfull(gamma,theta,lam);
[A,B,C,D,E,F,G,H,I]=getconsistvars(gamma,theta,1);


% (4.5)
load('scheme1si.mat')
lam=3/20;
[S,flag]=checkunconditionallystabilityfull(gamma,theta,lam);
[A,B,C,D,E,F,G,H,I]=getconsistvars(gamma,theta,1);

% (4.6)
load('scheme1125si.mat')
lam=1/235;
[S,flag]=checkunconditionallystabilityfull(gamma,theta,lam);
[A,B,C,D,E,F,G,H,I]=getconsistvars(gamma,theta,1);
