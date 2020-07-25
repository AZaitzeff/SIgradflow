%%
%1d Allen-Cahn

order=2;


%error1 is Table 2
onedac


%%
%2d Allen-Cahn (takes a long time to run)
addpath('allencahn/')

order=2;
sieac %Creates 'true' solution
mssiac

%error1 is Table 3
createerrortable


%%
%Cahn-Hillard (takes a really long time to run)
addpath('cahnhilliard/')

order=2;
siech %Creates 'true' solution
mssich

%error1 is Table 4
createerrortablemlch
%%
% Pourous medium equation
addpath('PME/')
trbdf2 % creates 'true solution'


order = 2;
%error1 is Table 5
PMEsigf

%%
addpath('heamobility/')


order = 2;
%error1 is Table 6
mshqwm


%%
% Pourous medium equation
addpath('PME/')
trbdf2 % creates 'true solution'


order = 2;
%error1 is Table 7
PMEwm53

%%
%1d cahn hillard (takes time to run)
addpath('cahnmobility1d/') 
siechwm1d % creates 'true solution'
order = 2;

msch1d
%error1 is Table 8
errortable1d
