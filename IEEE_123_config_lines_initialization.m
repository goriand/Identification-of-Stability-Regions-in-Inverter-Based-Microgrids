clear all

load('linedata.mat')

% % transform
% 135 115
% 149 116
% 151 117
% 152 118
% 160 119
% 197 120
% 250 121
% 300 122
% 450 123



linedata.NodeA(114) = 115;
linedata.NodeA(115) = 116;
linedata.NodeB(51) = 117;
linedata.NodeA(116) = 118;
linedata.NodeA(117) = 119;
linedata.NodeA(118) = 120;
linedata.NodeB(32) = 121;
linedata.NodeB(108) = 122;
linedata.NodeB(99) = 123;

% N_lines = size(linedata(:,1),1);
% for i=1:N_lines
%     nabla(i,graph(i,1)) = 1;
%     nabla(i,graph(i,2)) = -1;
% end

%% Switches
% Adding the switches connected - new 5 lines (exluding node 610)

% load('switchdata.mat')
% switchdata = table2array(switchdata([1:3,5:6],1:2));
% nabla = [nabla; zeros(5,N_nodes)];
% Index=sub2ind(size(nabla),N_lines+1:N_lines+5,switchdata(:,1)');
% nabla(Index) = 1;
% Index=sub2ind(size(nabla),N_lines+1:N_lines+5,switchdata(:,2)');
% nabla(Index) = -1;

% switchdata =
% 
%     13   152
%     18   135
%     60   160
%     97   197

% nabla(:,13) = nabla(:,13) + nabla(:,118);
% nabla(:,118) = [];
% nabla(:,18) = nabla(:,18) + nabla(:,115);
% nabla(:,115) = [];
% nabla(:,60) = nabla(:,60) + nabla(:,119);
% nabla(:,119) = [];
% nabla(:,97) = nabla(:,97) + nabla(:,120);
% nabla(:,120) = [];

linedata.NodeA(116) = 13;
linedata.NodeA(114) = 18;
linedata.NodeA(117) = 60;
linedata.NodeA(118) = 97;

graph = table2array(linedata(:,1:2));
Gt = graph(:,1:2);
G = digraph(Gt(:,1)', Gt(:,2)');
nabla = incidence(G)';




% % Making the incidence matrix nabla
N_nodes = size(nabla,2);
N_lines = size(nabla,1);

% nabla = zeros(N_lines,N_nodes);
% % Index=sub2ind(size(nabla),1:N_lines,graph(:,1)');
% % nabla(Index) = 1;
% % Index=sub2ind(size(nabla),1:N_lines,graph(:,2)');
% % nabla(Index) = -1;

% graph(graph==135) = 115;
% graph(graph==149) = 116;
% graph(graph==151) = 117;
% graph(graph==152) = 118;
% graph(graph==160) = 119;
% graph(graph==197) = 120;
% graph(graph==250) = 121;
% graph(graph==300) = 122;
% graph(graph==450) = 123;

%% Adding the loads
load('spotloadsdata')
% take spot loads locations
ind_loads = table2array(spotloadsdata(1:end-1,1));
% take real power load for each phase abc and approximate as balanced load
P1_load = table2array(spotloadsdata(1:end-1,3));
P2_load = table2array(spotloadsdata(1:end-1,5));
P3_load = table2array(spotloadsdata(1:end-1,7));
Pload_balanced = 1/3*(P1_load + P2_load + P3_load);
% take reactive power load for each phase abc and approximate as balanced load
Q1_load = table2array(spotloadsdata(1:end-1,4));
Q2_load = table2array(spotloadsdata(1:end-1,6));
Q3_load = table2array(spotloadsdata(1:end-1,8));
Qload_balanced = 1/3*(Q1_load + Q2_load + Q3_load);
% add the self addmitances to the nabla
M_loads = length(ind_loads);
nabla_loads = zeros(M_loads,N_nodes);
for i=1:M_loads
    nabla_loads(i,ind_loads(i)) = 1;
end
nabla = [nabla;nabla_loads];



%% Line Impedances
% abc to sequences transformation
a= exp(1i*2*pi/3);
A = [1 1 1; 1 a a^2; 1 a^2 a];
% kms per mile
c = 1.60934;

%config 1
Zp_1 = 0.7680 + 1.6908i;
Zp_1 = c*Zp_1;

%config 2
Zp_2 = 0.4725 + 0.8005i;
Zp_2 = c*Zp_2;
%config 3
Z=[0.4615+1j*1.0651,0.1535+1j*0.3849,0.1580+1j*0.4236
0,0.4576+1j*1.0780,0.1560+1j*0.5017
0,0,0.4666+1j*1.0482];
Z = Z + transpose(Z);
Z_tr = A*Z*inv(A);
Zp_3 = Z_tr(2,2);
Zp_3 = c*Zp_3;

%config 4
 Z=[0.4615+1j*1.0651,0.1580+1j*0.4236,0.1535+1j*0.3849
 0,0.4666+1j*1.0482,0.1560+1j*0.5017
 0,0,0.4576+1j*1.0780];
 Z = Z + transpose(Z);
 Z_tr = A*Z*inv(A);
 Zp_4 = Z_tr(2,2);
Zp_4 = c*Zp_4;
 %config 5
 Z=[0.4666+1j*1.0482,0.1560+1j*0.5017,0.1580+1j*0.4236
 0,0.4576+1j*1.0780,0.1535+1j*0.3849
 0,0,0.4615+1j*1.0651];
 Z = Z + transpose(Z);
 Z_tr = A*Z*inv(A);
 Zp_5 = Z_tr(2,2);
Zp_5 = c*Zp_5;
%config 6
 Z=[0.4576+1j*1.0780,0.1535+1j*0.3849,0.1560+1j*0.5017
 0,0.4615+1j*1.0651,0.1580+1j*0.4236
 0,0,0.4666+1j*1.0482];
Z = Z + transpose(Z);
 Z_tr = A*Z*inv(A);
 Zp_6 = Z_tr(2,2);
Zp_6 = c*Zp_6;

%config 7
 Z=[0.4576+1j*1.0780,0.0000+1j*0.0000,0.1535+1j*0.3849
 0,0.0000+1j*0.0000,0.0000+1j*0.0000
 0,0,0.4615+1j*1.0651];
Z = Z + transpose(Z);
 Z_tr = A*Z*inv(A);
 Zp_7 = Z_tr(2,2);
Zp_7 = c*Zp_7; 

%config 8
 Z=[0.4576+1j*1.0780,0.1535+1j*0.3849,0.0000+1j*0.0000
 0,0.4615+1j*1.0651,0.0000+1j*0.0000
 0,0,0.0000+1j*0.0000];
Z = Z + transpose(Z);
 Z_tr = A*Z*inv(A);
 Zp_8 = Z_tr(2,2);
Zp_8 = c*Zp_8;
%config 9
 Z=[1.3292+1j*1.3475,0.0000+1j*0.0000,0.0000+1j*0.0000
 0,0.0000+1j*0.0000,0.0000+1j*0.0000
 0,0,0.0000+1j*0.0000];
Z = Z + transpose(Z);
 Z_tr = A*Z*inv(A);
 Zp_9 = Z_tr(2,2);
Zp_9 = c*Zp_9;
%config 10
 Z=[0.0000+1j*0.0000,0.0000+1j*0.0000,0.0000+1j*0.0000
 0,1.3292+1j*1.3475,0.0000+1j*0.0000
 0,0,0.0000+1j*0.0000];
Z = Z + transpose(Z);
 Z_tr = A*Z*inv(A);
 Zp_10 = Z_tr(2,2);
Zp_10 = c*Zp_10;
%config 11
 Z=[0.0000+1j*0.0000,0.0000+1j*0.0000,0.0000+1j*0.0000
 0,0.0000+1j*0.0000,0.0000+1j*0.0000
 0,0,1.3292+1j*1.3475];
Z = Z + transpose(Z);
 Z_tr = A*Z*inv(A);
 Zp_11 = Z_tr(2,2);
Zp_11 = c*Zp_11;
 %config 12
 Z=[1.5209+1j*0.7521,0.5198+1j*0.2775,0.4924+1j*0.2157
 0,1.5329+1j*0.7162,0.5198+1j*0.2775
 0,0,1.5209+1j*0.7521];
Z = Z + transpose(Z);
 Z_tr = A*Z*inv(A);
 Zp_12 = Z_tr(2,2);
Zp_12 = c*Zp_12;

Z_line_per_km = [Zp_1; Zp_2; Zp_3; Zp_4; Zp_5; Zp_6; Zp_7; Zp_8; Zp_9; Zp_10; Zp_11; Zp_12];
Z_line = 0.0003048*table2array(linedata(:,3)).*Z_line_per_km(table2array(linedata(:,4)));
clear Zp_1 Zp_2 Zp_3 Zp_4 Zp_5 Zp_6 Zp_7 Zp_8 Zp_9 Zp_10 Zp_11 Zp_12



% for i=1:N_lines
%     nabla(i,graph(i,1)) = 1;
%     nabla(i,graph(i,2)) = -1;
% end



