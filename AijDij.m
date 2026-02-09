clear
clc

%% Aij
SCN1 = xlsread('mic_indiv_scn1.csv');
SCN2 = xlsread('mic_indiv_scn2.csv');
SCN3 = xlsread('mic_indiv_scn3.csv');
SCN4 = xlsread('mic_indiv_scn4.csv');
SCN5 = xlsread('mic_indiv_scn5.csv');
%MIC = [0.949,0.935,0.9895,0.968,0.969],

[a,b]=size(SCN1);
A1=zeros(a,b);
for i=1:a
    for j=1:b
        if SCN1(i,j)>=0.949
            A1(i,j)=1;
        end
    end
end

[a,b]=size(SCN2);
A2=zeros(a,b);
for i=1:a
    for j=1:b
        if SCN2(i,j)>=0.935
            A2(i,j)=1;
        end
    end
end

[a,b]=size(SCN3);
A3=zeros(a,b);
for i=1:a
    for j=1:b
        if SCN3(i,j)>=0.9895
            A3(i,j)=1;
        end
    end
end

[a,b]=size(SCN4);
A4=zeros(a,b);
for i=1:a
    for j=1:b
        if SCN4(i,j)>=0.968
            A4(i,j)=1;
        end
    end
end

[a,b]=size(SCN5);
A5=zeros(a,b);
for i=1:a
    for j=1:b
        if SCN5(i,j)>=0.969
            A5(i,j)=1;
        end
    end
end

%% dij

coordinatesPath = fullfile('..\data',['scn_',ticker, '_coordinates.txt']);
node_coordinates = readmatrix(coordinatesPath, 'Range', 2); 

nodes = node_coordinates(:,1)+1;
theta = node_coordinates(:,3);
r = node_coordinates(:,4);
x = node_coordinates(:,5);
y = node_coordinates(:,6);

d_H_all = zeros(size(nodes));
d_E_all = zeros(size(nodes));

for i = 1:length(nodes)-1
    for j = i+1:length(nodes)
        d_H_all(i,j) = hyp_dist(r(i,1), theta(i,1), r(j,1), theta(j,1));
        d_H_all(j,i) = d_H_all(i,j);
        d_E_all(i,j) = sqrt((x(j,1) - x(i,1))^2 + (y(j,1) - y(i,1))^2);
        d_E_all(j,i) = d_E_all(i,j);
    end
end