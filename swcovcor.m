%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%     SPATIALLY-WEIGHTED      %%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%%%   COVARIANCE & CORRELATION    %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%           MATRIX            %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Author: Saeid ESMAEILOGHLI
% February 25, 2019
% ------------------------
% Institut de Recherche en Astrophysique et Planétologie (IRAP) 
% Centre National de la Recherche Scientifique (CNRS),
% Observatoire Midi-Pyrénées (OMP),
% Université Paul Sabatier (UPS),
% 14 Avenue Edouard Belin,
% 31400 Toulouse,
% France.
% ------------------------
% Department of Mining Engineering,
% Isfahan University of Technology (IUT),
% University Boulevard, Esteghlal Square,
% 84156-83111 Isfahan,
% Iran.
% ------------------------
% Phone: (+33) 6 25 47 67 03
% E-mail: saeid@mail.fr
% ORCID: 0000-0002-7786-657X
% Website: http://esmaeiloghli.mining.iut.ac.ir
% GoogleScholar: http://scholar.google.com/citations?user=EZKZcwQAAAAJ&hl=en
% ResearchGate: https://www.researchgate.net/profile/Saeid_Esmaeiloghli
% ACADEMIA: https://cnrs.academia.edu/SaeidEsmaeiloghli


% -------------------------------------------------------------------------

% Author's statements:
%
% This program was designed to calculate the spatially-weighted covariance
% (correlation) matrix according to a distance-based spatio-geological
% weighting function. A work similar in application but different in
% mathematics can be found at:
% Cheng, Qiuming, Greame Bonham-Carter, Wenlei Wang, Shengyuan Zhang,
% Wenchang Li, and Xia Qinglin. "A spatially weighted principal component
% analysis for multi-element geochemical data for mapping locations of
% felsic intrusions in the Gejiu mineral district of Yunnan, China."
% Computers & Geosciences 37, no. 5 (2011): 662-669.
%
% The input parameters include:
%
% data: an Excel file comprising the m * n elemental data matrix, where
% each row (i = 1,...,m) represents a geochemical sample while each column
% (j = 1,...,n) expresses a centered (zero-mean) log- or ilr-transformed
% geochemical element;
%
% locs: an Excel file comprising the sample locations; it is a three-column
% matrix, XYZ, denoting the longitude, latitude, and elevation of the
% geochemical samples, respectively;
%
% refpoints: an Excel file including the 3D spatial locations of potential
% sources; it is a three-column matrix, XYZ, denoting the longitude,
% latitude, and elevation of geochemical reference points, respectively;
%
% wl: an Excel file containing a 1 * m vector to represent the
% geologist-defined ore productivity scores assigned to samples to account
% for the spatial complexities of various lithological units in the study
% area;
%
% wh: an Excel file containing a 1 * m vector to represent the
% geologist-defined ore productivity scores assigned to samples to account
% for the spatial complexities of various hydrothermal alterations zones in
% the study area;
%
% a: an user-defined specific value to tune the relative importance of
% lithology scores;
%
% b: an user-defined specific value to tune the relative importance of
% hydrothermal alteration scores;
%
% eps: an user-specified constant subject to eps > 1.
%
% The program provides the following outputs:
%
% figure 1: a plot indicating the templete of shortest distance factor;
%
% figure 2: a plot indicating the lithological weighting templete;
%
% figure 3: a plot indicating the alteration weighting templete;
%
% figure 4: a plot indicating the templete of overall spatial weighting
% function;
%
% W: an Excel file including a vector which represents the spatial weights;
%
% SWS: an Excel file containing spatially-weighted covariance matrix;
%
% SWR: an Excel file containing spatially-weighted correlation matrix.
%
% To cite this file, this would be an appropriate format:
% Esmaeiloghli, S., Tabatabaei, S.H., Carranza, E.J.M., Hosseini, S.,
% Deville, Y. "Spatially-weighted factor analysis for extraction of
% source-oriented mineralization feature in 3D coordinates of surface
% geochemical signal." Journal of Geochemical Exploration (2020).


%% Load Inputs

clear all
close all
clc

data = xlsread('data');

locs = xlsread('coordinates');

refpoints = xlsread('refpoints');

wl = xlsread('wl');

wh = xlsread('wh');

a = 1;

b = 1;

eps = 2;

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%        NOW        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   YOU CAN RUN THE   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%      PROGRAM      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% A. Calculating the Spatially Weighting Factor (W)

X = locs(:,1);
Y = locs(:,2);
Z = locs(:,3);

rp = refpoints;
n_dim = length(X);
n_r = length(rp);

for i = 1:n_dim
    
    for j = 1:n_r
        
        alldist(i,j) = sqrt(((rp(j,1)-X(i))^2)+((rp(j,2)-Y(i))^2)+...
            ((rp(j,3)-Z(i))^2)); 
        
    end
    
end

for i = 1:n_dim
    
    mindist(i) = min(alldist(i,:));
    
end

for i = 1:n_dim
    
    ws(i) = 1./(log(mindist(i))+eps);% Non-scaled spatial weighting factor.
    
end

ws = ws';
ns_w = ws.*(wl.^a).*(wh.^b);

for i = 1:n_dim
        
    s_w(i) = (ns_w(i)-min(ns_w))./(max(ns_w)-min(ns_w));
        
    w = s_w';                       % Scaled spatial weighting factor.
    
end

xlswrite('W.xls',w);

%% B. Plot Weighting Factors

figure(1)                           % Distance weighting factor plot.

scatter(X,Y,50,ws,'filled');
colorbar
title('Shortest Distance','fontname','Cambria','FontSize', 16)
xlabel('Easting','fontname','Cambria','FontSize', 14)
ylabel('Northing','fontname','Cambria','FontSize', 14)
xlim([408500 411150])
ylim([3303850 3306480])
get(gca,'fontname')
set(gca,'fontname','Cambria')
set(gcf,'color','w')
box on

figure(2)                           % Lithological weighting factor plot.

scatter(X,Y,50,wl,'filled');
colorbar
title('Lithology Weight','fontname','Cambria','FontSize', 16)
xlabel('Easting','fontname','Cambria','FontSize', 14)
ylabel('Northing','fontname','Cambria','FontSize', 14)
xlim([408500 411150])
ylim([3303850 3306480])
get(gca,'fontname')
set(gca,'fontname','Cambria')
set(gcf,'color','w')
box on

figure(3)                           % Alteration weighting factor plot.

scatter(X,Y,50,wh,'filled');
colorbar
title('Hydrothermal Alteration Weight','fontname','Cambria','FontSize', 16)
xlabel('Easting','fontname','Cambria','FontSize', 14)
ylabel('Northing','fontname','Cambria','FontSize', 14)
xlim([408500 411150])
ylim([3303850 3306480])
get(gca,'fontname')
set(gca,'fontname','Cambria')
set(gcf,'color','w')
box on

figure(4)                           % Spatial weighting factor (W) plot.

scatter(X,Y,50,w,'filled');
colorbar
title('Spatial Weighting Function','fontname','Cambria','FontSize', 16)
xlabel('Easting','fontname','Cambria','FontSize', 14)
ylabel('Northing','fontname','Cambria','FontSize', 14)
xlim([408500 411150])
ylim([3303850 3306480])
get(gca,'fontname')
set(gca,'fontname','Cambria')
set(gcf,'color','w')
box on

%% C. Calculating the Spatially Weighted Covariance Matrix

% C.1. Calculating the weighted mean value for each variable

j_dim = length(data(1,:));
i_dim = length(data(:,1));

for j = 1:j_dim

    for i = 1:i_dim
       
        a(i,j) = w(i).*data(i,j);
        nu(j) = sum(a(:,j))/sum(w); % The weighted mean value of variables.
    
    end
   
end

% C.2. Calculating the weighted corvariance matrix

for j = 1:j_dim
       
    for k = 1:j_dim
            
        c(j,k) = sum(w(:,1).*(data(:,j)-nu(j)).*(data(:,k)-nu(k)));
        d(j,k) = sum(w(:,1));
            
        s(j,k) = c(j,k)./d(j,k);
            
    end     
    
end

xlswrite('SWS.xls',s);

%% D. Calculating the Spatially Weighted Correlation Matrix
 
for j = 1:j_dim
       
    for k = 1:j_dim
            
        e(j,k) = sqrt(sum(w(:,1).*((data(:,j)-nu(j)).*(data(:,j)-nu(j)))));
        f(j,k) = sqrt(sum(w(:,1).*((data(:,k)-nu(k)).*(data(:,k)-nu(k)))));
            
        m(j,k) = e(j,k).*f(j,k);
        r(j,k) = c(j,k)./m(j,k);
            
    end     
    
end

xlswrite('SWR.xls',r);
