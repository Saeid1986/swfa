function fabycovcor(X)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%     FACTOR ANALYSIS     %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%       BY COVARIANCE       %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%       & CORRELATION       %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%          MATRIX         %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Author: Saeid ESMAEILOGHLI
% February 26, 2019
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
% This program was designed to solve the factor analysis (FA) model by
% using principal component solution and based on covariance or correlation
% matrices. The possibility of adopting retainable factors in the FA model
% was provided according to the latent root criterion, while a varimax
% factor rotation is executable if the user chooses.
% The main reference is:
% Rencher, A. C. "Methods of Multivariate Analysis." 2nd. edition,
% New-Jersey: John Wiley & Sons. Chapter 13, (2002): 408-450.
%
% This program was modified after the code created by:
% A. Trujillo-Ortiz, R. Hernandez-Walls, A. Castro-Perez and M.
% Rodriguez-Ceja, "ANFACTPCWOD: Factor Analysis by the Principal Components
% Method Without Data. A MATLAB file. [WWW document]".
% URL http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?
% objectId=10602).
%
% Syntax: function fabycovcor(X)
%      
% The input data include:
%
% X: a correlation matrix. One can also input a covariance matrix, that
% through the standardization of any input matrix it assure a correlation
% matrix of it and needed for the procedure.
%
% The program provides the following outputs:
%
% Complete factor analysis results such as:
% - Table of primary extracted factors;
% - Table of unrotated principal components of the factor analysis;
% - Proportion of total (standardized) sample variance;
% - Table of cumulative proportion of total (standardized) sample variance;
% 
% as well as optional results such as:
% - Table of varimax rotated factor loadings;
% - Residual Matrix.
%
% To cite this file, this would be an appropriate format:
% Esmaeiloghli, S., Tabatabaei, S.H., Carranza, E.J.M., Hosseini, S.,
% Deville, Y. "Spatially-weighted factor analysis for extraction of
% source-oriented mineralization feature in 3D coordinates of surface
% geochemical signal." Journal of Geochemical Exploration (2020).


%% Input Covariance (Correlation) Matrix

X = [0.647885333	0.042610651	-0.087393774	-0.138437935	-0.041927208	0.065960929	0.07755487	-0.135361399	-0.059996214	0.096511842	-0.061946771	-0.080724486	0.086176433	-0.054389586	-0.134040138	-0.054565316	-0.09403313	-0.060752916	-0.013131184
0.042610651	0.747702431	-0.110392691	-0.406335398	-0.105822815	-0.150329321	0.013072184	-0.281437636	-0.192452653	0.33235212	-0.102203218	0.250006343	0.445646088	-0.077265773	-0.12207287	0.094054549	-0.084055367	-0.112008211	-0.181068413
-0.087393774	-0.110392691	0.124968992	0.124986201	0.083219706	-0.104525698	-0.260353346	0.077596658	0.129759389	-0.069640876	-0.057580575	-0.059695858	-0.11091334	0.0221461	0.146218625	0.033746037	0.075400954	0.004251749	0.038201746
-0.138437935	-0.406335398	0.124986201	1.004680001	0.195301363	-0.042188304	-0.824394394	0.61660967	0.429423991	-0.406705381	-0.111920933	-0.504370407	-0.376206124	-0.128450875	0.4056518	-0.339507166	0.177870863	-0.157509442	0.48150247
-0.041927208	-0.105822815	0.083219706	0.195301363	0.313764293	-0.201417643	-0.431957835	0.060300393	0.162078114	-0.100200846	-0.054825015	-0.106592967	-0.0619293	-0.007343774	0.149602625	0.007259423	0.090061315	-0.026266313	0.076696484
0.065960929	-0.150329321	-0.104525698	-0.042188304	-0.201417643	0.789061847	0.777055451	0.024754706	-0.262205737	-0.191078299	0.147482938	-0.230748078	-0.094577956	0.026145628	-0.334306163	-0.087443608	-0.149692762	0.03632825	-0.01827618
0.07755487	0.013072184	-0.260353346	-0.824394394	-0.431957835	0.777055451	2.190943015	-0.451320937	-0.778599941	-0.068022917	0.375298673	0.214831816	0.147133709	0.239544561	-0.864536768	0.143302783	-0.341843837	0.320874578	-0.478581663
-0.135361399	-0.281437636	0.077596658	0.61660967	0.060300393	0.024754706	-0.451320937	0.530604832	0.238421339	-0.268736076	-0.077961852	-0.346185665	-0.256149487	-0.088168122	0.244500815	-0.211899579	0.114874116	-0.108618259	0.318176483
-0.059996214	-0.192452653	0.129759389	0.429423991	0.162078114	-0.262205737	-0.778599941	0.238421339	0.450705682	-0.089670775	-0.132704207	-0.205194623	-0.214165217	-0.069417906	0.402736479	-0.106833306	0.142588724	-0.065782758	0.221309619
0.096511842	0.33235212	-0.069640876	-0.406705381	-0.100200846	-0.191078299	-0.068022917	-0.268736076	-0.089670775	0.720899929	-0.141662663	0.198999357	0.270557631	-0.111300806	-0.021229825	0.128705307	-0.10779009	-0.082242155	-0.089745478
-0.061946771	-0.102203218	-0.057580575	-0.111920933	-0.054825015	0.147482938	0.375298673	-0.077961852	-0.132704207	-0.141662663	0.536496363	-0.025433246	-0.075364995	0.056043541	-0.171108821	-0.031388649	-0.018375632	0.058486605	-0.111331542
-0.080724486	0.250006343	-0.059695858	-0.504370407	-0.106592967	-0.230748078	0.214831816	-0.346185665	-0.205194623	0.198999357	-0.025433246	0.919444349	0.195748814	0.034300171	-0.088840731	0.149128072	-0.033457939	0.07231413	-0.35352905
0.086176433	0.445646088	-0.11091334	-0.376206124	-0.0619293	-0.094577956	0.147133709	-0.256149487	-0.214165217	0.270557631	-0.075364995	0.195748814	0.534936588	-0.055353333	-0.178867855	0.085176759	-0.114377941	-0.062225664	-0.165244808
-0.054389586	-0.077265773	0.0221461	-0.128450875	-0.007343774	0.026145628	0.239544561	-0.088168122	-0.069417906	-0.111300806	0.056043541	0.034300171	-0.055353333	0.183008422	-0.095944569	0.105897144	0.029127227	0.120354749	-0.128932801
-0.134040138	-0.12207287	0.146218625	0.4056518	0.149602625	-0.334306163	-0.864536768	0.244500815	0.402736479	-0.021229825	-0.171108821	-0.088840731	-0.178867855	-0.095944569	0.51669947	-0.091166552	0.167588899	-0.104165681	0.173281259
-0.054565316	0.094054549	0.033746037	-0.339507166	0.007259423	-0.087443608	0.143302783	-0.211899579	-0.106833306	0.128705307	-0.031388649	0.149128072	0.085176759	0.105897144	-0.091166552	0.282622294	-0.016024498	0.071670926	-0.16273462
-0.09403313	-0.084055367	0.075400954	0.177870863	0.090061315	-0.149692762	-0.341843837	0.114874116	0.142588724	-0.10779009	-0.018375632	-0.033457939	-0.114377941	0.029127227	0.167588899	-0.016024498	0.163642352	-0.029250395	0.027747142
-0.060752916	-0.112008211	0.004251749	-0.157509442	-0.026266313	0.03632825	0.320874578	-0.108618259	-0.065782758	-0.082242155	0.058486605	0.07231413	-0.062225664	0.120354749	-0.104165681	0.071670926	-0.029250395	0.253695441	-0.129154634
-0.013131184	-0.181068413	0.038201746	0.48150247	0.076696484	-0.01827618	-0.478581663	0.318176483	0.221309619	-0.089745478	-0.111331542	-0.35352905	-0.165244808	-0.128932801	0.173281259	-0.16273462	0.027747142	-0.129154634	0.49481517
];

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%       MAIN       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%      SCRIPT      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%

error(nargchk(1,1,nargin))

[r c] = size(X);

if r ~= c, 
    error('Input covariance/correlation matrix it is not square!');
elseif ~all(all(X==X')),
    error('Input covariance/correlation matrix it is not symmetric!');
elseif any(diag(X) <= 0),
    error('The covariance/correlation matrix must be positive definite!');
else
end

dX = diag(X);
D = 1./sqrt(dX);
R = X.*(D*D');
p = r;
[A L A] = svd(R,0);
l = diag(L);
P = (l/sum(l))*100;
CP = cumsum(P);
ft = [1:p]';

disp(' ')
disp('Table of extracted primary factors:')
fprintf('-------------------------------------------------------------\n');
disp('                            Percent of       Cummulative       ');
disp(' Factors     Eigenvalue      Variance     Percent of Variance  '); 
fprintf('-------------------------------------------------------------\n');
fprintf('    %d      %10.4f      %10.4f        %10.4f\n',[ft,l,P,CP].');
fprintf('-------------------------------------------------------------\n');

B = A*sqrt(L);
f = L >= 1.0;
m = sum(sum(f));

disp('  ');
fprintf('According to the latent root citerion, the suggested number for retainable factors is %.i\n', m);
disp('  ');
ask = input('Do you want to work with this number of factors? (y/n): ','s');
if~strcmp(ask,'y'),
   m = input('Enter the number of factors you need: ');
   while (m > p)
       disp(' ');
       fprintf('Error: The number of factors requested is too large for the number of the observed variables. It must be equal or lesser than %.i\n', p);
       disp(' ');
       m = input('Enter the number of factors you need: ');
   end
end

F = B(:,1:m);
pt = sum(F.^2)/p;
pp = [pt sum(pt)];

C = F.^2;
C = sum(C,2);
Factors = [F C];

disp(' ')
disp('Table of unrotated factor model:')
disp('-----------------------------------------------------------------------------');
Factors
fprintf('-----------------------------------------------------------------------------\n');
fprintf('On factors, Factor 1 = column 1 and so forth to %.i\n', m);
disp('Last column is the communality');
fprintf('On variates, Variate 1 = first row and so forth to %.i\n', p);
disp(' ')
disp('Proportion of total (standardized) sample variance:');
disp('-----------------------------------------------------------------------------');
pp
fprintf('-----------------------------------------------------------------------------\n');

sv = 1 - C;
Re = F*F' + diag(sv);
Rm = R - Re;

disp('  ');
rt = input('Do you want to execute a varimax rotation on the loading matrix? (y/n): ','s');
if rt == 'y'
    loadings = F;
    b = loadings;
    [n,nf] = size(loadings);
    hjsq = diag(loadings*loadings');
    hj =sqrt(hjsq);

    
    for iter = 1:10; 
        
        for i = 1:nf-1,
            jl = i+1;
            for j = jl:nf,
                xj = loadings(:,i)./hj;
                yj = loadings(:,j)./hj;
                uj = xj.*xj - yj.*yj;
                vj = 2*xj.*yj;
                A = sum(uj);
                B = sum(vj);
                C = uj'*uj - vj'*vj;
                D = 2*uj'*vj;
                num = D - 2*A*B/n;
                den = C - (A^2 - B^2)/n;
                tan4p = num/den;
                phi = atan2(num,den)/4;
                angle = phi*180/pi;
                if abs(phi) > 0.00001;
                    Xj = cos(phi)*xj + sin(phi)*yj;
                    Yj = -sin(phi)*xj + cos(phi)*yj;
                    bj1 = Xj.*hj;
                    bj2 = Yj.*hj;
                    b(:,i) = bj1;
                    b(:,j) = bj2;
                    loadings(:,i) = b(:,i);
                    loadings(:,j) = b(:,j);
                end
            end
        end
        loadings = b;
    end
    F = loadings;
    pt = sum(F.^2)/p;
    pp = [pt sum(pt)];
    
    
    C = F.^2;
    C = sum(C,2);
    Factors = [F C];
    
    disp(' ')
    disp('Table of varimax rotated factor model:')
    disp('-----------------------------------------------------------------------------');
    Factors
    fprintf('-----------------------------------------------------------------------------\n');
    fprintf('On factors, Factor 1 = column 1 and so forth to %.i\n', m);
    disp('Last column is the communality');
    fprintf('On variates, Variate 1 = first row and so forth to %.i\n', p);
    disp(' ')
    disp('Table of cumulative proportion of total (standardized) sample variance:');
    disp('-----------------------------------------------------------------------------');
    pp
    fprintf('-----------------------------------------------------------------------------\n');
    
    sv = 1 - C;
    Re = F*F' + diag(sv);
    Rm = R - Re;
    
else
end

disp(' ');
rm = input('Do you need to output the residual matrix? (y/n): ','s');
disp(' ');
if rm == 'y',
    disp('Residual matrix:');
    Rm
else
end

return

