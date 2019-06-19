%%%% 2D MINIMUM-VOLUME TOPOLOGY OPTIMIZATION CODE %%%%
%%%% WITH APPROXIMATE REANALYSIS %%%%
%%%% NO MATRIX ASSEMBLY FOR REANALYSIS %%%%
% General parts of the code are ased on 88-line code (Andreassen et al.)
% This code is freely distributed as complementary material to the article:
% Amir O, "Revisiting approximate reanalysis in topology optimization: on
% the advantages of recycled preconditioning in a minimum weight procedure"
% http://link.springer.com/article/10.1007/s00158-014-1098-7
% To reproduce the testcase in section 3.2 of the article, run:
% minVreMF(200,100,3,2,2,136.1,'trial')
function minVreMF(nelx,nely,penal,rmin,ft,compconst,filename)
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
ndof = 2*(nelx+1)*(nely+1);
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% Define loads and supports (cantilever with two point loads)
F = sparse(2*(nelx+1)*(nely+1)-1:2*(nelx+1)*(nely+1),1,[-0.5 -1],2*(nely+1)*(nelx+1),1);
fixeddofs = 1:2*(nely+1);
U = zeros(2*(nely+1)*(nelx+1),1);
% Null space elimination of supports
N = ones(ndof,1); N(fixeddofs) = 0; Null = spdiags(N,0,ndof,ndof);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
x = ones(nely,nelx);
xPhys = x;
loop = 0;
change = 1;
lam = 1e10;
pcgr = 1;
count = 9;
Stats = zeros(200,4);
%% START ITERATION
while (change > 1e-3 && loop < 200)
    loop = loop + 1;
    %% FE-ANALYSIS
    if (count == 9 || pcgr > 1e-2)
        sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
        K = sparse(iK,jK,sK); K = (K+K')/2;
        K = Null'*K*Null - (Null-speye(ndof,ndof));
        [Ufac,~,S] = chol(K);
        Lfac = Ufac';ST = S';
        U = S * (Ufac \ (Lfac \ (ST*F)));
        pcgi = 0; pcgr = 0; count = 0;
    else
        [pcgi,pcgr,U] = mypcgMatrixFree(F,U,Lfac,Ufac,S,ST,1e-3,5,...
            xPhys,KE,Emin,E0,penal,edofMat,Null,1);
        count = count + 1;
    end
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    comp = full(sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce)));
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dv = ones(nely,nelx)/nelx/nely;
    if (loop == 1 && nargin < 7)
        compconst = 2*comp;
    end
    f0val = full(mean(xPhys(:)));
    fval = comp/compconst - 1;
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
    end
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    l1 = 0; l2 = lam*2; move = 0.05; i = 0;
    while (l2-l1)/(l1+l2) > 1e-6
        % Check if uniform reduction violates linearized constraint
        if (comp - compconst - sum(sum(dc))*move < 0)
            xnew = x - move;
            lam = 1e10;
            break;
        end
        % Non-uniform design change
        i = i + 1;
        lam = 0.5*(l2+l1);
        xnew = max(1e-10,max(x-move,min(1,min(x+move,x.*sqrt(-lam*dc./dv)))));
        if comp - compconst + dc(:)'*((xnew(:)-x(:)).*(x(:)./xnew(:))) > 0, l1 = lam; else l2 = lam; end
    end
    change = full(max(abs(xnew(:)-x(:))));
    x = xnew;
    if ft == 1
        xPhys = xnew;
    elseif ft == 2
        xPhys(:) = (H*xnew(:))./Hs;
    end
     %% PRINT RESULTS
    fprintf(' It.:%5i Con.:%11.3e Vol.:%11.3e Lam.: %11.3e ch.:%7.3f innerit: %3i pcgi: %3i pcgr: %11.3e\n',...
        loop,fval,f0val,lam,change,i,pcgi,pcgr);
    %% PLOT DENSITIES
    colormap(jet); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
    %% Keep stats
    Stats(loop,1) = pcgi;
    Stats(loop,2) = pcgr;
    Stats(loop,3) = f0val;
    Stats(loop,4) = fval;
end
%% Save results
save(filename,'f0val','xPhys','Stats');
end
%% Matrix-free PCG with standard & alternative convergence criteria
function [i,relres,u] = mypcgMatrixFree(b,u,Lfac,Ufac,S,ST,tol,maxiter,...
    xPhys,KE,Emin,E0,penal,edofMat,Null,toltype)
xPhys8 = repmat(xPhys(:),1,8);
Evals = Emin + xPhys8.^penal*(E0-Emin);
Auvals = u(edofMat)*KE;
Auvals = Evals.*Auvals;
Au = Null*sparse(edofMat,1,Auvals);
r = b - Au;
res0 = norm(b);
uAu1 = u'*Au;
for i = 1:1e6
    z = S * (Ufac \ (Lfac \ (ST*r)));
    rz = r'*z;
    if i==1
        p = z;
    else
        beta = rz/rzold;
        p = beta*p + z;
    end
    Apvals = p(edofMat)*KE;
    Apvals = Evals.*Apvals;
    Ap = Null*sparse(edofMat,1,Apvals);
    pAp = p'*Ap;
    alpha = rz/pAp;
    u = u + alpha*p;
    Au = Au + alpha*Ap;
    r = r - alpha*Ap;
    rzold = rz;
    relres = norm(r)/res0;
    uAu2 = u'*Au;
    altres = abs(uAu2-uAu1)/uAu2;
    if ((toltype==1 && relres < tol) || i >= maxiter)
        break
    elseif ((toltype==2 && altres < tol) || i >= maxiter)
        break
    end
    uAu1 = uAu2;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Oded Amir,                                     %
% Faculty of Civil & Environmental Engineering,                           %  
% Technion - Israel Institute of Technology.                              %
%                                                                         %
% Contact: odedamir@technion.ac.il                                     %
%                                                                         %
% Details are discussed in the paper:                                     %
% "Revisiting approximate reanalysis in topology optimization: on the     %
% advantages of recycled preconditioning in a minimum weight procedure",  %
% http://link.springer.com/article/10.1007/s00158-014-1098-7              %
%                                                                         %
% Disclaimer:                                                             %
% The author reserves all rights but does not guarantee that the code is  %
% free from errors. Furthermore, the author shall not be liable in any    %
% event caused by the use of the program.                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

