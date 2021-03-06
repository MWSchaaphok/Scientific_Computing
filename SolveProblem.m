%% Scientific Computing  
% Group assignment for course Scientific Computing 
% By: Nerine Usman & Marianne Schaaphok
% Date: 25-11-2018

function [u,u_ex,err,tF,tS, fill_ratio, resid,rrf,m] = SolveProblem(p,dimension,iter, solver,redsc,m_max)
%% Parameters

%dimension = 3; 
tol = 1e-4;
maxnorm = size(1,9);
%p = 2; 
n = 2^p; 
h = 1/n;

%% Construct matrices A
 
% Help matrices for construction of A 
H_1 = spdiags([h^2; zeros(n-1,1); h^2],0,n+1,n+1);
D_1 = spdiags([0; ones(n-1,1);0], 0, n+1,n+1);
D_2 = kron(D_1,D_1);
if dimension == 3
    D_3 = kron(D_1,D_2);
end 
T_1 = spdiags(-1*[0;ones(n-2,1); 0; 0],-1,n+1,n+1) ...
    + spdiags(-1*[0;0;ones(n-2,1);0],1,n+1,n+1);
I_1 = spdiags(ones((n+1),1),0,(n+1),(n+1));
I_2 = spdiags(ones((n+1)^2,1),0,(n+1)^2,(n+1)^2);


% Boundary neighbours of interior points
B_1 = sparse(zeros(n+1)); 
B_1(2,1) = 1; 
B_1(n,n+1) = 1; 
B_2 = kron(D_1,B_1) + kron(B_1,D_1);

if dimension == 3 
    B_3 = kron(D_1,B_2) + kron(B_1,D_2);
end 

% 1D matrix A 
A_1 = H_1+2*D_1 + T_1;

% 2D matrix A
A_2 = kron(H_1,I_1) + kron(D_1,A_1+2*D_1) + kron(T_1,D_1);

% 3D matrix A
if dimension == 3
    A_3 = kron(H_1,I_2) + kron(D_1,A_2+2*kron(D_1,D_1)) + kron(T_1,D_2);
end 

A_1 = 1/h^2 * A_1;
A_2 = 1/h^2 * A_2; 
if dimension == 3
    A_3 = 1/h^2 * A_3; 
end 
%imagesc(A_2)

%% Construct vector f 2D
if (dimension == 2)
    % Construct mesh
    x = 0:h:1; 
    y = 0:h:1; 
    [X,Y] = meshgrid(x,y); 


    % Construct f for internal points 2D problem
    F2 = f2(X,Y); 
    F2 = reshape(F2,[(n+1)^2,1]);
    f_int2 = D_2*F2; 

    % Add f for boundary points
    U2 = exact(X,Y,1);
    U2 = reshape(U2, [(n+1)^2,1]);
    f_boun2 = (I_2-D_2)*U2; 

    % Correct f for symmetric A
    f_nb2 = B_2*(U2/h^2);

    f_2 = f_int2 + f_boun2 + f_nb2; 
    
    A = A_2; 
    f = f_2; 
    u_ex = U2;
%% Construct vector f 3D
elseif (dimension == 3)
    I_3 = spdiags(ones((n+1)^3,1),0,(n+1)^3,(n+1)^3);
    
    % Construct mesh
    x = 0:h:1; 
    y = 0:h:1; 
    z = 0:h:1;
    [X,Y,Z] = meshgrid(x,y,z); 


    % Construct f for internal points 2D problem
    F3 = f3(X,Y,Z); 
    F3 = reshape(F3,[(n+1)^3,1]);
    f_int3 = D_3*F3; 

    % Add f for boundary points
    U3 = exact(X,Y,Z);
    U3 = reshape(U3, [(n+1)^3,1]);
    f_boun3 = (I_3-D_3)*U3; 

    % Correct f for symmetric A
    f_nb3 = B_3*(U3/h^2);

    f_3 = f_int3 + f_boun3 + f_nb3; 
    
    A = A_3; 
    f = f_3; 
    u_ex = U3;
else 
    fprintf('Please choose dimension 2 or 3')
end

%% Solve the system
% Use of bandwith reduction scheme
if redsc == 1 
     s = symamd(A);
%      figure;
%      spy(A); 
%      title('A before matrix reordering')
     A = A(s,s); 
%      figure
%      spy(A)
%      title('A after matrix reordering')
     f = f(s);
end 
 
 % Define variables 
 R = sparse(size(A));
 %R = zeros(size(A));
 u = zeros(size(f));
 resid = zeros(m_max,1); 
 normR = zeros(m_max,1);
 rrf = zeros(5,1); 
 tF = 0; 
 tS = 0; 
 r = 1;
 nf = norm(f);
 m=1;
 
 % Use given solver 
 if strcmp(solver,'Cholesky')
    tic; 
    R = chol(A,'lower'); 
    tF = toc;

%     Cholesky Decomposition
%     for k = 1:size(A,2)
%         j = max(k-(n+1),1);
%         A(k,k) = sqrt(A(k,k) - A(k,j:k-1)*A(k,j:k-1)');
%         for i = k+1:min(k+n+1,size(A,2))
%             l = max(i-(n+1),1);
%            A(i,k) = 1/(A(k,k))*(A(i,k) - A(i, l:k-1)*A(k,l:k-1)');
%         end
%     end 
%     R = tril(A);
     
    
    % Direct solving algorithm
    tic; 
    y = R\f; 
    u = R'\y;
    tS = toc; 
    
    if redsc == 1 
    u(s) = u; 
    end 

elseif strcmp(solver,'SSOR')
     omega = 1.5; 
     resid(1) = 1; 
     u = zeros(size(u_ex));
     tic;
     while m<m_max && norm(r)/nf>10^-10
         %tic;
         for i = 1:length(u)
            sigma = u(i); 
            %u(i) = (f(i) - A(i,1:i-1)*u(1:i-1) - A(i,i+1)*u(i+1,n))/A(i,i);
            u(i) = 0; 
            u(i) = (f(i)-A(i,:)*u)/A(i,i);
            u(i) = (1-omega)*sigma + omega*u(i);
         end

         for i = length(u):-1:1
            sigma = u(i); 
            u(i) = 0; 
            u(i) = (f(i)-A(i,:)*u)/A(i,i);
            u(i) = (1-omega)*sigma + omega*u(i);
         end
         %tS = toc;
         r = f-A*u;
         normR(m) = norm(r); 
         resid(m) = normR(m)/nf;
         m=m+1; 
         
     end 
     
     tS = toc; 
     rrf = normR(max(m-5,1):m-1)./normR(max(m-6,1):m-2);
     
  
 elseif strcmp(solver,'PCG')
     omega = 1.5; 
     % Initialize
     u = zeros(size(u_ex));
     
     % Decompose A 
     D = diag(diag(A));
     E = D - tril(A); 
     F = E';
     
     % Preconditioning matrix
     M = (D-omega*E)*spdiags(1./diag(D),0,size(D,1),size(D,2))*(D-omega*F)/(omega*(2-omega));
     
     % Algorithm
     tic; 
     r = f; 
     
     while m<=m_max && norm(r)/nf>10^-10
        z = M\r;
        if  m==1
            p = z;  
            current = r'*z;
        else
            current = r'*z;
            beta = current/prev;
            p = z + beta*p;
        end
        temp = A*p; 
        alpha = current/(p'*temp);
        u = u + alpha*p;
        r = r - alpha*temp;
        prev = current;
        
        resid(m) = norm(r)/nf;
        m = m+1;
        %tS = toc;        
     end 
     tS = toc; 
     
 end 

%  figure 
%  spy(R)
%  title('Cholesky matrix after reordering')
%  hold off; 
 fill_ratio = nnz(R)/nnz(A);
 err = norm(u-u_ex, 'inf');

% if dimension==2
%     u_pl = reshape(u,[(n+1),(n+1)]);
%     u_ex1 = reshape(U2,[n+1,n+1]);
% 
%     figure; 
%     surf(X,Y,u_pl)
%     hold on;
%     surf(X,Y,u_ex1)
%     hold off; 
% end 


%% Functions
function [f] = f2(x,y)
    f = (x.^2 + y.^2).*sin(x.*y);
end 

function [f] = f3(x,y,z)
    f = ((x.*y).^2 + (y.*z).^2 + (x.*z).^2).*sin(x.*y.*z);
end 

function [ex] = exact(x,y,z)
    ex = sin(x.*y.*z); 
end
end
