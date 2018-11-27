%% Scientific Computing  
% Group assignment for course Scientific Computing 
% By: Nerine Usman & Marianne Schaaphok
% Date: 25-11-2018

clear;
%% Parameters
p=2; 
n = 2^p; 
h = 1/n; 
dimension = 2; 
tol = 1e-4;

T_1 = -1*diag(ones(n,1),-1) + 4*eye(n+1) -1*diag(ones(n,1),1);

%% Matrix and vector construction in 2D
if dimension == 2 
    % Construct matrix A without elimination 
    A(1:(n+1),01:(n+1)) = T_1; 
    A(1:n+1, n+2: 2*n+2) = -eye(n+1); 

    for i = 2:n 
        A((i-1)*(n+1)+1:i*(n+1),(i-2)*(n+1)+1:(i-1)*(n+1)) = -eye(n+1); 
        A((i-1)*(n+1)+1:i*(n+1),(i-1)*(n+1)+1:i*(n+1)) = T_1; 
        A((i-1)*(n+1)+1:i*(n+1),(i)*(n+1)+1:(i+1)*(n+1)) = -eye(n+1);
    end 

    A((n)*(n+1)+1:(n+1)*(n+1),(n-1)*(n+1)+1:(n)*(n+1)) = -eye(n+1); 
    A((n)*(n+1)+1:(n+1)*(n+1),(n)*(n+1)+1:(n+1)*(n+1)) = T_1; 

    % Construct vector f  without elimination 
    f = ones((n+1)^2,1);
    for i = 1:(n+1)
       for j = 1:(n+1)
            f((i-1)*(n+1)+j) = sin(((i-1)*h)*((j-1)*h));
       end 
    end

%% Matrix and vector construction in 3D
elseif dimension == 3 
    
    fprintf('Not yet implemented\n')
    
else 
    fprintf('Please choose dimension 2 or 3\n')
end 

%% Direct Solver
% Compute Cholesky decomposition
R = chol(A); 

% Firect solving algorithm for norm(r)<tolerance 
u = 0;
r = f;
i = 0;
while norm(r)>tol
    y = R\r; 
    du = R'\y; 
    u = u+du; 
    r = f-A*u; 
    i = i+1; 
    if i>50
       fprintf('Could not converge within 50 iteration\n') 
       break
    end
end
