%% Scientific Computing  
% Group assignment for course Scientific Computing 
% By: Nerine Usman & Marianne Schaaphok
% Date: 25-11-2018

clear;
%% Parameters

dimension = 3; 
tol = 1e-4;
maxnorm = size(1,9);
p = 2; 
n = 2^p; 
h = 1/n;

%% Construct matrices A
 
% Help matrices for construction of A 
H_1 = spdiags([h^2; zeros(n-1,1); h^2],0,n+1,n+1);
D_1 = spdiags([0; ones(n-1,1);0], 0, n+1,n+1);
D_2 = kron(D_1,D_1);
D_3 = kron(D_1,D_2); 
T_1 = spdiags(-1*[0;ones(n-2,1); 0; 0],-1,n+1,n+1) ...
    + spdiags(-1*[0;0;ones(n-2,1);0],1,n+1,n+1); 
I_1 = eye(n+1);
I_2 = eye((n+1)^2);
I_3 = eye((n+1)^3);

% Boundary neighbours of interior points
B_1 = zeros(n+1); 
B_1(2,1) = 1; 
B_1(n,n+1) = 1; 
B_2 = kron(D_1,B_1) + kron(B_1,D_1);
B_3 = kron(D_1,B_2) + kron(B_1,D_2);

% 1D matrix A 
A_1 = H_1+2*D_1 + T_1;

% 2D matrix A
A_2 = kron(H_1,I_1) + kron(D_1,A_1+2*D_1) + kron(T_1,D_1);

% 3D matrix A

A_3 = kron(H_1,I_2) + kron(D_1,A_2+2*kron(D_1,D_1)) + kron(T_1,D_2);

A_1 = 1/h^2 * A_1;
A_2 = 1/h^2 * A_2; 
A_3 = 1/h^2 * A_3; 



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
    R = chol(A); 

    % Firect solving algorithm for norm(r)<tolerance 
    r = f;
    u = zeros(size(f)); 
    i = 0;
    while norm(r)>tol
        y = R\r; 
        du = R'\y; 
        u = u+du; 
        r = f-A*u; 
        i = i+1; 
        if i>1000
           fprintf('Could not converge within 50 iteration\n') 
           break
        end
    end
err = norm(u-u_ex)

if dimension==2
    u_pl = reshape(u,[(n+1),(n+1)]);
    u_ex1 = reshape(U2,[n+1,n+1]);

    figure; 
    surf(X,Y,u_pl)
    hold on;
    surf(X,Y,u_ex1)
end

figure;
imagesc(A)

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
