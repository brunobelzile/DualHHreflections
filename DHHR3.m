function [x,x_o,A,A_o,b,b_o,H,H_o]=DHHR3(A,A_o,b,b_o)
%% Dual Householder Reflections
% Bruno Belzile
% CIM, McGill University
% June 2019

% Input: matrices and vectors representing the system of linear equations
% Output: least-square solution, transformed matrices before
% back substitution

% Please refer to the following paper:
% Belzile, B. and Angeles, J., 2019, 
% “Reflections Over the Dual Ring---Applications to Kinematic Synthesis,”
% ASME Journal of Mechanical Design, Vol. 141, Issue 7, 
% pp. 072302-1--072302-9, DOI: 10.1115/1.4043204.

%% Preallocation of the variables
[m,n] = size(A);
H = eye(m); H_o = zeros(m);
x = zeros(n,1); x_o = zeros(n,1);

%% Error detection
if ~and(size(A_o,1)==m,size(A_o,2)==n) 
    error('Matrices A and A_o must have the same dimensions.');
elseif size(b,1)~=m 
    error('The numbers of rows in A and b must be equal.');
elseif size(b_o,1)~=m
    error('Vectors b and b_o must have the same length.');
elseif n>m
    error('Underdetermined systems are not yet supported.');
elseif rank(A)<min(size(A))
    error('Rank-deficient systems (primal part) are not yet supported.');
end

if n<m
    A = transpose(A); A_o = transpose(A_o);
end
    
%% Computation of the n individual Householder reflections (matrices)
for i=1:n
    a_ii = A(i,i); a_oii = A_o(i,i);
    a_i = A(:,i); a_oi = A_o(:,i);
    alpha = norm(a_i(i:end)); % primal
    if a_ii<0; alpha = -alpha; end
    alpha_o = (transpose(a_i(i:end))*a_oi(i:end))/alpha; % dual
    
    % primal
    u = zeros(m,1); 
    u(i,1) = a_i(i) + alpha;
    u(i+1:m,1) = a_i(i+1:m);
    
    % dual
    u_o = zeros(m,1);
    u_o(i,1) = a_oi(i) + alpha_o;
    u_o(i+1:m,1) = a_oi(i+1:m);
    
    lambda = norm(u)^2/2; % primal
    lambda_o = 2*alpha_o*alpha+a_ii*alpha_o+a_oii*alpha; % dual
    
    X_i = u*transpose(u)/lambda; % primal
    X_oi = ((-u*transpose(u)*lambda_o+(u*transpose(u_o)+...
        u_o*transpose(u))*lambda)/lambda^2); % dual
    
    H_i = eye(m)-X_i; % primal
    H_oi = -X_oi; % dual
    
    H_o = H_i*H_o+H_oi*H; % dual
    H = H_i*H; % primal
    
    A_o = H_i*A_o+H_oi*A; % dual
    A = H_i*A; % ptimal
end

%% Application of the Householder reflections
b_o = H*b_o+H_o*b; % dual
b = H*b; % primal

%% Back substitution
for i=n:-1:1
    v = A(i,:)*x(:,1); % primal
    v_o = A(i,:)*x_o(:,1)+A_o(i,:)*x(:,1); % dual
    x(i,1) = (b(i)-v)/A(i,i); % primal
    x_o(i,1) = (b_o(i)-A_o(i,i)*x(i,1)-v_o)/A(i,i); % dual
end
end