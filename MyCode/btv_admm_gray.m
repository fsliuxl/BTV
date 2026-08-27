function xrec = btv_admm_gray(y,A,opts)
% solve min_x 1/2||y-Ax||_2^2 + lambda * ||Dx||_{2,1}
% y observed signal
% A measurement operator
mu      = 1e-2;
rho     = 1.1;
maxiter = 1000;
tol     = 1e-8;
lambda  = 5e-1;
if isfield(opts,'mu'); mu = opts.mu; end
if isfield(opts,'rho'); rho = opts.rho; end
if isfield(opts,'maxiter'); maxiter = opts.maxiter; end
if isfield(opts,'tol'); tol = opts.tol; end
if isfield(opts,'lambda'); lambda = opts.lambda; end
if isfield(opts,'block_size'); block_size = opts.block_size; end

Block_Size = [block_size,1];
N = opts.N;
% Initialization
x0 = A'*y;
z  = zeros([N,2]);
Dualz = z;

for iter = 1:maxiter
    % Update x
    temp = A'*y;
    rhs = temp(:) + mu*diffT2(z - Dualz/mu,N);
    x1   = myPCG_btv(x0(:),A,rhs(:),N,mu);
    x1   = reshape(x1,N);
    %     lhs  = A'*A + mu*(Dif'*Dif);
    %     x1   = lhs\rhs;
    % Update z
    z   = shrink_block(reshape(diff2(x1,N),[N,2]) + Dualz/mu, lambda/mu, Block_Size);
    z   = reshape(z,[N,2]);
    % Update Dualz
    Dualz = Dualz + mu*(reshape(diff2(x1,N),[N,2]) - z);
    mu    = min(rho*mu,1e8);
    if stopcretia(x0,x1) < tol; break; end
    x0 = x1;
end
xrec = reshape(x1,N);
end

%% block shrink
function z = shrink_block(x, r, Block_Size)
x = x(:);
process_func = @(block_struct) sqrt(sum(block_struct.data.^2, 'all'));
result = blockproc(x, Block_Size, process_func);
t_value = kron(result, ones(Block_Size));
t_value = t_value(1:size(x,1), 1:size(x,2), :);
z = max(0, 1 - r./(t_value + 1e-16)) .* x;
end
%% stop cretia
function y = stopcretia(x1,x2)
x1 = x1(:); x2 = x2(:);
y = norm(x1-x2)/max(norm(x1),1);
end


