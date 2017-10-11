%% test tridiag mex
%% check value against \ and ir_apply_tridiag

nrep = 10;
for jj = 1:nrep
    rng(jj);
    N = 200;
    M = 100;
    scale = 10;
    d = scale*randn(N,M);
    d = d + 1i*scale*rand(N,M);
    a = scale*rand(N-1,1)-scale/2;
    b = scale*rand(N,1)-scale/2;
    c = scale*rand(N-1,1)-scale/2;
    
    d = single(d);
    a = single(a);
    b = single(b);
    c = single(c);
    
    T = diag(a,-1) + diag(b) + diag(c,1);
    x0 = T\d;
    
    x1 = ir_apply_tridiag_inv(a, b, c, d);
    
    % compile as needed
    % mex -O CFLAGS="\$CFLAGS -std=c99" -I./def/ tridiag_inv_mex.c
    
    
    try
        x2 = tridiag_inv_mex(a, b, c, d);
    catch
        display('tridiag_inv_mex.c failed');
    end
    if norm(x0-x2)/N > 1e-3
        %     keyboard;
    end
    norm(x0-x2)/N
end

% norm(x0-x2)
norm(x1-x2)
%equivs(x0, x2)
%equivs(x1, x2)

%% speed test with \ and ir_apply_tridiag

rng(0);
nrep = 4;
N = 200;
M = 100;
scale = 10;
d = scale*randn(N,M);
d = single(d + 1i*scale*rand(N,M));
a = single(scale*rand(N-1,1)-scale/2);
b = single(scale*rand(N,1)-scale/2);
c = single(scale*rand(N-1,1)-scale/2);
for ii = 1:nrep
    cpu etic 
    T = diag(a,-1) + diag(b) + diag(c,1);
    x0 = T\d;
    cpu etoc backslash1
end
for ii = 1:nrep
	cpu etic
	x1 = ir_apply_tridiag_inv(a, b, c, d);
	cpu etoc apply1    
end
for ii = 1:nrep
	cpu etic
    x3 = tridiag_inv_mex(a, b, c, d);
   	cpu etoc pthr_mex1
end

return;
%% bad inputs

% mixed row and col vectors for a, b, c OK
x3 = tridiag_inv_mex(a, b', c, d);

% bad sizes
x3 = tridiag_inv_mex(a(3:end), b, c, d);
x3 = tridiag_inv_mex(a, b, [c; 1; -1], d);
x3 = tridiag_inv_mex(a, b, c, scale*rand(N+1, M));

% bad types
x3 = tridiag_inv_mex(double(a), b, c, d);
x3 = tridiag_inv_mex(a, b, c, double(d));
x3 = tridiag_inv_mex(a, int16(b), c, d);