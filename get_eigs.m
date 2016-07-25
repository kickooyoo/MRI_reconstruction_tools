function [eigvs,Q] = get_eigs(A, numblocks, varargin)
% function [eigs,Q] = get_eigs(A, numblocks)
% gets eigenvalues of BCCB matrix
% assume if it has multiple blocks, they are all identical
arg.Nx = 0;
arg.Ny = 0;
arg = vararg_pair(arg, varargin);

if isfield(A.arg, 'blocks') % legacy code
        % assume it was constructed using block_fatrix (or construct_fourierfat.m)
        
        if (numblocks > 1)
                B = A.arg.blocks{1};
                % if blocks are different, should iterate over A.arg.blocks{ii}
        elseif (numblocks == 1)
                B = A;
        else
                error('invalid number of blocks');
        end
        
        if isa(B,'fatrix2')
                nx = B.odim(1);
                ny = B.odim(2);
                %display('stop here');
                %keyboard;
        else
                nx = size(B,1);
                ny = size(B,2);
        end
        
        BB = B'*B;
        Q = Gdft('mask',true(nx,ny));
        if ~isa(B,'fatrix2')
                display('figure out this bug');
                Qtest = Gdft('mask',true(nx,1));
                keyboard;
        end
        eigvs = repmat(Q*BB(:,1),numblocks,1); %*sqrt(size(B,1)); ,
        
        %test code
        if (numblocks > 1)
                
                Qcells = repmat({Q},numblocks,1);
                Qbig = block_fatrix(Qcells);
                %    %test1 = (1/(N1*N2))*Qbig*diag(eigs)*Qbig';
                %    test2 = (1/(nx*ny))*Qbig'*diag(eigs)*Qbig;
                Q = Qbig;
                %else
                %    %test1 = (1/(N1*N2))*Q*diag(eigs)*Q';
                %    test2 = (1/(nx*ny))*Q'*diag(eigs)*Q;
        end

else % assume constructed with GsplineF

        AA = A'*A;
        %Q = Gdft('mask',true(A.arg.Nx,A.arg.Ny));
        %Qcells = repmat({Q},A.arg.Nc,1);
        %Qbig = block_fatrix(Qcells);
        if isfield(A.arg, 'Nx') && isfield(A.arg, 'Ny')
                Q = GsplineF(A.arg.Nx, A.arg.Ny, 1, A.arg.Nc); % no real difference between block_fatrix method and GsplineF
                e0 = zeros(A.arg.Nx, A.arg.Ny, A.arg.Nc);
        else
                Q = GsplineF(arg.Nx, arg.Ny, 1,numblocks);
                e0 = zeros(arg.Nx, arg.Ny, numblocks);
        end
        % Q = F_par(A.arg.Nx, A.arg.Ny, A.arg.Nc); way slower :'(

	
	e0(1,1,:) = 1;
        eigvs = Q * (AA * e0(:));
end
