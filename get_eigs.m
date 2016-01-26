function [eigs,Q] = get_eigs(A,numblocks)
% function [eigs,Q] = get_eigs(A,numblocks)
% gets eigenvalues of BCCB matrix
% assume if it has multiple blocks, they are all identical

if isfield(A.arg, 'blocks')
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
        eigs = repmat(Q*BB(:,1),numblocks,1); %*sqrt(size(B,1)); ,
        
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
        Q = Gdft('mask',true(A.arg.Nx,A.arg.Ny));
        Qcells = repmat({Q},A.arg.Nc,1);
        Qbig = block_fatrix(Qcells);
        eigs = Qbig*AA(:,1);
        eigblock = eigs(1:A.arg.Nx*A.arg.Ny);
        eigs = repmat(eigblock, A.arg.Nc, 1);
        Q = Qbig;
        
        % can I just take a shortcut and say it's Nx*Ny*ones(Nx*Ny*Nx,1)?
end
