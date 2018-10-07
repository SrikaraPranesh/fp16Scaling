function varargout = svd(A)
    if nargout <= 1
        varargout{1} = fp16(svdtx(A));
    else
        [U,s,V] = svdtx(A);
        varargout{1} = fp16(U);
        varargout{3} = fp16(V);
        [n,p] = size(A);
        if n >= p
            varargout{2} = fp16([diag(s); zeros(n-p,p)]);
        else
            varargout{2} = fp16([diag(s) zeros(n,p-n)]);
        end                
    end
end
% fp16/svd