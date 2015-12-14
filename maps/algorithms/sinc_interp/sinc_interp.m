% From http://phaseportrait.blogspot.com/2008/06/sinc-interpolation-in-matlab.html

% Ideally "resamples" x vector from s to u by sinc interpolation
function y = sinc_interp(x,s,u)
    % Interpolates x sampled sampled at "s" instants
    % Output y is sampled at "u" instants ("u" for "upsampled")
    % (EXPECTS x, s, and u to be ROW VECTORS!!)

    % Find the period of the undersampled signal
    T = s(2)-s(1);

    % When generating this matrix, remember that "s" and "u" are
    % passed as ROW vectors and "y" is expected to also be a ROW
    % vector. If everything were column vectors, we'd do.
    %
    % sincM = repmat( u, 1, length(s) ) - repmat( s', length(u), 1 );
    %
    % So that the matrix would be longer than it is wide.
    % Here, we generate the transpose of that matrix.
    sincM = repmat( u, length(s), 1 ) - repmat( s', 1, length(u) );

    % Equivalent to column vector math:
    % y = sinc( sincM'/T )*x';
    y = x*sinc( sincM/T );
end

n_iline = length(unique(ilxl_read(:,1)));
n_xline = length(unique(ilxl_read(:,2)));
n_samp = size(traces,1);

traces = reshape(n_xline,n_iline,size(traces,1));

%bsxfun(fun,A,B)

y = interpft(x,n,dim)
