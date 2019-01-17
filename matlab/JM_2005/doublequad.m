function Q = doublequad(intfcn,xmin,xmax,ymin,ymax,tol,quadf,varargin)

intfcn = fcnchk(intfcn);

trace = [];
Q = feval(quadf, @innerintegral, ymin, ymax, tol, trace, intfcn, ...
           xmin, xmax, tol, quadf, varargin{:}); 

%---------------------------------------------------------------------------

function Q = innerintegral(y, intfcn, xmin, xmax, tol, quadf, varargin) 

% Evaluate the inner integral at each value of the outer variable. 

Q = zeros(size(y)); 
trace = [];
for i = 1:length(y) 
    Q(i) = feval(quadf, intfcn, xmin, xmax, tol, trace, y(i), varargin{:});
end 

Q = max(Q,0);