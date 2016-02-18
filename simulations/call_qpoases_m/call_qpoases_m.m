function [ x,fval,exitflag,iter,lambda ] = call_qpoases_m( H,g,A,lb,ub,lbA,ubA )
%#codegen
%#eml

    assert( isa(H,'double') );
    assert( isa(g,'double') );
    assert( isa(A,'double') );
    assert( isa(lb,'double') );
    assert( isa(ub,'double') );
    assert( isa(lbA,'double') );
    assert( isa(ubA,'double') );

%     coder.varsize( 'H',  [30,30] );
%     coder.varsize( 'g',  [30,1] );
%     coder.varsize( 'A',  [91,30] );
%     coder.varsize( 'lb', [30,1] );
%     coder.varsize( 'ub', [30,1] );
%     coder.varsize( 'lbA',[91,1] );
%     coder.varsize( 'ubA',[91,1] );
    
    
    eml.varsize( 'x',        [100,1] );
    eml.varsize( 'fval',     [1,1] );
    eml.varsize( 'exitflag', [1,1] );
    eml.varsize( 'iter',     [1,1] );
    eml.varsize( 'lambda',   [100,1] );
    
    [nC,nV] = size(A);
    AA = A';
    x = zeros(nV,1);
    fval = 0;
    exitflag = 0;
    iter = 0;
    lambda = zeros(nC+nV,1);
    

    eml.ceval(    'call_qpoases', ...
                    nV, nC, ...
                    eml.rref(H), ...
                    eml.rref(g), ...
                    eml.rref(AA), ...
                    eml.rref(lb), ...
                    eml.rref(ub), ...
                    eml.rref(lbA), ...
                    eml.rref(ubA), ...
                    eml.ref(x), ...
                    eml.ref(fval), ...
                    eml.ref(exitflag), ...
                    eml.ref(iter), ...
                    eml.ref(lambda) ...
                    );

end
