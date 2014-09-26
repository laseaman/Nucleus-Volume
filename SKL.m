function Div = SKL(P, Q)
%Symmetrized Kullback-Leibler Divergence
% works for vectors or matricies of the same size adds up all elements.
    KLPQ = P.*log(P./Q);
    KLQP = Q.*log(Q./P);
    Div = 1/2*( nansum(KLQP(:)) + nansum(KLPQ(:)) );
    
    Div1 =  1/2*nansum(log(P./Q).*P)+1/2*nansum(log(Q./P).*Q);
end
