function B = PDEMatrixB(Q)
% Non-conservative product matrix
global pdetype gamma g nVar
switch(pdetype)
    case(0)
        B = zeros(nVar,nVar);
    case(1)
        B = zeros(nVar,nVar);
    case(2)
        B = zeros(nVar,nVar);
        B(2,3) = g*Q(1);
end

end