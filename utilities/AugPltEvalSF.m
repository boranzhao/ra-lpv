%% evaluate a generalized plant by plugging in theta value 
function Gaev = AugPltEvalSF(Ga,Par)
%% Ga is augmented plant with symbol variable;
%% Ganv is augmented plant evaluated
if length(Par) == 1
    Gaev.A = double(subs(Ga.A,Par));   
    Gaev.B1 = double(subs(Ga.B1,Par));   
    Gaev.B2 = double(subs(Ga.B2,Par));   
    Gaev.C1 = double(subs(Ga.C1,Par));   
    Gaev.D1 = double(subs(Ga.D1,Par));   
    Gaev.D2 = double(subs(Ga.D2,Par));   
elseif length(Par) == 2
    pbar_s = Par(1); veloc_s= Par(2);
    Gaev.A = double(subs(Ga.A));   
    Gaev.B1 = double(subs(Ga.B1));   
    Gaev.B2 = double(subs(Ga.B2));   
    Gaev.C1 = double(subs(Ga.C1));   
    Gaev.D1 = double(subs(Ga.D1));   
    Gaev.D2 = double(subs(Ga.D2));       
else
    error('Cannot deal with more than two schduling parameters');  
end

