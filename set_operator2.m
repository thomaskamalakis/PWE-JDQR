function S0=set_operator2(S);

if strcmp(S.opt.bandtype,'all')
    % Define preconditioner operator
    S.funP=@(Vs) invMVblock(Vs,S);           

    % Define eigenproblem operator
    S.funM=@(Vs) MVblock(Vs,S);       

% or just TE?
elseif strcmp(S.opt.bandtype,'TE')

    % Define preconditioner operator for TE
    S.funP=@(Vs) invMVblockTE(Vs,S);           

    % Define eigenproblem operator for TE
    S.funM=@(Vs) MVblockTE(Vs,S);

% or just TM
elseif strcmp(S.opt.bandtype,'TM')

    % Define preconditioner operator for TM
    S.funP=@(Vs) invMVblockTM(Vs,S);           

    % Define eigenproblem operator for TM
    S.funM=@(Vs) MVblockTM(Vs,S);       
else
    fprintf(1,'Not valid string for bandtype field in structure');
    quit
end
    
S0=S;