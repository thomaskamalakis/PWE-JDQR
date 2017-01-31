function S0=st_orth_vecs(k,S);

[S.vec_a,S.vec_b,S.kplus_M]=orth_vecs(k,S.Gx,S.Gy,S.Gz); 
S0=S;
