function rst = dRdt(J_vec,R_vec,alpha_vec,conc)
    rst=1/conc*(sum(J_vec.*R_vec.*alpha_vec.*[1,-1,-1,1])+R_vec(3)*sum(J_vec.*[1,-1,-1,1]));
end