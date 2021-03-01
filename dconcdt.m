function rst = dconcdt(J_vec)
% J_vec denotes : [jin_plus,jin_minus,jout_plus,jout_minus]
    rst=dot(J_vec,[1,-1,-1,1]);
end