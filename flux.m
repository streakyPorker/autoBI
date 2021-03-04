function J = flux(k,E,concr_vec,Kr_vec,nr_vec,concp_vec,Kp_vec,np_vec)
    coef=k*E;
    numerator=prod((concr_vec./Kr_vec).^nr_vec);
    denominator=1+numerator+prod((concp_vec./Kp_vec).^np_vec);
    J=coef*numerator/denominator;
end