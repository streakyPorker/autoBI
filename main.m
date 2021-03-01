[metabolites,transport,reactions]=resolveTXT();

transport_item={};
for i=1:length(transport)
    transport_item{end+1}=transport{i,1};
end

[t,y] =ode45(sumary_eqns,[0,1e3],[0;0;0]);
function rst = conc_eqns()
    % generate 
end

function rst = flux_eqns()
    % generate 
end

function rst = sumary_eqns()
    % generate 
end

