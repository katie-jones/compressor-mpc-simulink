function results = read_cpp_results(fname_res, fname_info, n_states, n_inputs, n_outputs)

f = fopen(fname_res,'r');

A = fscanf(f, '%f',[n_states+n_inputs+n_outputs+1,inf])';

results.t = A(:,1);
results.x = A(:,2:1+n_states);
results.y = A(:,2+n_states:1+n_states+n_outputs);
results.u = A(:,2+n_states+n_outputs:end);

fclose(f);

f = fopen(fname_info,'r');

A = fscanf(f,'%f');
results.uwt = reshape(A(1:n_inputs*n_inputs),n_inputs,n_inputs);
results.ywt = reshape(A(n_inputs*n_inputs+(1:n_outputs*n_outputs)),n_outputs,n_outputs);
results.yref = A(n_inputs*n_inputs+n_outputs*n_outputs+(1:n_outputs));

end