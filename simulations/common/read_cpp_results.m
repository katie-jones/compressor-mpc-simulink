function results = read_cpp_results(filename, n_states, n_inputs, n_outputs)

f = fopen(filename,'r');

A = fscanf(f, '%f',[n_states+n_inputs+n_outputs+1,inf])';

results.t = A(:,1);
results.x = A(:,2:1+n_states);
results.y = A(:,2+n_states:1+n_states+n_outputs);
results.u = A(:,2+n_states+n_outputs:end);

fclose(f);

end