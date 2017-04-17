function q_v = unwrap_from_ref(p_v,ref)

N = numel(p_v);

p_v1 = p_v(1:ref);
p_v2 = p_v(ref:N);

q_v1 = flip(unwrap(flip(p_v1)));
q_v2 = unwrap(p_v2);

q_v = zeros(size(p_v));
q_v(1:ref) = q_v1;
q_v(ref:N) = q_v2;
