% Tests for tri_inverse2

load tri_inverse2_fail;
V = tri_inverse2(a+z, b, cslogb, idx);
assert(~any(isnan(V(:))));
disp('PASSED');
