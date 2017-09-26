function mtx_v_norm = norm_mtx_vec(mtx_v)

dd = diag(sqrt(mtx_v*mtx_v'));
mtx_v_norm = mtx_v./repmat(dd,1,size(mtx_v,2));
