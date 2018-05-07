find_t_cutoff <- function(phi_U, t_phi_U, n){
	
	# Find smallest t > 0 at which phi_U reaches its largest local maximum
	ind <- which(t_phi_U >= 0)
	d_phi_U <- phi_U[ind[2:length(ind)]] - phi_U[ind[1:(length(ind) - 1)]]
	first_min_ind <- ind[min(which(d_phi_U >= 0))]
	print(first_min_ind)
	phi_U_threshold <- max(phi_U[ind[ind>first_min_ind]])
	print(phi_U_threshold)

	# Find first t > 0 that falls below this threshold
	tmp <- t_phi_U[phi_U < phi_U_threshold]
	if (length(tmp) > 1) {
		t_cutoff <- min(tmp[tmp > 0])	
	} else {
		t_cutoff <- Inf
	}

# Find first t > 0 that falls below this threshold
	tmp <- t_phi_U[phi_U < n^(-1/4)]
	if (length(tmp) > 1) {
		t_cutoff <- min(tmp[tmp > 0])	
	} else {
		t_cutoff <- Inf
	}

	print(t_cutoff)
	t_cutoff
}