find_t_cutoff <- function(phi_U, t_phi_U){
	
	# Find smallest t > 0 at which phi_U reaches its largest local maximum
	ind <- which(t_phi_U >= 0)
	d_phi_U <- phi_U[ind[2:length(ind)]] - phi_U[ind[1:(length(ind) - 1)]]

	if (length(which(d_phi_U >= 0)) == 0) {
		# phi_U is always decreasing for positive t
		t_cutoff <- t_phi_U[length(t_phi_U)]
	} else {
		first_min_ind <- ind[min(which(d_phi_U >= 0))]
		phi_U_threshold <- max(phi_U[ind[ind>first_min_ind]])
		tmp <- t_phi_U[phi_U < phi_U_threshold]
		t_cutoff <- min(tmp[tmp > 0])	
	}	

	t_cutoff
}