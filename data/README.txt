filename: 'alpha_Ftilde.mat'
	objective: nondimensional topographic criticality parameter \alpha
			 and Baines body force \F_{tilde} 
			over a range of pycnocline widths (len=np) and depths (len=nd);
			same topography/domain set-up as in the manuscript
	variables:
		mupyc_vec[1,np]: vector of pycnocline widths (m)
		Zpyc_vec[1,nd]: vector of pycnocline depths (m)
		alpha_30deg[nd,np]: topographic criticality parameter \alpha computed at f=0.73e-4
					(30 degree latitude)
		alpha_40deg[nd,np]: topographic criticality parameter \alpha computed at f=0.93e-4
					(40 degree latitude)
		alpha_50deg[nd,np]: topographic criticality parameter \alpha computed at f=1.11e-4
					(50 degree latitude)
		F_Baines[nd,np]: Baines body force (s^-2)


filename: 'coriolis_sweep.mat'
	objective: find depth-integrated cross-shore energy flux at the off-shore boundary
			as a function of Coriolis parameter (f) (len=nf)
			for non-constant stratification (pycnocline depths: 300m, 400m, 600m, 800m)
			and for constant stratification (N^2 = 1.22e-5 s^-2);
			same topography/domain set-up as in the manuscript
	variables:
		**all non-constant stratification cases are structured the same**
		**below is example for pycnocline depth of 300m**
		**for each pycnocline depth, there are two different pycnocline widths (np=2)**
		f_vec[nf,1] = Coriolis parameter f (s^-1)
		alpha_300_fsweep[np,nf]: topographic criticality parameter \alpha
		mupyc_300[np,1]: pycnocline widths (m)
		Xflux_300_fsweep[np,nf]: vertically-integrated cross-shore energy flux 
						at the off-shore boundary (W/m)

		alpha_constN2_fsweep[nf,1]: topographic criticality parameter \alpha
						for constant stratification
		Xflux_constN2_fsweep[nf,1]: vertically-integrated cross-shore energy flux 
						at the off-shore boundary (W/m)
						for constant stratification


filename: 'shelfdepth_sweep.mat'
	objective: find depth-integrated cross-shore energy flux at the off-shore boundary
			as a function of shelf depth (h_s) (len=nh)
			and coastal depth (h_c) = [0.5, 1]*h_s
			for non-constant stratification (pycnocline depths: 300m, 400m, 600m, 800m)
			and for constant stratification (N^2 = 1.22e-5 s^-2);
			same topography/domain set-up as in the manuscript
	variables:
		**all non-constant stratification cases are structured the same**
		**below is example for pycnocline depth of 300m**
		**for each pycnocline depth, there are two different pycnocline widths (np=2)**
		**computed for two different coastal depths (nhc=2)**
		hs_vec[1,nh]: shelf depth (m)
		hc_vec[1,nhc]: coastal depth/shelf depth ([0.5, 1])
		alpha_300_hssweep[np,nhc,nh]: topographic criticality parameter \alpha
		mupyc_300[np,1]: pycnocline widths (m)
		Xflux_300_hssweep[np,nhc,nh]: vertically-integrated cross-shore energy flux 
						at the off-shore boundary (W/m)

		alpha_constN2_hssweep[nh,nhc]: topographic criticality parameter \alpha
						for constant stratification
		Xflux_constN2_hssweep[nh,nhc]: vertically-integrated cross-shore energy flux 
						at the off-shore boundary (W/m)
						for constant stratification


filename: 'shelfwidth_sweep.mat'
	objective: find depth-integrated cross-shore energy flux at the off-shore boundary
			as a function of shelf width (x_s) (len=nxs)
			and coastal depth (h_c) = [0.5, 1]*150
			for non-constant stratification (pycnocline depths: 300m, 400m, 600m, 800m)
			and for constant stratification (N^2 = 1.22e-5 s^-2);
			same topography/domain set-up as in the manuscript
	variables:
		**all non-constant stratification cases are structured the same**
		**below is example for pycnocline depth of 300m**
		**for each pycnocline depth, there are two different pycnocline widths (np=2)**
		**computed for two different coastal depths (nhc=2)**
		xs_vec[1,nxs]: shelf width (m)
		hc_vec[1,nhc]: coastal depth/shelf depth ([0.5, 1]); shelf depth(h_s) = 150m
		alpha_300_xssweep[np,nhc,nxs]: topographic criticality parameter \alpha
		mupyc_300[np,1]: pycnocline widths (m)
		Xflux_300_xssweep[np,nhc,nxs]: vertically-integrated cross-shore energy flux 
						at the off-shore boundary (W/m)

		alpha_constN2_xssweep[nxs,nhc]: topographic criticality parameter \alpha
						for constant stratification
		Xflux_constN2_xssweep[nxs,nhc]: vertically-integrated cross-shore energy flux 
						at the off-shore boundary (W/m)
						for constant stratification

filename: 'slopewidth_sweep.mat'
	objective: find depth-integrated cross-shore energy flux at the off-shore boundary
			as a function of slope width (x_W) (len=nxw)
			and coastal depth (h_c) = [0.5, 1]*150
			for non-constant stratification (pycnocline depths: 300m, 400m, 600m, 800m)
			and for constant stratification (N^2 = 1.22e-5 s^-2);
			same topography/domain set-up as in the manuscript
	variables:
		**all non-constant stratification cases are structured the same**
		**below is example for pycnocline depth of 300m**
		**for each pycnocline depth, there are two different pycnocline widths (np=2)**
		**computed for two different coastal depths (nhc=2)**
		xW_vec[1,nxw]: slope width (m)
		hc_vec[1,nhc]: coastal depth/shelf depth ([0.5, 1]); shelf depth(h_s) = 150m
		alpha_300_xWsweep[np,nhc,nxw]: topographic criticality parameter \alpha
		mupyc_300[np,1]: pycnocline widths (m)
		Xflux_300_xWsweep[np,nhc,nxw]: vertically-integrated cross-shore energy flux 
						at the off-shore boundary (W/m)

		alpha_constN2_xWsweep[nxw,nhc]: topographic criticality parameter \alpha
						for constant stratification
		Xflux_constN2_xWsweep[nxw,nhc]: vertically-integrated cross-shore energy flux 
						at the off-shore boundary (W/m)
						for constant stratification


filename: 'Zpyc400_mupyc200_param_sweep.mat'
	objective: find depth-integrated cross-shore energy flux at the off-shore boundary
			for pycnocline depth of 400m and width of 200m
			as a function of Coriolis parameter f and
				1) shelf depth h_s (fix h_c/h_s, x_s, x_W)
			        2) shelf width x_s (fix h_c/h_s, h_s, x_W)
				3) slope width x_W (fix h_c/h_s, h_s, x_s)
			same topography/domain set-up as in the manuscript
	variables:
		f_vec[1,nf]: Coriolis parameter vector (s^-1)
		hs_vec[1,nhs]: shelf depth vector (m)
		xs_vec[1,nxs]: shelf width vector (m)
		xW_vec[1,nxw]: slope width vector (m)
		alpha_hs_sweep[nf,nhs]: topographic criticality parameter \alpha as a function of
						Coriolis parameter and shelf depth
		alpha_xs_sweep[nf,nxs]: topographic criticality parameter \alpha as a function of
						Coriolis parameter and shelf width
		alpha_xW_sweep[nf,nxw]: topographic criticality parameter \alpha as a function of
						Coriolis parameter and slope width
		Xflux_hs_sweep[nf,nhs]: vertically-integrated cross-shore energy flux 
						at the off-shore boundary (W/m) as a function of
						Coriolis parameter and shelf depth
		Xflux_xs_sweep[nf,nxs]: vertically-integrated cross-shore energy flux 
						at the off-shore boundary (W/m) as a function of
						Coriolis parameter and shelf width
		Xflux_xW_sweep[nf,nxw]: vertically-integrated cross-shore energy flux 
						at the off-shore boundary (W/m) as a function of
						Coriolis parameter and slope width


filename: 'example_data_400_200.mat'
	objective: example output from a simulation with pycnocline depth=400m, width=200m,
			 \omega=1.41e-4 1/s (M2 tide), and topography defined in manuscript
	variables:
		AFlux[1,Nx/2]: vertically-integrated along-shore energy flux (W/m)
		XFlux[1,Nx/2]: vertically-integrated cross-shore energy flux (W/m)
		fv[Nx/2,Nz]: along-shore energy flux (W/m^2)
		fu[Nx/2,Nz]: cross-shore energy flux (W/m^2)
		xx[Nx,Nz]: x-coordinates
		z[Nx,Nz]: z-coordinates
		uModeF[Nx/2,Nz]: cross-shore velocity (m/s)
		vModeF[Nx/2,Nz]: along-shore velocity (m/s)
		pModeF[Nx/2,Nz]: pressure (m^2/s^2)
		N2[Nx,Nz]: buoyancy frequency N^2 (1/s^2)
		h[Nx,1]: depth at each cross-shore location (m)
		hx[Nx,1]: horizontal gradient of depth, dh/dx
		F[Nx,Nz]: forcing in pressure (odd x-points) and u (even x-points)
		up_modal_offshore[Nz,Nz]: cross-shore energy flux (W/m^2) 
				in each vertical mode 
				at the off-shore boundary as a function of depth
		up_modal_all[Nz,Nz]: vertically-integrated cross-shore energy flux
		             (W/m) in each mode and cross-mode
		up_modal_integrated[Nz,1]: vertically-integrated cross-shore energy
				flux (W/m) in each vertical mode at the off-shore boundary

		Nx: grid points in x
		Nz: grid points in z
		R: rigid lid approx
		Zpyc : pycnocline depth
		mupyc : pycnocline width
		force_type : Baines body force
		l : along-shore wavenumber; small value close to 0
		L : cross-shore domain extent (m)
		h0 : max depth (m)
		W : slope width (m)
		xs : shelf width (m)
		hc : depth at the coast (m)
		hs : depth at shelf break
		f : Coriolis parameter (1/s)
		sigma : forcing frequency (1/s), M2 tide
		rho0 : background density (kg/m^3)
		g : gravity (m/s^2)


filename: 'example_data_constantN2.mat'
	objective: example output from a simulation constant stratification,
			 \omega=1.41e-4 1/s (M2 tide), and topography defined in manuscript
	variables:
		same format as in 'example_data_400_200.mat'


filename: 'resonance_scan.mat'
	objective: example resonance scan 
			over forcing frequencies (len=nf) and along-shore wavenumbers (len=nl)
			 for a simulation with pycnocline depth=600m, width=200m,
			 and topography defined in manuscript
	variables:
		l_vec2[nl,1]: range of along-shore wavenumbers considered
		SI[nf,1]: range of forcing frequencies considered
		P0p_sd[nl,nf]: pressure response to forcing
		P0u_sd[nl,nf]: cross-shore velocity response to forcing

		Nx: grid points in x
		Nz: grid points in z
		R: rigid lid approx
		Zpyc : pycnocline depth
		mupyc : pycnocline width
		force_type : Baines body force
		l : along-shore wavenumber; small value close to 0
		L : cross-shore domain extent (m)
		h0 : max depth (m)
		W : slope width (m)
		xs : shelf width (m)
		hc : depth at the coast (m)
		hs : depth at shelf break
		f : Coriolis parameter (1/s)
		sigma : forcing frequency (1/s), M2 tide
		rho0 : background density (kg/m^3)
		g : gravity (m/s^2)
