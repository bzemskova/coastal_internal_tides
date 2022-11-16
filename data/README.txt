filename: 'hcrit_forcing_sweep.mat'
	objective: find topographic criticality parameter (hcrit) and Baines body forcing (max and min values)
		over a range of pycnocline widths (len=np) and depths (len=nd)
		both for \omega=1.25f=1.16e-4 1/s and \omega=1.41e-4 1/s (M2 tide)
		same topography/domain set-up as in the manuscript
	variables:
		f: Coriolis parameter
		xx[Nx,Nz]: x-coordinates
		z[Nz,Nz]: z-coordinates
		mupyc_vec[np,1]: vector of pycnocline widths
		zpyc_vec[nd,1]: vector of pycnocline depths
		Fmax_125f[nd,np]: max Baines forcing for \omega=1.25f=1.16e-4 1/s
		Fmin_125f[nd,np]: min Baines forcing for \omega=1.25f=1.16e-4 1/s
		Fmax_M2[nd,np]: max Baines forcing for \omega=1.41e-4 1/s (M2 tide)
		Fmin_M2[nd,np]: min Baines forcing for \omega=1.41e-4 1/s (M2 tide)
		hcrit_125f[nd,np]: topographic criticality parameter for \omega=1.25f=1.16e-4 1/s
		hcrit_M2[nd,np]: topographic criticality parameter for \omega=1.41e-4 1/s (M2 tide)
		hcrit_M2_600_200[1,Nx]: topographic criticality as a function of cross-shore distance 
			for \omega=1.41e-4 1/s (M2 tide)
			and pycnocline width=200m, depth=600m

filename: 'up_sweep_M2.mat'
	objective: find cross-shore energy flux using Baines body force for
			 \omega=1.41e-4 1/s (M2 tide)
			over a range of pycnocline widths (len=np) and depths (len=nd);
			same topography/domain set-up as in the manuscript
	variables:
		XFlux_offshore_sweep[nd,np]: value of vertically-integrated cross-shore energy flux (W/m) 
			at the off-shore boundary
		XFlux_sweep[Nx/2,nd,np]: value of vertically-integrated cross-shore energy flux (W/m)
			in cross-shore direction (at each x)
		up_mode_sweep[Nz-1,nd,np]: value of vertically-integrated cross-shore energy flux (W/m)
			at the off-shore boundary for each vertical mode (1:Nz-1)
		mode1_sweep[nd,np]: value of vertically-integrated cross-shore energy flux (W/m)
			at the off-shore boundary for vertical mode 1
			(same for mode2_sweep, mode3_sweep, mode4_sweep, mode5_sweep, mode6_sweep)
		max_mode_sweep[nd,np]: mode number of the vertical mode with the largest value
			of cross-shore energy flux (W/m) at the off-shore boundary
		max_mode_val_sweep[nd,np]: value of vertically-integrated cross-shore energy flux (W/m)
			at the off-shore boundary for the vertical mode with the largest flux value
		Nx: grid points in x
		Nz: grid points in z
		xx[Nx,Nz]: x-coordinates
		z[Nz,Nz]: z-coordinates


filename: 'up_sweep_125f.mat'
	objective: find cross-shore energy flux using Baines body force for
			 \omega=1.25f=1.16e-4 1/s
			over a range of pycnocline widths (len=np) and depths (len=nd);
			same topography/domain set-up as in the manuscript
	variables:
		same as 'up_sweep_M2.mat'


filename: 'example_data_600_200.mat'
	objective: example output from a simulation with pycnocline depth=600m, width=200m,
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
