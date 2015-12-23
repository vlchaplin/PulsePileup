
;  Must call link_ppu_library before this routine compiles!!!


;
;		The procedure ppu_fit_example.pro  'fits' some example data.
;		The data were generated from a band function. The example as is tries a set of rates for the input detectors, for fixed shape params.
;		The way to try different parameters is described below.
;




function BandsGRBSpectrum, alpha, beta, ep, energy
	
	
	if n_elements(energy) eq 1 then energy=[energy]
	
	iLt = where( energy lt (ep*(alpha - beta)/(2.0 + alpha)), nLt, COMP=iGt, NCOMP=nGt )
	
	nchan = n_elements(energy)
	
	result=dblarr(nchan)
	
	if ( nLt gt 0 ) then begin
		result[iLt] = exp( - energy[iLt]*(2.0 + alpha)/ep ) * (energy[iLt])^alpha
	end
	
	if ( nGt gt 0  ) then begin
		result[iGt] = exp(beta-alpha) * (energy[iGt]^beta) * (ep*(alpha-beta)/(2.0+alpha))^(alpha-beta) 
	end
	
	return, result

end

function BandsGRB_NumberSpectrum, alpha,beta,ep, channel_edges, channel_dE

	nchan=n_elements(channel_edges)

	dNdE = BandsGRBSpectrum(alpha, beta, ep, channel_edges)

	;trapezoidal rule
	number_spectrum = channel_dE/2.0*(dNdE[0:nchan-2] + dNdE[1:nchan-1]) 
	
	return, number_spectrum
end

function normal_dist, X, g_mean, g_sigma

	return,exp( - (X - g_mean)^2 / (2*g_sigma)^2 )/(g_sigma*Sqrt(2*!pi))
end

function Cstat, model, data


	LL = 0.0d


	; IDL forces simple things to be convoluted
	
	stirling_idx = where(data ge 90.0, use_stirling, COMP=pois_idx,NCOMP=use_pois)

	if use_stirling gt 0 then begin
		
		
		mu = model[stirling_idx]
		d = data[stirling_idx]
	
		iz = where(d eq 0,nz, COMP=idx_nonz, NCOMP=nnonz)
		
		if nz gt 0 then LL += Total( -mu[iz] )
		
		if nnonz gt 0 then LL += Total( -mu[idx_nonz] + d[idx_nonz]*Alog(mu[idx_nonz]) + d[idx_nonz]*(1.0 - Alog(d[idx_nonz])) )
		
		;stop
	end
	
	if use_pois gt 0 then begin
	
		mu = model[pois_idx]
		d = data[pois_idx]
	
		iz = where(d eq 0, nz, COMP=idx_nonz, NCOMP=nnonz)
	
		if nz gt 0 then LL += Total( -mu[iz] )
		
		if nnonz gt 0 then begin
		
			P = exp(-mu[idx_nonz]) * (mu[idx_nonz])^(d[idx_nonz]) /  factorial(d[idx_nonz])
		
		
			nonzero = where(P gt 0,c)
			if c gt 0 then LL += Total( Alog(P[nonzero]) )
		
		end
	
	end
	
	return, -2.0 * LL

end



pro ppu_fit_example, stat_grid, H5=h5file

	;h5file = '/gbm/TGF/PPU/C_PPU/PPUprob_Volts.h5'
		
	if n_elements(h5file) eq 0 then $
		h5file = dialog_pickfile(FILTER=['*.h5'], TITLE='Select the HDF5 pileup response file')

	print, "The pileup response file is: "+h5file
	
	status = PPUSetMaxOrder(8)
	status = PPULoadH5file(h5file)
	if status ne 0 then begin
		print, "ERROR: unable to load HDF5 file: "+h5file
		return
	end

	;determine the upper and lower bounds intrinsic to the file
	status = PPUGetH5Channels(h5file, h5edges)
	h5nChan=n_elements(h5edges)-1

	h5lowerBound = h5edges[0]        ;  should be <= the LLT. Our file has 0.005 volts
	h5upperBound = h5edges[h5nChan]  ;  should be larger than 5V.  our file has 6V. 

	;data_phas = ['/gbm/TGF/PPU/C_PPU/TestingData/b1_data.pha']

	if n_elements(data_phas) eq 0 then $
		data_phas = dialog_pickfile(FILTER=['*.pha'], TITLE='Select data files: single-spectrum (type-1 PHA) fits file', /Multiple)
	
	data_phas = [data_phas]
	
	if data_phas[0] eq '' then return
	
	print, "The files contain detector data that we will try to 'fit'"
	print, data_phas
	
	
	; First, loop over the input PHA data and DRMs to load them
	nfiles = n_elements(data_phas) 

	drms = strarr(nfiles)
	numSetsToCalc=0	
	
	spec_data=PTRARR(nfiles)         ; spectra 
	chan_volts=PTRARR(nfiles)		; channels for spectrum i converted to volts
	nchans =LONARR(nfiles)        ; number channels for spectrum i
	integration_times_us = dblarr(nfiles)  ;exposures in microseconds
	
	drmdata = OBJARR(nfiles)
	
	
	
	;exposure in the optical sense -- binsize, not livetime
	exposure = 1.0
	
	for i=0, nfiles-1 do begin
	
		
		pha = data_phas[i]
		
		; Add a reference to this detector or file into the PPU instance hash. If one already exists this does no harm...
		status = PPUCreateInstance( data_phas[i], h5file )
		
		if status ne 0 then begin
			print, "error making instance"
			continue
		end
		
		
		;GUI to select DRM
		drm = dialog_pickfile(FILTER=['*.rsp*'], TITLE='Choose DRM for '+file_basename(pha))
	
		drms[i] = drm
		
		if drm eq '' then continue
		
		
		
		spectable = mrdfits( pha, 2, header, /use_colnum, /fscale)
		
		
		spectype = strtrim( sxpar(header, 'HDUCLAS3'), 2)
		filetype = sxpar(header, 'HDUCLAS4',COUNT=ok)
		
		if ok && strtrim(filetype,2) eq 'TYPEII' then have_type2=1 else have_type2=0 
		
		if have_type2 then begin 
		
			; PHA II read (assumes second bin contains source+bg, first bin only background...
			
			; use column numbers so we don't have to figure out whether it's called COUNTS or RATES
			spectrum_column = (spectable.c1)
						
			background = spectrum_column[*,0] ; not used in this example
			spectrum = spectrum_column[*,1]
			
			; do this to resolve the pha2 data structure... sometimes the stat_err column is present
			; stupid mrdfits and stupid idl
			
			spectable = mrdfits( pha, 2, header, /fscale)
			
			integrations = spectable.endtime - spectable.time 
			
			;store the time for the src+bg bin in microseconds
			integration_times_us[i] = 1e06*integrations[1]
			
			if spectype eq 'RATE' then begin
			
				;convert to counts 
				spectrum *= integrations[1] 
			
			end
			
			
		end else begin
			
			; PHA I read
			spectrum = spectable.c2
		
		
			if spectype eq 'RATE' then begin
		
				;convert to counts
			
				tstart = strtrim( sxpar(header, 'TSTART'), 2)
				tstop = strtrim( sxpar(header, 'TSTOP'), 2)
				
				integration = tstop-tstart
				
				print, 'Converted input spectrum to counts using TSTOP-TSTART'
				spectrum *= integration
				
				integration_times_us[i] = 1e06*integration
			end
		
		end
		
		
		spec_data[i] = ptr_new(spectrum)
		
		
		;object to read .rsp file
		drmOBJECT = obj_new('drmdata', drm)
		
		drmdata[i] = drmOBJECT
		
		numSetsToCalc++
	
		; Get DRM channel edges
		; convert channel edges to volts
		chans_kev = double(drmdata[i]->ch(/edges))
		nchan = n_elements(chans_kev)-1
		
		LLTvolts = 5.0*24.0/4096.0
		LLTkeV = chans_kev[0]
		ULTkeV = chans_kev[nchan-1] ; this is so the upper-level threshold hits the start of the last channel
		
		;linear map into voltage range
		edges_volts = double(LLTvolts + (5.0 - LLTvolts)*(chans_kev - LLTkev)/(ULTkev - LLTkev))
		
		
		chan_volts[i] = ptr_new(edges_volts)
		nchans[i] = nchan
		
		; Initialize with a dummy spectrum just to get arrays allocated...
		
		counts=dblarr(nchan)
		status = PPUInitZerothOrderCounts(data_phas[i], edges_volts, counts, nchan)
	
		if status ne 0 then begin
			print, "ERROR initializing"
			return
		end
		
		
	
	endfor
		
	
	;set-up the sample grid
	
	;1 element means these are fixed
	alphas = [-0.5d]
	betas = [-2.4d]
	
	; map of Epeak parameter
	minEp=800.0d
	maxEp = 1200.0d
	step = 50.0d
	epeaks = dindgen(floor((maxEp-minEp)/step)+1)*step + minEp
	
	epeaks = [1000.0d] ;just fix this one for initial example
	
	; map of rate parameter
	; If using multiple detectors, each will have a different rate (similar to effective area correction)
	; Thus it might be good to have one grid per detector...
	minR=0.10
	maxR =0.35
	step =0.005	
	rates = dindgen(floor((maxR-minR)/step)+1)*step + minR
	
	;rates=[0.215d] 
	
	Nal = n_elements(alphas)
	Nbe = n_elements(betas)
	Nep = n_elements(epeaks)
	Nrt = n_elements(rates)
	
	;array to hold c-stat values from the parameter grid
	; If using multiple detectors, each will have a different rate (similar to effective area correction)
	; Thus it might be good to have one grid per detector... such as:
	
	stat_grid = dblarr(nfiles, Nal, Nbe, Nep, Nrt)
	;it can be recombined later
	
	
	;loop over epeaks
	for e=0, Nep-1 do begin
	
		ep = epeaks[e]
	
		;loop over beta
		for b=0, Nbe-1 do begin
		
			beta = betas[b]
		
			; loop over alpha
			for a=0,Nal-1 do begin
			
				alpha = alphas[a]
			
				; loop over detectors
				for i=0, nfiles-1 do begin
				
					; photedge_ptrs[0] - pointer to photon edges
					; photedge_ptrs[1] - pointer to photon bin-sizes (dE)
					
					photedge_ptrs = drmdata[i]->ph()
				
					;this band function isn't normalized
					phot_model = BandsGRB_NumberSpectrum( alpha, beta, ep, *photedge_ptrs[0], *photedge_ptrs[1] )
					
					zeroth_order_count_model = drmdata[i]->mat() # phot_model
					
					zeroth_order_count_model /= total(zeroth_order_count_model)
					
					edges_volts = *(chan_volts[i])
					
					status = PPUInitZerothOrderCounts( data_phas[i], edges_volts, zeroth_order_count_model, nchans[i] )
					
					
					; now loop over rates and compute the model, then compare to data...
					
					for r=0, Nrt-1 do begin
						
						rate_per_us = rates[r]
						
						mean_counts = integration_times_us[i] * rate_per_us 
						
						status = PPUComputeSpectrum(data_phas[i], rate_per_us, edges_volts, nchan, ppu_count_model)
						
						
						count_data =  round(*(spec_data[i]))
						ppu_count_model *= mean_counts
						
						
						
						;evaluate fit statistic for channels 3 through 124
						C = Cstat( ppu_count_model[3:124], count_data[3:124]  )
						
						stat_grid(i,a,b,e,r) += C
						
						;stop
					
					endfor
			
				endfor	
				
			endfor
		endfor
	endfor
	

	;free all memory associated with the pileup model and h5file. To free only certain objects, use PPUDeleteInstance() instead.
	s=PPUDeleteAll()
	
	; hopefully the parameter grid has a well-defined minimum
	; stat_grid(i,a,b,e,r)
	
	;in this case we just did a 1-d search for the correct rate, so it can be plotted simply:
	
	plot, rates, stat_grid[*]
	
	
	
end