
; MUST CALL link_ppu_library BEFORE THIE ROUTINE COMPILES FIRST TIME!!

;  Read the comments to get an idea what's happening

pro ppu_example, edges_kev, inputSpec, outputSpec

	;use the volts units file... path to the h5 (aka the HDF5) file to use
	h5file = './PPUProb_Volts.h5'



	;set the maximum order BEFOR CALLING PPULoadH5File().  For as long as h5file is loaded,
	; no approximation orders greater than this can be used.
	status = PPUSetMaxOrder(5)


	
	;load the h5 file
	status = PPULoadH5File(h5file)
	
	if status ne 0 then begin
		print, "ERROR loading file"
		return
	end


	
	;;  To make the file persist until it is explicitly deleted, do this:
	status = PPUCreateInstance(h5file, h5file)  
	
	

	; load the intrinsic binning of the h5 pileup response.  since we're using 
	;    the _Volts.h5 file, these edges are in volts
	;
	;  this is necessary so our results will be properly normalized-- the pileup "responses"
	;  are computed over a certain voltage (or energy) range, and are normalized over this range.
	;
	
	status = PPUGetH5Channels(h5file, h5edges)
	h5nChan=n_elements(h5edges)-1

	h5upperBound = h5edges[0]        ;  should be <= the LLT. Our file has 0.005 volts
	h5lowerBound = h5edges[h5nChan]  ;  should be larger than 5V.  our file has 6V. 




	;create model instance named after a detector
	detector="some_detector_name"
	status = PPUCreateInstance(detector, h5file)
	
	;set the calculation order for this detector.. different detectors can have different calculation orders
	status = PPUSetApproxOrder(detector, 2)
	
	;set up the voltage llt
	LLTVolts = 24.0d / 4096.0d
	
	;  since we have the PPUprob_Volts file loaded, the spectral channels must be in volts.
	;
	;  must be DOUBLE type!!!
	
	; make a logarithmic binning in voltage such that the last channel BEGINS at 5.0 V (ie., it's the overflow channel).
	; in practice this would be replaced by the detector's EBOUNDS table converted to volts
	
	nchan = 128
	edges_volts = LLTVolts* (5.0 / LLTVolts)^(dindgen(nchan+1)/(nchan-1))
	
	;;;  This defines a crude voltage-to-energy conversion.  It is a linear mapping from LLT to 5V onto a specified keV range

	;input the (known from calibration) equivalent keV of the LLT
	LLTinKeV = 4.0d
	
	;known kev of the overflow channel (=5.0 V)
	OVERkeV = 1000.0d
	
	;Linear map from voltage to energy. Result will be DOUBLE type
	edges_kev = LLTinKeV + (edges_volts - LLTVolts)*(OVERkeV - LLTinKeV)/(5.0d - LLTVolts)
	
	dE = edges_kev[1 : nchan] - edges_kev[0 : nchan-1]

	;... code to convert photons to counts from drm...
	;
	;   For this example, just define the counts spectrum (just a power law for now). Trapezoidal integration.
	;   The normalization is not important here
	counts = dE*((edges_kev[0 : nchan-1] / 100.0d)^(-2.0d)  +   (edges_kev[1 : nchan] / 100.0d)^(-2.0d)) / 2.0d


	
	
	; Initialize the model.. be sure to pass the volts channels 
	status = PPUInitZerothOrderCounts(detector, edges_volts, counts, nchan)

	if status ne 0 then begin
		print, "ERROR initializing"
		return
	end

	rate = 0.1d ; counts per microsecond (0.1 = 100,000 cps)




	;return spectrum as a probability function (probability per channel…must be renormalized).
	;
	; outputSpec will contain a probability function (a probability spectrum... not a PDF!!).
	;    It normalizes <1.0 because of PPU+DT losses.
	
	status = PPUComputeSpectrum(detector, rate, edges_volts, nchan, outputSpec)




	;re-normalize the output to be in terms of the rate
	
	outputSpec *= rate   ; this is counts per microsecond... multiply by 1e06 to convert to cps
	
	
	
	
	; bin-centers
	centers = (edges_kev[0:nchan-1] + edges_kev[1:nchan]) / 2.0




	;re-normalize the input spectrum unless it was already done
	inputSpec = rate * counts/Total(counts)





	;show results

	plot_oo, centers, inputSpec
	oplot, centers, outputSpec


	print, '% change [output - input / input] = ', (outputSpec - inputSpec) / inputSpec
stop

	; Free memory. Once instances have been created they must all be deleted which will then free the 
	;   h5file object
	
	status = PPUDeleteInstance(detector)

	; If there was only one associate to h5file, deleting it would caus the file to be unloaded too
	; However, we made bound associations in the above code:
	;   PPUCreateInstance(detector, h5file) and PPUCreateInstance(h5file, h5file)
	; so the second one has to be deleted as well:
	
	status = PPUDeleteInstance(h5file)
	; Now all memory has been freed.
	;
	; If we didn't delete these objects would persist outside of this procedure, and could be accessed from 
	;   another context using the same PPU.. routines passing the same filename or detector string
	
end