
pro link_ppu_library

	; path of compiled library... be sure to include its directory into $DYLD_LIBRARY_PATH
	library = '/gbm/TGF/PPU/C_PPU/libIDLppu.dylib'
	
	
	
	LINKIMAGE, 'PPUSetMaxOrder', library, 1, 'PPUSetMaxOrder', MAX_ARGS=1, MIN_ARGS=1
	
	LINKIMAGE, 'PPULoadH5File', library, 1, 'PPULoadH5File', MAX_ARGS=1, MIN_ARGS=1
	
	LINKIMAGE, 'PPUCreateInstance', library, 1, 'PPUCreateInstance', MAX_ARGS=2, MIN_ARGS=2
	
	LINKIMAGE, 'PPUSetApproxOrder', library, 1, 'PPUSetApproxOrder', MAX_ARGS=2, MIN_ARGS=2
	
	LINKIMAGE, 'PPUGetH5NumChannels', library, 1, 'PPUGetH5NumChannels', MAX_ARGS=2, MIN_ARGS=2
	
	LINKIMAGE, 'PPUGetH5Channels', library, 1, 'PPUGetH5Channels', MAX_ARGS=2, MIN_ARGS=2
	
	LINKIMAGE, 'PPUInitZerothOrderCounts', library, 1, 'PPUInitZerothOrderCounts', MAX_ARGS=4, MIN_ARGS=4
	
	LINKIMAGE, 'PPUComputeSpectrum', library, 1, 'PPUComputeSpectrum', MAX_ARGS=5, MIN_ARGS=5
	
	LINKIMAGE, 'PPUDeleteInstance', library, 1, 'PPUDeleteInstance', MAX_ARGS=1, MIN_ARGS=1
	
	LINKIMAGE, 'PPUDeleteAll', library, 1, 'PPUDeleteAll', MAX_ARGS=0
end