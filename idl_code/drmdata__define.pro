FUNCTION binplace, v,i, edges, nedges
	
	IF v GE edges[nedges-1] THEN RETURN,nedges-1
	IF v LT edges[0] THEN RETURN, 0
	IF i+1 GT nedges THEN i=binplace(v,i-1,edges,nedges)
	
	IF v GE edges[i] && v LT edges[i+1] THEN RETURN,i
	
	IF v LT edges[i] THEN i=binplace(v,i-1,edges,nedges) $
	ELSE IF v GT edges[i] THEN i=binplace(v,i+1,edges,nedges)
	
	RETURN,i
END


FUNCTION drmdata::init, file, _REF_EXTRA=ex

	IF N_ELEMENTS(file) NE 1 THEN RETURN, 1
	
	self->loadrsp, file, _EXTRA=ex
	
	RETURN,1
END

FUNCTION drmdata::getFile
	return, self.file
END

FUNCTION drmdata::getDet
	return, self.det
END

PRO drmdata::loadrsp, file, _REF_EXTRA=ex
	ON_ERROR,2
	
	
	IF N_ELEMENTS(file) EQ 1 THEN self.file = file
	
	prm = mrdfits(self->file(), 0, phdr, _EXTRA=ex)
	ebounds = mrdfits(self->file(), 1, ebhdr, _EXTRA=ex)
	rspstr = mrdfits(self->file(), 2, mhdr, _EXTRA=ex)
	
	IF !err EQ -1 THEN BEGIN
		print, 'mrdfits failed with file', file 
		RETURN
	ENDIF
	
	self.det = strtrim( sxpar(phdr, 'DETNAM'), 2)
	
	
	PTR_FREE,self.hdrs
	self.hdrs[0] = ptr_new( ebhdr )
	self.hdrs[1] = ptr_new( mhdr )
	
	dims = SIZE(rspstr.matrix, /DIMENSIONS)
	
	ch_de = ebounds.e_max - ebounds.e_min
	ph_de = rspstr.energ_hi - rspstr.energ_lo
;	en_area = ch_de # ph_de
	
	ch_min = (ebounds.e_min)[0]
	ch_max = (ebounds.e_max)[dims[0]-1]
	ph_min = (rspstr.energ_lo)[0]
	ph_max = (rspstr.energ_hi)[dims[1]-1]
	
	evar = [ [ch_min, ph_min], [ch_max, ph_max] ]
	
	
	self->set, IJB=evar, DIM=dims, MAT=rspstr.matrix, $
		CH=[PTR_NEW([ebounds.e_min,ch_max]),PTR_NEW(ch_de)], PH=[PTR_NEW([rspstr.energ_lo,ph_max]),PTR_NEW(ph_de)]

END

FUNCTION drmdata::sxpar, extno, key

	if extno lt 1 or extno gt 2 then return, ''
	return, sxpar( *(self.hdrs[extno-1]), key)

END

PRO drmdata::transpose, MODE=mode
	*self.mat = transpose(*self.mat)
	self.dim = REVERSE(self.dim)
	self.ijb = REVERSE(self.ijb)
	self.var = NOT(self.var)
END

FUNCTION drmdata::rsp_var
	ON_ERROR,2
	RETURN,{v:self.var, ch:self.dim[self.var], ph:self.dim[~ self.var]}
END

PRO drmdata::set,FILE=file,VAR=var,DIM=dim,IJB=ijb,MAT=mat,PH=ph,CH=ch
	ON_ERROR,2
	IF N_ELEMENTS(file) NE 0 THEN self.file=file
	IF N_ELEMENTS(var) NE 0 THEN self.var=var
	IF N_ELEMENTS(dim) NE 0 THEN self.dim=dim
	IF N_ELEMENTS(ijb) NE 0 THEN self.ijb=ijb
	IF N_ELEMENTS(mat) NE 0 THEN BEGIN
		IF PTR_VALID(self.mat) THEN PTR_FREE, self.mat
		self.mat = PTR_NEW(mat,/NO_COPY)
	ENDIF
	IF N_ELEMENTS(ph) NE 0 THEN BEGIN
		i=WHERE(PTR_VALID(self.ph),c)
		IF c GT 0 THEN PTR_FREE,(self.ph)[i]
		self.ph = ph
	ENDIF
	IF N_ELEMENTS(ch) NE 0 THEN BEGIN
		i=WHERE(PTR_VALID(self.ch),c)
		IF c GT 0 THEN PTR_FREE,(self.ch)[i]
		self.ch = ch
	ENDIF
END

FUNCTION drmdata::di
	IF self.var THEN RETURN, self->ph(/DE) ELSE RETURN, self->ch(/DE)
END

FUNCTION drmdata::dj
	IF self.var THEN RETURN, self->ch(/DE) ELSE RETURN, self->ph(/DE)
END

FUNCTION drmdata::i_vals
	IF self.var THEN RETURN, self->ph(/EDGES) ELSE RETURN, self->ch(/EDGES)
END

FUNCTION drmdata::j_vals
	IF self.var THEN RETURN, self->ch(/EDGES) ELSE RETURN, self->ph(/EDGES)
END

FUNCTION compdate::num_phot_edges
	RETURN, self.dim[~self.var]+1
END

FUNCTION compdate::num_chan_edges
	RETURN, self.dim[self.var]+1
END

; Fold a differiental dN/(dA*dE) flux density spectrum through the matrix using one of the integration methods (only TRAPEZOIDAL is supported now)
; to integrate over energy in each channel.  If input spectrum is a quantity n/(cm^2*keV), 
; output is counts. Similiarly if the input spectrum is a quantity n/(s*cm^2*keV), 
; output is rates. The output is not differential.
;
; First argument is the name of function which returns df/dE for an array of input energies.  The second argument
; is an array of parameters which are passed to the function, along with the photon energies of the DRM.
; E.g., if the function being called is named 'crab', then spec = crab( [params], energies )
; If FORM1 is set, then the arguments are switched, so spec = crab( energies, params ).  The function
; should know how to parse the arguments.
; 
;
; Alternatively, with model as a null string, use the SPEC keyword to pass a discrete flux spectrum that has the same number of elements as the photon side 
; of the response matrix.  In this case the keyword DIFF sets whether the input spectrum is differential ( has partial unit 1/keV ).  If this is 
; set, the input is first converted to a number spectrum. 
; The spectrum elements must be dN/dA values, or dN/(dA*dE) with the DIFF keyword.
;
;
FUNCTION drmdata::fold, model, params, SPEC=in_spec, DIFF=isdiff, FORM1=form1, TRAPEZOIDAL=trapezoidal, _REF_EXTRA=ex
	ON_ERROR,2
	
	IF  ( n_elements(model) gt 0 && strlen(model) gt 0) && $
		( N_ELEMENTS(in_spec) EQ 0 ) THEN BEGIN
		ph_energy = self->ph(/EDGES)
		dE = self->ph(/DE)
		if keyword_set(isdiff) then begin
			if keyword_set(form1) then spec = CALL_FUNCTION(model, ph_energy, params, _EXTRA=ex) $
			else spec = CALL_FUNCTION(model, params, ph_energy, _EXTRA=ex)
				
			np_edeges = self.dim[~self.var]+1
			spec = dE/2.0 * (spec[0:np_edeges-2] + spec[1:*])
		endif else begin
		
			if keyword_set(form1) then spec = CALL_FUNCTION(model, ph_energy[1:*]-dE, ph_energy[1:*], params, _EXTRA=ex) $
			else spec = CALL_FUNCTION(model, params, ph_energy[1:*]-dE, ph_energy[1:*], _EXTRA=ex)
		
		endelse
		
	ENDIF ELSE BEGIN
		IF N_ELEMENTS(in_spec) NE self.dim[~ self.var] THEN BEGIN
			print, 'SPEC Length Must equal # of photon energies in response.'
			print, 'Call obj->rsp_var() to get the dimensions'
			RETURN,[-1]
		ENDIF
		
		if keyword_set(isdiff) then begin
			dE = self->ph(/DE)
			spec = in_spec * dE
		endif else spec = in_spec
		
	ENDELSE

	conv = (*self.mat) # spec	

	RETURN,conv
END

FUNCTION drmdata::specconv, model, params, SPEC=spec
	ON_ERROR,2
	
	IF N_ELEMENTS(spec) EQ 0 THEN BEGIN
		ph_energy = self->ph(/EDGES)
		spec = CALL_FUNCTION(model, params, ph_energy)
	ENDIF ELSE $
		IF N_ELEMENTS(spec) NE self.dim[~ self.var] THEN BEGIN
			print, 'SPEC Length Must equal # of photon energies in response.'
			print, 'Call obj->rsp_var() to get the dimensions'
			RETURN,[-1]
		ENDIF
	
	;kev2 = self->en_area()
	
	conv = (*self.mat) # spec
	RETURN,conv
END

FUNCTION drmdata::integrate, PH=phrange, CH=chrange, PH_E=ph_e, CH_E=ch_e
	mat = *self.mat
	
	vari = {ch: self.var, ph: ~ self.var, $
			ij:[PTR_NEW(),PTR_NEW()], $
			di: [0,0], ax:STRARR(2)}
	
	vari.ax[vari.ph] = 'ph'
	vari.ax[vari.ch] = 'ch'
	
	IF N_ELEMENTS(phrange) GT 0 AND N_ELEMENTS(chrange) GT 0 THEN $
		RETURN, {tot:TOTAL(mat),$
				calcdim:vari.ax[0]+'X'+vari.ax[1],$
				di:vari.di}
	
	CASE N_ELEMENTS(phrange) OF
		1: BEGIN
			phrange = INDGEN(phrange[0])
			vari.di[vari.ph] = N_ELEMENTS(phrange)
			vari.ij[vari.ph] = PTR_NEW(phrange,/NO_COPY)
			END
		2: BEGIN
			phrange = INDGEN(1 + phrange[1] - phrange[0]) + phrange[0]
			vari.di[vari.ph] = N_ELEMENTS(phrange)
			vari.ij[vari.ph] = PTR_NEW(phrange,/NO_COPY)
			END
	ENDCASE
	
	CASE N_ELEMENTS(chrange) OF
		1: BEGIN
			chrange = INDGEN(chrange[0])
			vari.di[vari.ch] = N_ELEMENTS(chrange)
			vari.ij[vari.ch] = PTR_NEW(chrange,/NO_COPY)
			END
		2: BEGIN
			chrange = INDGEN(1 + chrange[1] - chrange[0]) + chrange[0]
			vari.di[vari.ch] = N_ELEMENTS(chrange)
			vari.ij[vari.ch] = PTR_NEW(chrange,/NO_COPY)
			END
	ENDCASE
	
	vari.ax[vari.ph] = 'ph'
	vari.ax[vari.ch] = 'ch'
	
	
	r=WHERE(vari.di GT 0,result_dims)
	CASE result_dims OF
		0:RETURN, {tot:TOTAL(mat),calcdim:vari.ax[0]+'X'+vari.ax[1]}
		1: BEGIN
			tot = FLTARR(1)
			axis = vari.ax[r]
			END
	ENDCASE
			
	
	;FOR i=0,1 DO BEGIN
		
	
	sum1 = TOTAL(mat[*,j],2)
	effa = TOTAL(sum1[i])
	
END

FUNCTION drmdata::en_area
	subscr_alph = self->di()
	subscr_beta = self->dj()
	
	RETURN, subscr_alph # subscr_beta
END

FUNCTION drmdata::r_i
	RETURN, self.ijb[0,1]-self.ijb[0,0]
END

FUNCTION drmdata::r_j
	RETURN, self.ijb[1,1]-self.ijb[1,0]
END

FUNCTION drmdata::ebox_xy4 
	coordrect = [self.ijb[0,0],self.ijb[0,1],self.ijb[1,0],self.ijb[1,1] ]
	RETURN, coordrect
END

FUNCTION drmdata::var 
	RETURN,self.var 
END 

FUNCTION drmdata::ph, EDGES=edges, DE=de
	IF KEYWORD_SET(edges) THEN BEGIN	
		ph = self.ph
		RETURN, *ph[0]
	ENDIF
	IF KEYWORD_SET(de) THEN BEGIN	
		ph = self.ph
		RETURN, *ph[1]
	ENDIF
	RETURN,self.ph 
END

FUNCTION drmdata::ch, EDGES=edges, DE=de
	IF KEYWORD_SET(edges) THEN BEGIN	
		ch = self.ch
		RETURN, *ch[0]
	ENDIF
	IF KEYWORD_SET(de) THEN BEGIN	
		ch = self.ch
		RETURN, *ch[1]
	ENDIF
	RETURN,self.ch 
END

FUNCTION drmdata::dim & RETURN,self.dim & END
FUNCTION drmdata::ijb & RETURN,self.ijb & END
FUNCTION drmdata::mat & RETURN,*self.mat & END
FUNCTION drmdata::file & RETURN,self.file & END

PRO drmdata::cleanup
	PTR_FREE,self.mat
	PTR_FREE,self.ph,self.ch
	PTR_FREE,self.hdrs
	OBJ_DESTROY,self
END

PRO drmdata__define
	tmp = {drmdata, $
				file:'',$
				det: '',$
				var:0b,$
				dim:LONARR(2),$
				ijb:FLTARR(2,2),$
				mat:PTR_NEW(),$
				hdrs:PTRARR(2),$
				ph:PTRARR(2),$
				ch:PTRARR(2) $
			}
END