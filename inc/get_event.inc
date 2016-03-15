	integer nbins_basic
	real*8 x_basic(0:maxnbin)
	common /basicbin/ x_basic,nbins_basic

	real*8 fxsec_CCQE,fxsec_CCRes,fxsec_pirs,fxsec_pico,fxsec_r
	common /fxsec/ fxsec_CCQE,fxsec_CCRes,fxsec_pirs,fxsec_pico,fxsec_r