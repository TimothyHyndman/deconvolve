#Fourier transform of commonly used kernel function in deconvolution. (a default setting for the input)

phiK2 <-function(tt)
{
	y=(1-tt^2)^3
	y[abs(tt)>1]=0
	y
}
