## Author: Georgi Tadeus <georgi.tadeus@gmail.com>
## This program is under CC-SA license.

## -*- texinfo -*-
## @deftypefn {} {} DeconpH (FILE, STIMULATIONFRAME, SEARCHRANGEFORPEAK, ENDFRAMEFORFIT)
## Deconvolves pHluorin Trace in FILE, and saves RELEASE_RATE and CUM_RELEASE in the directory of FILE
##
## @end deftypefn

function DeconpH(file, stimulationframe, peakrange, fitend)
	if nargin == 4
	[dir, name, ext, ver] = fileparts (file);	% extract directory from filename
	fulltrace = load(file);					% load file
	
	% find peak of first stimulation and location of peak (indmax)
	[maxfirststim, indmax] = max(fulltrace(stimulationframe:(stimulationframe+peakrange)));
	indmax+=(stimulationframe-1);

	% normalize (1,0), set first peak to 1; and global minimum to 0;
	% it is important for later deconvolution that there are no negative values
	normtrace = (fulltrace - min(fulltrace)) / ( maxfirststim - min(fulltrace) );

	% prepare equidistant x values
	xfit = round((linspace( 1, fitend-indmax+1,  fitend-indmax+1 ))');
	% guess coefficients
	p = [normtrace(indmax)-normtrace(fitend) 0.1 normtrace(fitend)];
	% do single exponential fit, using the following function (seperate file)
	% function y = SingleExp(x, p)
	%    y = p(1)* exp ( - p(2) * x ) + p(3);
	% endfunction
	[yfit pfit cvg iter] = leasqr(xfit, normtrace(indmax:fitend), p, "SingleExp");

	% print tau
	tau = pfit(2)
	
	N=1024;		% important for deconvolutionl powers of 2 are useful values
	% pad normtrace to N
	y = padarray(normtrace,[N-size(normtrace, 1)], 0, "Post");
	% calculate filter for deconvolution based on fit
	h=(exp(-(1:(N/2)).*tau))';

	% deconvolution
	Lx=length(y)-length(h)+1;  
	Lx2=pow2(nextpow2(Lx));   		 % Find smallest power of 2 that is > Lx
	Y=fft(y, Lx2);		  			
	H=fft(h, Lx2);		  			 
	X=Y./H;        		  			 
	x=real(ifft(X, Lx2));      		% Inverse fast Fourier transform
	x=x(1:1:Lx);               		% Take just the first N elements
	release_rate=x/max(abs(x));		% Normalize the output

	% cut off padded trails
	release_rate = release_rate(1:length(normtrace));

	% integrate release_rate
	cum_release = zeros(length(release_rate));
	cum_release = release_rate(1);

	for i = 2:length(release_rate)
		cum_release(i)=cum_release(i-1)+release_rate(i);
	end
	
	% save output
	cum_release=cum_release';
	save("-ascii", strcat(dir, "/", "release_rate.txt"), "release_rate");
	save("-ascii", strcat(dir, "/", "cum_release.txt"), "cum_release");
	
	% plot output
	subplot(2,2,1);plot(fulltrace);title('fulltrace');subplot(2,2,2);plot(normtrace);title('normtrace & fit');hold all;plot(xfit+indmax-1, yfit, '@');
	subplot(2,2,3);plot(release_rate);title('release_rate');subplot(2,2,4);plot(cum_release);title('cum_release')
  else
    help DeconpH;
  end
endfunction
