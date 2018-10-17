2018/10/17

NOTES:

"Run 2" used optimization parameters: step_size=0.002 max_iter=50000 b1=0.99 eps_f=1e-10

Stan run with "--ntasks=1 --cpus-per-task=4"

Mongrel run with "--ntasks=1 --cpus-per-task=1"

Mongrel eigendecomp and Cholesky fail at largest dimensions of D attempted (N = 100, D >= 250, Q = 5) with "negative length vectors are not allowed" error in optimMongrelCollapsed
	"negative length" sounds like overflow/memory issue?
	allocated 64GB RAM, could try at 128?
	this was largest N*(D-1) attempted

Stan gives treedepth warnings at a few largish configuration but doens't seem like a big issue
	JOBID	N	D	Q
	3644251	100	75	5
	3644257	100	250	5
	3644278	100	30	20
	3644281	100	30	50
	3644284	100	30	75
	3644287	100	30	100
	3644290	100	30	250
	3644293	100	30	500

random seed fixed: Stan, Mongrel get same input matrices for each N/D/Q configuration

PLOTS:

*_SpES_plot.png
	seconds per effective sample size; seconds for Stan = execution time of longest running chain

*_MSE_plot.png
	mean squared error of each Lambda sample with its corresponding Lambda_true value

*_95CI_plot.png
	percent of Lambda_true values not within posterior (0.025, 0.975) quantiles

*_zero_plot.png
	percent zero counts in Y as a function of increasing N, D, Q
