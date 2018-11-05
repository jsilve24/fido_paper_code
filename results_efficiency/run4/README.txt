Default optimization parameters on these runs:
	step_size = 0.002
	max_iter = 50000
	b1 = 0.99
	eps_f = 1e-10

These cases had to be run with INCREASED step_size (=0.004)

	model=ME N=100 D=30 Q=100 R=2
	model=ME N=100 D=30 Q=75 R=2

	model=MC N=100 D=30 Q=100 R=2
	model=MC N=100 D=30 Q=75 R=2

	model=MC N=100 D=500 Q=5 R=1
	model=MC N=100 D=500 Q=5 R=2
	model=MC N=100 D=500 Q=5 R=3

These Mongrel cases DID NOT FINISH in 48 hrs.!

	model=ME N=100 D=500 Q=5 R=1
	model=ME N=100 D=500 Q=5 R=2
	model=ME N=100 D=500 Q=5 R=3



Interesting things in N runs:

(1) Stan collapsed doesn't finish any runs at N=1000 D=30 Q=5 in the 48-hour limit but Stan
    uncollapsed /does/

Interesting things in D runs:

(1) Stan uncollapsed doesn't finish any runs at N=100 D=250 Q=5 or N=100 D=500 Q=5 in the
    48-hour limit
