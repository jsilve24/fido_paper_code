use strict;
use warnings;

# -------------------------------------------------------------------------------------------------------------------------
# generate SLURM jobs to test algorithmic efficiency (as seconds per effective sample size); each job is a distinct
# set of dimensions (N,D,Q); probably overkill to submit each as a separate job but it's def faster
# -------------------------------------------------------------------------------------------------------------------------

my @N_vals;
my @D_vals;
my @Q_vals;

my $vary = 'Q';

open(my $job_listing, '>>', 'job_listing_'.$vary.'.txt');

if($vary eq 'test') {
	# trivial case to test script generation
	@N_vals = qw(1000);
	@D_vals = qw(1000);
	@Q_vals = qw(500);
} elsif($vary eq 'N') {
	# varying N
	@N_vals = qw(10 20 30 50 100 250 500 750 1000);
	@D_vals = qw(30);
	@Q_vals = qw(5);
} elsif($vary eq 'D') {
	# varying D
	@N_vals = qw(100);
	@D_vals = qw(3 5 10 25 50 75 100 250 500);
	@Q_vals = qw(5);
} elsif($vary eq 'Q') {
	# varying Q
	@N_vals = qw(100);
	@D_vals = qw(30);
	@Q_vals = qw(2 4 10 20 50 75 100 250 500);
}

#my @methods = qw(me mc mcp sc su);
my @methods = qw(su);
# 1 : Mongrel (eigendecomposition)
# 2 : Mongrel (Cholesky)
# 3 : Mongrel (Cholesky, partial)
# 4 : Stan (collapsed)
# 5 : Stan (uncollapsed)
# 6 : naive (conjugate linear model)

my $replicates = 1;

print("Generating ".(($#N_vals+1)*($#D_vals+1)*($#Q_vals+1)*($#methods+1)*$replicates)." slurm scripts...\n");

my $simulation_ret = -1;
my $sys_response = '';
for my $N (@N_vals) {
	for my $D (@D_vals) {
		for my $Q (@Q_vals) {

			for (my $rep=1; $rep <= $replicates; $rep++) {
				# separately generate data; eventually we'll want more replicates
				chdir '..';
				my $simulation_ret = system("Rscript generate_data.R $N $D $Q $rep simulated_data");
				chdir 'batch_jobs';
				if($simulation_ret != 0) {
					print("Problem generating simulated data for $N, $D, $Q case (varying $vary)!\n");
					exit(1);
				}

				for my $m_idx (@methods) {

					my $filename = 'scripts/'.$vary.'-varying_N'.$N.'_D'.$D.'_Q'.$Q.'_'.$m_idx.'.slurm';
					my $logfile = 'run_'.$vary.'-varying_2018-10-25.log';

					open(my $fh, '>', $filename);

					print $fh '#!/bin/bash'."\n";
					print $fh '#SBATCH -J '.uc($m_idx).'_'.$N.'_'.$D.'_'.$Q."\n";
					print $fh '#SBATCH --mem=64GB'."\n";
					print $fh '#SBATCH --get-user-env'."\n";
					print $fh '#SBATCH --time=48:00:00'."\n";
					print $fh '#'."\n\n";

					print $fh 'module add R/3.4.2-fasrc01'."\n";
					print $fh 'module add gcc/5.3.0-fasrc01'."\n\n";

					print $fh 'cd /data/mukherjeelab/Mongrel/mongrel_paper_code'."\n\n";

					print $fh 'srun Rscript simulate_efficiency.R '.$N.' '.$D.' '.$Q.' '.$rep.' '.$m_idx.' '.$logfile.' 0.002 50000 0.99 1e-10'."\n";

					close $fh;

					if($vary ne 'test') {
						if($m_idx eq 'sc' || $m_idx eq 'su') {
							$sys_response = `sbatch --ntasks=1 --cpus-per-task=4 $filename`;
						} else {
							$sys_response = `sbatch --ntasks=1 --cpus-per-task=1 $filename`;
						}
						sleep(1);
						chomp($sys_response);
						print $job_listing substr($sys_response,20,length($sys_response))."\tmodel=".uc($m_idx)."\tN=".$N."\tD=".$D."\tQ=".$Q."\n";
					}
				}
			}

		}
	}
}

