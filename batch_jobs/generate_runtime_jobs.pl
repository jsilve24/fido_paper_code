use strict;
use warnings;

# -------------------------------------------------------------------------------------------------------------------------
# generate SLURM jobs to test algorithmic efficiency (as seconds per effective sample size); each job is a distinct
# set of dimensions (N,D,Q); probably overkill to submit each as a separate job but it's def faster
# -------------------------------------------------------------------------------------------------------------------------

my @N_vals;
my @D_vals;
my @Q_vals;

my $vary = 'N';

if($vary eq 'test') {
	# trivial test case
	@N_vals = qw(10);
	@D_vals = qw(10);
	@Q_vals = qw(5);
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
	@Q_vals = qw(1 2 3 5 10 20 50 75 100 250 500);
}

my @methods = qw(1 2 3);
# 1 : Stan
# 2 : Mongrel (eigendecomposition)
# 3 : Mongrel (Cholesky decomposition)

print("Generating ".(($#N_vals+1)*($#D_vals+1)*($#Q_vals+1)*($#methods+1))." slurm scripts...\n");

for my $N (@N_vals) {
	for my $D (@D_vals) {
		for my $Q (@Q_vals) {
			for my $m_idx (@methods) {
				my $m = 'stan';
				my $m_short = 'S';
				if($m_idx == 2) {
					$m = 'mongrel_eigen';
					$m_short = 'ME';
				} elsif($m_idx == 3) {
					$m = 'mongrel_cholesky';
					$m_short = 'MC';
				}

				my $filename = 'scripts/'.$vary.'-varying_N'.$N.'_D'.$D.'_Q'.$Q.'_'.$m_short.'.slurm';
				my $logfile = 'run_'.$vary.'-varying.log';

				open(my $fh, '>', $filename);

				print $fh '#!/bin/bash'."\n";
				print $fh '#SBATCH -J '.$m_short.'_'.$N.'_'.$D.'_'.$Q."\n";
				print $fh '#SBATCH --mem=64GB'."\n";
				print $fh '#SBATCH --get-user-env'."\n";
				print $fh '#SBATCH --time=48:00:00'."\n";
				print $fh '#'."\n\n";

				print $fh 'module add R/3.4.2-fasrc01'."\n";
				print $fh 'module add gcc/5.3.0-fasrc01'."\n\n";

				print $fh 'cd /data/mukherjeelab/Mongrel/mongrel_paper_code'."\n\n";

				print $fh 'srun Rscript simulate_efficiency.R '.$N.' '.$D.' '.$Q.' '.$m_idx.' '.$logfile.' 0.002 50000 0.99 1e-10'."\n";

				close $fh;

				if($m_idx < 2) {
					`sbatch --ntasks=1 --cpus-per-task=4 $filename`;
				} else {
					`sbatch --ntasks=1 --cpus-per-task=1 $filename`;
				}
				sleep(1);
				print("Submitting $filename...\n");
			}

		}
	}
}

