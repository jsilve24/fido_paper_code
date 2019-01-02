use strict;
use warnings;

# -------------------------------------------------------------------------------------------------------------------------
# generate SLURM jobs to test algorithmic efficiency (as seconds per effective sample size); each job is a distinct
# set of dimensions (N,D,Q); probably overkill to submit each as a separate job but it's def faster
# -------------------------------------------------------------------------------------------------------------------------

my @N_vals;
my @D_vals;
my @Q_vals;

if($#ARGV < 2) {
	print("Usage: perl generate_runtime_jobs.pl {vary} {cores} {MAP-only}\n");
}

my $vary = $ARGV[0];
my $cores = $ARGV[1];
my $use_MAP = "F";
if($ARGV[2] != 0) {
	$use_MAP = "T";
}
my $optim = 'lbfgs';
my $samples = 2000;

open(my $job_listing, '>>', 'job_listing_'.$vary.'.txt');

if($vary eq 'test') {
	# trivial case to test script generation
	@N_vals = qw(1000);
	@D_vals = qw(1000);
	@Q_vals = qw(500);
} elsif($vary eq 'N') {
	# varying N
	@N_vals = qw(1 3 5 10 20 30 50 100 250 500 750 1000);
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

# available methods: me mc sc su clm svbcm svbcf svbum svbuf
my @methods = qw(mc);

my @rseed = qw(1 2 3);

print("Generating ".(($#N_vals+1)*($#D_vals+1)*($#Q_vals+1)*($#methods+1)*($#rseed+1))." slurm scripts...\n");

my @usable_nodes = qw(x2-01-3 x2-04-4 x2-05-1 x2-06-3 x2-07-2 x2-01-2 x2-01-4 x2-02-1 x2-02-2 x2-02-3 x2-02-4 x2-03-1 x2-03-2 x2-03-3 x2-03-4 x2-04-1 x2-04-2 x2-04-3 x2-05-2 x2-05-3 x2-05-4 x2-06-1 x2-06-2 x2-06-4 x2-07-1x2-07-3 x2-07-4 x2-08-1 x2-08-2 x2-08-3 x2-08-4);

my $use_MKL = 1;

my $simulation_ret = -1;
my $sys_response = '';
my $rand_n = -1;
my $node_found = 0;

my $file_suffix = "";
if($use_MKL) {
	$file_suffix = $file_suffix."_MKL";
}
$file_suffix = $file_suffix."_".$cores;

for my $N (@N_vals) {
	for my $D (@D_vals) {
		for my $Q (@Q_vals) {
			for my $rep (@rseed) {

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

					# extra arguments will be ignored if not needed
					print $fh 'srun Rscript simulate_efficiency.R '.$N.' '.$D.' '.$Q.' '.$rep.' '.$m_idx.' '.$samples.' '.$use_MAP.' '.$optim.' 50000 1e-10 0.004 0.99 fitted_models_2018-12-30 '.$file_suffix."\n";

					close $fh;

					# find a usable node
					$node_found = 0;
					while(!$node_found) {
						$rand_n = int(rand($#usable_nodes+1));
						print("Trying node $usable_nodes[$rand_n]...\n");
						$sys_response = `sinfo -N $usable_nodes[$rand_n] -o "%N,%C,%e" | grep $usable_nodes[$rand_n]`;
						print("\t".$sys_response);
						# e.g.: x2-01-3,16/12/0/28,217480
						if($sys_response =~ /^.*?,\d+\/(\d+)\/.*,(\d+)$/) {
						        if($1 >= $cores && $2 >= 64000) {
								$node_found = 1;
						        }
						}
					}
					if($vary ne 'test') {
						my $call_str = "sbatch --ntasks=1 --cpus-per-task=$cores --nodelist=$usable_nodes[$rand_n] $filename";
						print("Calling: ".$call_str."\n");
						$sys_response = `$call_str`;
						sleep(1);
						chomp($sys_response);
						print $job_listing substr($sys_response,20,length($sys_response))."\tmodel=".uc($m_idx)."\tN=".$N."\tD=".$D."\tQ=".$Q."\tR=".$rep."\n";
					}
				}
			}
		}
	}
}


