package Config;

use strict;
use warnings;
use autodie;
use Cwd qw(abs_path);
use Sys::Hostname;
use File::Spec::Functions qw(catdir catfile);
use YAML::XS              qw(LoadFile);
use List::Util            qw(all);
use Types::Standard       qw(Str Int);
use Type::Utils           qw(enum);

# Define custom enum types for allowed parameter values
my $mode_type          = enum [qw(single cohort)];
my $pipeline_type      = enum [qw(wes wgs mit)];
my $genome_type        = enum [qw(hg19 hg38)];
my $organism_type      = enum [ 'Homo Sapiens',       'Mus musculus' ];
my $technology_type    = enum [ 'Illumina HiSeq',     'NovaSeq' ];
my $capture_type       = enum [ 'Agilent SureSelect', 'NimbleGen' ];
my $workflow_mode_type = enum [qw(bash nextflow snakemake)];

# Map each parameter key to its type constraint
my %param_types = (
    mode          => $mode_type,
    pipeline      => $pipeline_type,
    genome        => $genome_type,
    organism      => $organism_type,
    technology    => $technology_type,
    capture       => $capture_type,
    workflow_mode => $workflow_mode_type,
);

# Define default values
my %default = (
    mode          => 'single',
    sample        => undef,
    pipeline      => 'wes',
    genome        => 'hg19',
    organism      => 'Homo Sapiens',
    technology    => 'Illumina HiSeq',
    capture       => 'Agilent SureSelect',
    workflow_mode => 'bash'
);

sub read_param_file {
    my $yaml_file = shift;
    my $param     = LoadFile($yaml_file);

    # Merge provided parameters with defaults, and validate allowed values
    foreach my $key ( keys %$param ) {
        if ( exists $default{$key} ) {
            my $value = $param->{$key};

            # Only validate if an allowed type exists and value is defined
            if ( exists $param_types{$key} && defined $value ) {

                # This will die with a detailed error if validation fails
                $param_types{$key}->assert_valid($value);
            }
            $default{$key} = $value;
        }
        else {
            die "Parameter '$key' does not exist (typo?)\n";
        }
    }

    return wantarray ? %default : \%default;
}

sub set_config_values {
    my $param = shift;
    my $user  = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
    my $workflows_bash_dir =
      abs_path( catdir( $main::Bin, '..', 'workflows', 'bash' ) );
    my $workflows_snakemake_dir =
      abs_path( catdir( $main::Bin, '..', 'workflows', 'snakemake' ) );

    # We load %config with the default values
    my %config = (
        user            => $user,
        bash_parameters => catfile( $workflows_bash_dir, 'parameters.sh' ),
        bash_wes_single => catfile( $workflows_bash_dir, 'wes_single.sh' ),
        bash_wes_cohort => catfile( $workflows_bash_dir, 'wes_cohort.sh' ),
        bash_mit_single => catfile( $workflows_bash_dir, 'mit_single.sh' ),
        bash_mit_cohort => catfile( $workflows_bash_dir, 'mit_cohort.sh' ),
        bash_coverage   => catfile( $workflows_bash_dir, 'coverage.sh' ),
        bash_jaccard    => catfile( $workflows_bash_dir, 'jaccard.sh' ),
        bash_vcf2sex    => catfile( $workflows_bash_dir, 'vcf2sex.sh' ),
        smk_wes_single => catfile( $workflows_snakemake_dir, 'wes_single.smk' ),
        smk_wes_cohort => catfile( $workflows_snakemake_dir, 'wes_cohort.smk' ),
        smk_mit_single => catfile( $workflows_snakemake_dir, 'mit_single.smk' ),
        smk_mit_cohort => catfile( $workflows_snakemake_dir, 'mit_cohort.smk' ),
        smk_config     => catfile( $workflows_snakemake_dir, 'config.yaml' )
    );

    # Below are a few internal configaters that do not have (or we don't allow for) default values
    $config{id}   = time . substr( "00000$$", -5 );
    $config{date} = localtime();
    my $tmp_str = '/'
      . 'cbicall' . '_'
      . $param->{workflow_mode} . '_'
      . $param->{pipeline} . '_'
      . $param->{mode} . '_'
      . $config{id};    # User will make symbolic link to final folder
    $config{projectdir} = catdir( abs_path( $param->{sample} ), $tmp_str );
    my @tmp = split /\//, $param->{sample};
    $config{output_basename} = $tmp[-1];
    $config{hostname}        = hostname;
    $config{user}            = $user;
    chomp( my $threadshost = qx{/usr/bin/nproc} ) // 1;
    $config{threadshost} = 0 + $threadshost;                          # coercing it to be a number
    $config{threadsless} = $threadshost > 1 ? $threadshost - 1 : 1;
    my $str_threadsless = $config{threadsless};                       # We copy it (otherwise it will get "stringified" below and printed with "" in log.json)
    $config{zip} =
      ( -x '/usr/bin/pigz' )
      ? "/usr/bin/pigz -p $str_threadsless"
      : '/bin/gunzip';

    # Determine the architecture
    my $uname = `uname -m`;
    chomp($uname);
    my $arch =
      ( $uname eq 'x86_64' )
      ? 'x86_64'
      : ( $uname eq 'aarch64' ? 'arm64' : $uname );

    # Load arch to config
    $config{arch} = $arch;

    # Check if all required bash files exist and have +x permission
    die "You don't have +x permission for one or more Bash files"
      unless all { -x $config{$_} }
      qw(
      bash_parameters bash_wes_single bash_wes_cohort
      bash_mit_single bash_mit_cohort bash_coverage bash_jaccard bash_vcf2sex
      );

    return wantarray ? %config : \%config;

}
1;
