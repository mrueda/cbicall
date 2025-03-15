package DNAseq;

use strict;
use warnings;
use Carp;
use File::Spec::Functions qw(catfile);
use feature               qw(say);

sub new {

    my ( $class, $self ) = @_;
    bless $self, $class;
    return $self;
}

sub variant_calling {

    my $self          = shift;
    my $pipeline      = $self->{pipeline};
    my $workflow_mode = $self->{workflow_mode};
    my $dir           = $self->{projectdir};
    my $threads       = $self->{threads};
    my $mode          = $self->{mode};
    my $id            = $self->{id};
    my $debug         = $self->{debug};
    my $cmd;
    my $log;

    if ( $workflow_mode eq 'bash' ) {
        my $bash_str = 'bash_' . $pipeline . '_' . $mode;
        my $bash     = $self->{$bash_str};
        $log = $bash_str . '.log';
        $cmd = "cd $dir && $bash -t $threads > $log 2>&1";
    }
    elsif ( $workflow_mode eq 'snakemake' ) {

        # Example snakemake command; adjust as needed.
        my $smk_str = 'smk_' . $pipeline . '_' . $mode;
        my $smk     = $self->{$smk_str};
        $log = $smk_str . '.log';
        my $snakemake_cmd = 'snakemake --forceall all';
        $cmd = "cd $dir && $snakemake_cmd -s $smk --cores $threads > $log 2>&1";
    }
    else {
        die
"Invalid workflow_mode: $workflow_mode. Only 'bash' and 'snakemake' are implemented.\n";
    }

    submit_cmd( $cmd, $dir, $log, $id, $debug );
    return 1;
}

sub submit_cmd {

    my ( $cmd, $dir, $log, $id, $debug ) = @_;
    my $path_log = catfile( $dir, $log );
    my $msg      = "Failed to execute: $id\nPlease check this file:\n$path_log\n";
    {
        local $SIG{__DIE__} = 'DEFAULT';
        system($cmd) == 0 or ( $debug ? confess($msg) : die($msg) );
    }
    return 1;
}

1;
