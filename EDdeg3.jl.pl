#!/usr/bin/env perl

# Get the affine or projective ED degree of the marginal independence model
# for simplicial complexes on four vertices. The vertices must all appear
# in the complex and be labeled by {1, 2, 3, 4}.
#
# TODO: generic degrees.

use Modern::Perl 2018;
use utf8::all;
use Getopt::Long qw(:config bundling);

use Path::Tiny qw(tempfile);
use Tie::Select;
use IPC::Run3;

use List::Util qw(uniq any);
use Array::Set qw(set_diff);
use Algorithm::Combinatorics qw(subsets);

GetOptions(
    'A|affine'     => \my $affine,
    'P|projective' => \my $projective,
) or die 'failed parsing options';

# Default is to run affine computation.
$affine //= 1 if not $projective;

# Saturate a list of sets (facets) to the generated simplicial complex.
sub complexify {
    my $faces = shift;
    my @N = uniq sort map { @$_ } @$faces;
    [ grep {
            my $A = $_;
            any { set_diff($A, $_)->@* == 0 } @$faces
        } subsets(\@N)
    ]
}

# Run Julia on the given file and return its output, dying on error.
sub run_julia {
    my $file = shift;

    run3 ['julia', $file], \undef, \my $stdout, \my $stderr;

    my $code = $? >> 8;
    die "Julia terminated with exit code $code and message:\n" .
        "$stderr\nSource code:\n" . $file->slurp_utf8
        if $code != 0 or not defined $stdout or not length $stdout;

    $stdout
}

# Get a substitution list for the Julia program from a simplicial complex.
# The first argument is the simplicial complex, all other arguments are
# prepended to the subslist. Use it to set "q=>1" to get the affine degree.
sub subslist {
    my $complex = shift;

    join ', ', @_, map {
        my @A = sort @$_;
        my $name = join('', @A);
        'q'  . join('', @A) .
        '=>' . join('*', map "q$_", @A)
    } grep { $_->@* >= 2 } @$complex
}

# Get the affine and projective ED degree of the given complex.
sub ed_degree {
    use Data::Handle; # reusable <DATA>
    my $file = tempfile;
    my $data = Data::Handle->new(__PACKAGE__);
    $file->append_utf8(
        map { s/__SUBSLIST__/subslist(@_)/gre } <$data>
    );

    my $output = run_julia $file;
    my ($general, $real, $sample) = $output =~ /(\d+) (\d+) (.*)/;
    { ed_degree => $general, real_solutions => $real }
}

my @input = map { chomp; $_ } <<>>;
for my $item (@input) {
    chomp $item;
    my $complex = complexify eval($item =~ tr/{}/[]/r);

    my @msg;
    if ($affine) {
        my $res = ed_degree($complex, 'q=>1');
        push @msg, 'affine(' .
            join(', ', $res->@{'ed_degree', 'real_solutions'}) .
        ')';
    }
    if ($projective) {
        my $res = ed_degree($complex);
        push @msg, 'projective(' .
            join(', ', $res->@{'ed_degree', 'real_solutions'}) .
        ')';
    }
    say $item, ': ', join(', ', @msg);
}

__DATA__
using HomotopyContinuation

@var q q1 q2 q3 q12 q13 q23 q123

sample = randn(8)
p111,p112,p121,p122,
p211,p212,p221,p222 = sample

diffs = [
  -p111 + q123,
  -p112 + q12-q123,
  -p121 + q13-q123,
  -p122 + q23-q123,
  -p211 + q2-q12-q23+q123,
  -p212 + q3-q13-q23+q123,
  -p221 + q1-q12-q13+q123,
  -p222 + q-q1-q2+q12-q3+q13+q23-q123
]

dist = sum([d^2 for d in diffs])
subslist = [ __SUBSLIST__ ]
if length(subslist) > 0
  dist = subs(dist, subslist...)
end

vars = variables(dist)
eqns = differentiate(dist, vars)

R = solve(eqns, show_progress = false)
C = certify(eqns, R)
println(ndistinct_certified(C), " ", ndistinct_real_certified(C), " ", sample)
