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
# prepended to the subslist. Use it to set "z=>1" to get the affine degree.
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
        my $res = ed_degree($complex, 'z=>1');
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

@var q1 q2 q3 q4 q12 q13 q14 q23 q24 q34 q123 q124 q134 q234 q1234 z

sample = randn(16)
p0000,p0001,p0010,p0011,p0100,p0101,p0110,p0111,
p1000,p1001,p1010,p1011,p1100,p1101,p1110,p1111 = sample

diffs = [
  -p0000 + z*(q1234),
  -p1000 + z*(q234-q1234),
  -p0010 + z*(q124-q1234),
  -p0100 + z*(q134-q1234),
  -p0001 + z*(q123-q1234),
  -p0110 + z*(q14-q124-q134+q1234),
  -p0101 + z*(q13-q123-q134+q1234),
  -p0011 + z*(q12-q123-q124+q1234),
  -p1001 + z*(q23-q123-q234+q1234),
  -p1010 + z*(q24-q124-q234+q1234),
  -p1100 + z*(q34-q134-q234+q1234),
  -p1101 + z*(q3-q13-q23+q123-q34+q134+q234-q1234),
  -p1011 + z*(q2-q12-q23+q123-q24+q124+q234-q1234),
  -p0111 + z*(q1-q12-q13+q123-q14+q124+q134-q1234),
  -p1110 + z*(q4-q14-q24+q124-q34+q134+q234-q1234),
  -p1111 + z*(1-q1-q2+q12-q3+q13+q23-q123-q4+q14+q24-q124+q34-q134-q234+q1234)
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
