#!/usr/bin/env perl

# Get the affine or projective ED degree of the marginal independence model
# for simplicial complexes on four vertices. The vertices must all appear
# in the complex and be labeled by {1, 2, 3}.

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
    my ($general) = $output =~ /^(\d+)/;
    { ed_degree => $general }
}

my @input = map { chomp; $_ } <<>>;
for my $item (@input) {
    chomp $item;
    my $complex = complexify eval($item =~ tr/{}/[]/r);

    my @msg;
    if ($affine) {
        my $res = ed_degree($complex, 'z=>1');
        push @msg, "affine(@{[ $res->{ed_degree} ]})";
    }
    if ($projective) {
        my $res = ed_degree($complex);
        push @msg, "projective(@{[ $res->{ed_degree} ]})";
    }
    say $item, ': ', join(', ', @msg);
}

__DATA__
using HomotopyContinuation

@var q1 q2 q3 q12 q13 q23 q123 z

sample = randn(ComplexF64, 8)
p000,p001,p010,p011,
p100,p101,p110,p111 = sample

diffs = [
  -p000 + z*(q123),
  -p001 + z*(q12-q123),
  -p010 + z*(q13-q123),
  -p100 + z*(q23-q123),
  -p101 + z*(q2-q12-q23+q123),
  -p110 + z*(q3-q13-q23+q123),
  -p011 + z*(q1-q12-q13+q123),
  -p111 + z*(1-q1-q2+q12-q3+q13+q23-q123)
]

dist = sum([d^2 for d in diffs])
subslist = [ __SUBSLIST__ ]
if length(subslist) > 0
  dist = subs(dist, subslist...)
end

vars = variables(dist)
eqns = differentiate(dist, vars)

R = solve(eqns; show_progress = false,
  tracker_options = TrackerOptions(
    automatic_differentiation = 3,
    parameters = :conservative
  )
)

C = certify(eqns, R; show_progress = false)
println(ndistinct_certified(C))
