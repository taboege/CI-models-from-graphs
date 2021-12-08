#!/usr/bin/env perl

# Get the affine or projective ML degree of the marginal independence model
# for simplicial complexes on four vertices. The vertices must all appear
# in the complex and be labeled by {1, 2, 3}.
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

# Get the affine and projective ML degree of the given complex.
sub ml_degree {
    use Data::Handle; # reusable <DATA>
    my $file = tempfile;
    my $data = Data::Handle->new(__PACKAGE__);
    $file->append_utf8(
        map { s/__SUBSLIST__/subslist(@_)/gre } <$data>
    );

    my $output = run_julia $file;
    my ($general, $real, $sample) = $output =~ /(\d+) (\d+) (.*)/;
    { ml_degree => $general, real_solutions => $real }
}

my @input = map { chomp; $_ } <<>>;
for my $item (@input) {
    chomp $item;
    my $complex = complexify eval($item =~ tr/{}/[]/r);

    my @msg;
    if ($affine) {
        my $res = ml_degree($complex, 'z=>1');
        push @msg, 'affine(' .
            join(', ', $res->@{'ml_degree', 'real_solutions'}) .
        ')';
    }
    #if ($projective) {
    #    my $res = ml_degree($complex);
    #    push @msg, 'projective(' .
    #        join(', ', $res->@{'ml_degree', 'real_solutions'}) .
    #    ')';
    #}
    say $item, ': ', join(', ', @msg);
}

__DATA__
using HomotopyContinuation

@var q1 q2 q3 q12 q13 q23 q123 z u[1:8]

p = [
  z*(q123),
  z*(q12-q123),
  z*(q13-q123),
  z*(q23-q123),
  z*(q2-q12-q23+q123),
  z*(q3-q13-q23+q123),
  z*(q1-q12-q13+q123),
  z*(1-q1-q2+q12-q3+q13+q23-q123)
]

L = sum([u[i] * log(p[i]) for i = 1:length(u)])
subslist = [ __SUBSLIST__ ]
if length(subslist) > 0
  L = subs(L, subslist...)
end

vars = setdiff(variables(L), u)
eqns = System(differentiate(L, vars), parameters = u)

MR = monodromy_solve(eqns)

sparams = parameters(MR)
tparams = rand(1:1000, length(u))

R = solve(eqns, solutions(MR);
  start_parameters  = sparams,
  target_parameters = tparams
)

C = certify(eqns, R, tparams)
println(ndistinct_certified(C), " ", ndistinct_real_certified(C), " ", tparams)
