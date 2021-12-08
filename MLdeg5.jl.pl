#!/usr/bin/env perl

# Get the affine or projective ML degree of the marginal independence model
# for simplicial complexes on four vertices. The vertices must all appear
# in the complex and be labeled by {1, 2, 3, 4, 5}.
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

@var q1 q2 q3 q4 q5 q12 q13 q14 q15 q23 q24 q25 q34 q35 q45 q123 q124 q125 q134 q135 q145 q234 q235 q245 q345 q1234 q1235 q1245 q1345 q2345 q12345 z u[1:32]

p = [
  z*(q12345),
  z*(q1345-q12345),
  z*(q1235-q12345),
  z*(q2345-q12345),
  z*(q1234-q12345),
  z*(q1245-q12345),
  z*(q135-q1235-q1345+q12345),
  z*(q124-q1234-q1245+q12345),
  z*(q234-q1234-q2345+q12345),
  z*(q134-q1234-q1345+q12345),
  z*(q123-q1234-q1235+q12345),
  z*(q145-q1245-q1345+q12345),
  z*(q125-q1235-q1245+q12345),
  z*(q235-q1235-q2345+q12345),
  z*(q245-q1245-q2345+q12345),
  z*(q345-q1345-q2345+q12345),
  z*(q35-q135-q235+q1235-q345+q1345+q2345-q12345),
  z*(q34-q134-q234+q1234-q345+q1345+q2345-q12345),
  z*(q45-q145-q245+q1245-q345+q1345+q2345-q12345),
  z*(q15-q125-q135+q1235-q145+q1245+q1345-q12345),
  z*(q25-q125-q235+q1235-q245+q1245+q2345-q12345),
  z*(q14-q124-q134+q1234-q145+q1245+q1345-q12345),
  z*(q13-q123-q134+q1234-q135+q1235+q1345-q12345),
  z*(q12-q123-q124+q1234-q125+q1235+q1245-q12345),
  z*(q23-q123-q234+q1234-q235+q1235+q2345-q12345),
  z*(q24-q124-q234+q1234-q245+q1245+q2345-q12345),
  z*(q3-q13-q23+q123-q34+q134+q234-q1234-q35+q135+q235-q1235+q345-q1345-q2345+q12345),
  z*(q2-q12-q23+q123-q24+q124+q234-q1234-q25+q125+q235-q1235+q245-q1245-q2345+q12345),
  z*(q1-q12-q13+q123-q14+q124+q134-q1234-q15+q125+q135-q1235+q145-q1245-q1345+q12345),
  z*(q4-q14-q24+q124-q34+q134+q234-q1234-q45+q145+q245-q1245+q345-q1345-q2345+q12345),
  z*(q5-q15-q25+q125-q35+q135+q235-q1235-q45+q145+q245-q1245+q345-q1345-q2345+q12345),
  z*(1-q1-q2+q12-q3+q13+q23-q123-q4+q14+q24-q124+q34-q134-q234+q1234-q5+q15+q25-q125+q35-q135-q235+q1235+q45-q145-q245+q1245-q345+q1345+q2345-q12345)
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
tparams = randn(length(u))

R = solve(eqns, solutions(MR);
  start_parameters  = sparams,
  target_parameters = tparams
)

C = certify(eqns, R, tparams)
println(ndistinct_certified(C), " ", ndistinct_real_certified(C), " ", tparams)
