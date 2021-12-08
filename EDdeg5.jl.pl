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

@var q q1 q2 q3 q4 q5 q12 q13 q14 q15 q23 q24 q25 q34 q35 q45 q123 q124 q125 q134 q135 q145 q234 q235 q245 q345 q1234 q1235 q1245 q1345 q2345 q12345

sample = randn(32)
p11111,p11112,p11121,p11122, p11211,p11212,p11221,p11222,
p12111,p12112,p12121,p12122, p12211,p12212,p12221,p12222,
p21111,p21112,p21121,p21122, p21211,p21212,p21221,p21222,
p22111,p22112,p22121,p22122, p22211,p22212,p22221,p22222 = sample

diffs = [
  -p11111 + q12345
  -p11112 + q1345-q12345
  -p11121 + q1235-q12345
  -p11122 + q2345-q12345
  -p11211 + q1234-q12345
  -p11212 + q1245-q12345
  -p11221 + q135-q1235-q1345+q12345
  -p11222 + q124-q1234-q1245+q12345
  -p12111 + q234-q1234-q2345+q12345
  -p12112 + q134-q1234-q1345+q12345
  -p12121 + q123-q1234-q1235+q12345
  -p12122 + q145-q1245-q1345+q12345
  -p12211 + q125-q1235-q1245+q12345
  -p12212 + q235-q1235-q2345+q12345
  -p12221 + q245-q1245-q2345+q12345
  -p12222 + q345-q1345-q2345+q12345
  -p21111 + q35-q135-q235+q1235-q345+q1345+q2345-q12345
  -p21112 + q34-q134-q234+q1234-q345+q1345+q2345-q12345
  -p21121 + q45-q145-q245+q1245-q345+q1345+q2345-q12345
  -p21122 + q15-q125-q135+q1235-q145+q1245+q1345-q12345
  -p21211 + q25-q125-q235+q1235-q245+q1245+q2345-q12345
  -p21212 + q14-q124-q134+q1234-q145+q1245+q1345-q12345
  -p21221 + q13-q123-q134+q1234-q135+q1235+q1345-q12345
  -p21222 + q12-q123-q124+q1234-q125+q1235+q1245-q12345
  -p22111 + q23-q123-q234+q1234-q235+q1235+q2345-q12345
  -p22112 + q24-q124-q234+q1234-q245+q1245+q2345-q12345
  -p22121 + q3-q13-q23+q123-q34+q134+q234-q1234-q35+q135+q235-q1235+q345-q1345-q2345+q12345
  -p22122 + q2-q12-q23+q123-q24+q124+q234-q1234-q25+q125+q235-q1235+q245-q1245-q2345+q12345
  -p22211 + q1-q12-q13+q123-q14+q124+q134-q1234-q15+q125+q135-q1235+q145-q1245-q1345+q12345
  -p22212 + q4-q14-q24+q124-q34+q134+q234-q1234-q45+q145+q245-q1245+q345-q1345-q2345+q12345
  -p22221 + q5-q15-q25+q125-q35+q135+q235-q1235-q45+q145+q245-q1245+q345-q1345-q2345+q12345
  -p22222 + q-q1-q2+q12-q3+q13+q23-q123-q4+q14+q24-q124+q34-q134-q234+q1234-q5+q15+q25-q125+q35-q135-q235+q1235+q45-q145-q245+q1245-q345+q1345+q2345-q12345
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
