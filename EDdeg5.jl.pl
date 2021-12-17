#!/usr/bin/env perl

# Get the affine or projective ED degree of the marginal independence model
# for simplicial complexes on four vertices. The vertices must all appear
# in the complex and be labeled by {1, 2, 3, 4, 5}.

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

@var q1 q2 q3 q4 q5 q12 q13 q14 q15 q23 q24 q25 q34 q35 q45 q123 q124 q125 q134 q135 q145 q234 q235 q245 q345 q1234 q1235 q1245 q1345 q2345 q12345 z

sample = randn(ComplexF64, 32)
p00000,p00001,p00010,p00011, p00100,p00101,p00110,p00111,
p01000,p01001,p01010,p01011, p01100,p01101,p01110,p01111,
p10000,p10001,p10010,p10011, p10100,p10101,p10110,p10111,
p11000,p11001,p11010,p11011, p11100,p11101,p11110,p11111 = sample

diffs = [
  -p00000 + z*(q12345),
  -p01000 + z*(q1345-q12345),
  -p00010 + z*(q1235-q12345),
  -p10000 + z*(q2345-q12345),
  -p00001 + z*(q1234-q12345),
  -p00100 + z*(q1245-q12345),
  -p01010 + z*(q135-q1235-q1345+q12345),
  -p00101 + z*(q124-q1234-q1245+q12345),
  -p10001 + z*(q234-q1234-q2345+q12345),
  -p01001 + z*(q134-q1234-q1345+q12345),
  -p00011 + z*(q123-q1234-q1235+q12345),
  -p01100 + z*(q145-q1245-q1345+q12345),
  -p00110 + z*(q125-q1235-q1245+q12345),
  -p10010 + z*(q235-q1235-q2345+q12345),
  -p10100 + z*(q245-q1245-q2345+q12345),
  -p11000 + z*(q345-q1345-q2345+q12345),
  -p11010 + z*(q35-q135-q235+q1235-q345+q1345+q2345-q12345),
  -p11001 + z*(q34-q134-q234+q1234-q345+q1345+q2345-q12345),
  -p11100 + z*(q45-q145-q245+q1245-q345+q1345+q2345-q12345),
  -p01110 + z*(q15-q125-q135+q1235-q145+q1245+q1345-q12345),
  -p10110 + z*(q25-q125-q235+q1235-q245+q1245+q2345-q12345),
  -p01101 + z*(q14-q124-q134+q1234-q145+q1245+q1345-q12345),
  -p01011 + z*(q13-q123-q134+q1234-q135+q1235+q1345-q12345),
  -p00111 + z*(q12-q123-q124+q1234-q125+q1235+q1245-q12345),
  -p10011 + z*(q23-q123-q234+q1234-q235+q1235+q2345-q12345),
  -p10101 + z*(q24-q124-q234+q1234-q245+q1245+q2345-q12345),
  -p11011 + z*(q3-q13-q23+q123-q34+q134+q234-q1234-q35+q135+q235-q1235+q345-q1345-q2345+q12345),
  -p10111 + z*(q2-q12-q23+q123-q24+q124+q234-q1234-q25+q125+q235-q1235+q245-q1245-q2345+q12345),
  -p01111 + z*(q1-q12-q13+q123-q14+q124+q134-q1234-q15+q125+q135-q1235+q145-q1245-q1345+q12345),
  -p11101 + z*(q4-q14-q24+q124-q34+q134+q234-q1234-q45+q145+q245-q1245+q345-q1345-q2345+q12345),
  -p11110 + z*(q5-q15-q25+q125-q35+q135+q235-q1235-q45+q145+q245-q1245+q345-q1345-q2345+q12345),
  -p11111 + z*(1-q1-q2+q12-q3+q13+q23-q123-q4+q14+q24-q124+q34-q134-q234+q1234-q5+q15+q25-q125+q35-q135-q235+q1235+q45-q145-q245+q1245-q345+q1345+q2345-q12345)
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
