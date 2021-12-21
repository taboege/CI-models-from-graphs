#!/usr/bin/env perl

use Modern::Perl 2018;
use autodie;
use utf8;

use CInet::Base;
use CInet::ManySAT::Incremental;

use Path::Tiny;
use List::Util qw(uniq uniqstr any max sum0);
use Array::Set qw(set_diff set_union set_intersect);
use Algorithm::Combinatorics qw(subsets permutations);

use Text::Table;

my $prefix = shift // die 'need prefix';
my ($n) = $prefix =~ /-(\d+)$/;
my $N = [ 1 .. $n ];
my $cube = Cube($N);

sub getvar {
    $cube->pack([ [], shift ])
}

sub getset {
    $cube->unpack(0, shift)->[1]
}

sub fmtset {
    join('', sort shift->@*)
}

sub str {
    '[' . join(',', sort map fmtset($_), shift->@*) . ']'
}

sub act {
    my ($FF, $p) = @_;
    [ map { [ map { $p->[$_-1] } @$_ ] } @$FF ]
}

# SAT solver for checking the matroid property of a complex.
my $matroids = do {
    my $solver = CInet::ManySAT::Incremental->new;
    for my $A (subsets($N)) {
        my $a = getvar($A);
        for my $B (subsets($A)) {
            my $b = getvar($B);
            next if $a == $b;
            $solver->add([ -$a, +$b ]); # A => B
        }
    }
    # Simplicial complex should be covering the whole set $N
    $solver->add([ getvar([ $_ ]) ]) for @$N;

    # Steinitz property
    for my $A (subsets($N)) {
        my $a = getvar($A);
        for my $B (subsets($N)) {
            next unless $B->@* < $A->@*;
            my $b = getvar($B);
            my $D = set_diff($A, $B);
            my @c = map { getvar(set_union($B, [$_])) } @$D;
            # A and B => B \cup c for some c in A\B
            $solver->add([ -$a, -$b, @c ]);
        }
    }

    $solver
};

# Lookup table for the forbidden minors of graphic matroids.
# FIXME: For simplicity, we only include ONE of the FIVE forbidden
# minors because we tabulate only complexes on n <= 6.
my %nongraphs = do {
    # Forbidden minors:
    #   U_{2,4}    -> rank 2 on n=4
    #   F_7        -> rank 3 on n=7
    #   F_7^*      -> rank 4 on n=7
    #   K_5^*      -> rank 6 on n=10
    #   K_{3,3}^*  -> rank 4 on n=9
    my @minors = (
        [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]],
    );
    warn "Graphicness is only implemented for n <= 6" unless $n <= 6;

    my %nongraphs;
    for my $facets (@minors) {
        my $A = [ uniq sort map @$_, @$facets ];
        my @group = permutations($A);
        $nongraphs{0+ @$A} = [ uniqstr sort map str(act($facets, $_)), @group ];
    }
    %nongraphs
};

# Saturate a list of sets (facets) to the generated simplicial complex.
sub faces {
    my $facets = shift;
    my @N = uniq sort map { @$_ } @$facets;
    [ grep {
            my $A = $_;
            any { set_diff($A, $_)->@* == 0 } @$facets
        } subsets(\@N)
    ]
}

sub facets {
    my $faces = shift;
    my @max;
    F: for my $F (@$faces) {
        my $f = fmtset $F;
        for my $G (@$faces) {
            my $g = fmtset $G;
            next if $f eq $g;
            next F if set_diff($F, $G)->@* == 0;
        }
        push @max, $F;
    }
    \@max
}

# Enumerate all $k-minors of a given simplicial complex (the definition
# is that of matroids in terms of bases, i.e. facets).
sub minors {
    my ($facets, $k) = @_;

    sub reindex {
        my $facets = shift;
        my $A = [ uniq sort map @$_, @$facets ];
        my %lut = do {
            my $i = 1;
            map { $_ => $i++ } @$A
        };

        [ map { [ @lut{@$_} ] } @$facets ]
    }

    sub dual {
        my $facets = shift;
        my $A = [ uniq sort map @$_, @$facets ];
        [ map set_diff($A, $_), @$facets ]
    }

    sub restr {
        my ($facets, $T) = @_;
        facets [ map set_diff($_, $T), @$facets ]
    }

    sub contr {
        my ($facets, $T) = @_;
        dual restr(dual($facets), $T)
    }

    (
        map(reindex(restr($facets, $_)), subsets($N, $k)),
        map(reindex(contr($facets, $_)), subsets($N, $n-$k)),
    )
}

sub parse_complex {
    faces eval(shift =~ tr/{}/[]/r);
}

sub read_complexes {
    map parse_complex($_), path(shift)->lines_utf8
}

sub parse_betti {
    my %betti;
    my $str = '???';

    my $path = path(shift);
    return %betti unless $path->is_file;

    for ($path->lines_utf8) {
        if (m|^(\{.*\})$|) {
            $str = str facets parse_complex $1;
        }
        elsif (m|^dim/codim/deg: (\d+)/(\d+)/(\d+)|) {
            my ($codim, $degree) = ($2, $3);
            $betti{$str}->{pdim} = 2 ** $n - 1 - $codim;
            $betti{$str}->{codim} = $codim;
            $betti{$str}->{degree} = $degree;
        }
        elsif (m|^\s*(\d+): (.*)|) {
            my ($degree, $line) = ($1, $2);
            next if $degree == 0;
            push $betti{$str}->{mvector}->@*, sum0($line =~ /\d+/g);
        }
    }

    %betti
}

sub parse_ed {
    my %EDdeg;
    my %maxcount;

    my $path = path(shift);
    return %EDdeg unless $path->is_file;

    my $do_update = sub {
        my ($str, $key, $new) = @_;
        my $old = $EDdeg{$str}->{$key} // 0;
        my $max = max($old, $new);
        if ($old == $max) { $maxcount{$str}->{$key}++ }
        else { $maxcount{$str}->{$key} = 1 }
        $maxcount{$str}->{total}++;
        $EDdeg{$str}->{$key} = $max;
    };

    for ($path->lines_utf8) {
        my ($cpx, $aff, $proj) = /(\{.*\}): affine\((\d+)\), projective\((\d+)\)/;
        my $str = str facets parse_complex $cpx;
        $do_update->($str, affine => $aff);
        $do_update->($str, projective => $proj);
    }

    # Signalize uncertainty (maximum is attained less than half the time)
    for my $str (keys %EDdeg) {
        for my $key (qw(affine projective)) {
            my ($found, $total) = $maxcount{$str}->@{$key, 'total'};
            $EDdeg{$str}->{$key} .= '?' if 2*$found < $total / 2; # $total / 2 because of affine+projective
        }
    }

    %EDdeg
}

sub parse_ml {
    my %MLdeg;
    my %maxcount;

    my $path = path(shift);
    return %MLdeg unless $path->is_file;

    my $do_update = sub {
        my ($str, $key, $new) = @_;
        my $old = $MLdeg{$str}->{$key} // 0;
        my $max = max($old, $new);
        if ($old == $max) { $maxcount{$str}->{$key}++ }
        else { $maxcount{$str}->{$key} = 1 }
        $maxcount{$str}->{total}++;
        $MLdeg{$str}->{$key} = $max;
    };

    for ($path->lines_utf8) {
        my ($cpx, $ml) = /(\{.*\}): (\d+)/;
        my $str = str facets parse_complex $cpx;
        $do_update->($str, ml => $ml);
    }

    # Signalize uncertainty (maximum is attained less than half the time)
    for my $str (keys %MLdeg) {
        my ($found, $total) = $maxcount{$str}->@{'ml', 'total'};
        $MLdeg{$str}->{ml} .= '?' if 2*$found < $total;
    }

    %MLdeg
}

sub fvector {
    my $faces = shift;
    my @f;
    for my $i (0 .. $n) {
        my $k = grep { $_->@* == $i } @$faces;
        last if $k == 0;
        push @f, $k;
    }
    \@f
}

sub is_pure {
    my $facets = shift;
    my $dim = $facets->[0]->@*;
    for my $f (@$facets) {
        return 0 if $dim != $f->@*;
    }
    return 1;
}

sub is_matroid {
    my $faces = shift;
    my %lut = map { fmtset($_) => -1 } subsets($N);
    $lut{fmtset($_)} = +1 for @$faces;
    my @assump = map { $lut{fmtset($_)} * getvar($_) } subsets($N);
    defined $matroids->model(\@assump)
}

sub is_graph {
    my $faces = shift;
    return 0 unless is_matroid($faces);
    my $facets = facets $faces;

    for my $k (keys %nongraphs) {
        next unless $k <= $n;
        return 0 if grep {
                my $m = str $_;
                any { $m eq $_ }
                    $nongraphs{$k}->@*;
            } minors($facets, $k);
    }
    return 1;
}

my @complexes = read_complexes "$prefix.txt";
my %betti = parse_betti "$prefix-betti-mingens.txt";
my %EDdeg = parse_ed "$prefix-EDdeg-10x.txt";
my %MLdeg = parse_ml "$prefix-MLdeg-10x.txt";

sub YN {
    map { $_ ? 'Y' : 'N' } @_
}

sub NA {
    map { $_ // 'N/A' } @_
}

sub PAR {
    '(' . join(',', (shift // [])->@*) . ')'
}

sub make_columns {
    my @new;
    my $sep = {
        is_sep => 1,
        title  => ' ' x 4,
        body   => ' ' x 4,
    };
    for my $i (0 .. $#_) {
        push @new, {
            title => $_[$i],
            align => 'left',
        };
        push @new, $sep
            unless $i == $#_;
    }
    @new
}

my $tt = Text::Table->new(
    make_columns
    'complex', 'cdim', 'f-vector', 'faces', 'facets', 'pure', 'matroid', 'graph',
    'pdim', 'codim', 'degree', 'm-vector',
    'aED', 'pED', 'ML'
);

for my $cpx (@complexes) {
    my $faces = faces $cpx;
    my $facets = facets $cpx;
    my $str = str $facets;
    $tt->add(
        $str,
        max(map { $_->$#* } @$facets),
        PAR(fvector($faces)),
        0+ $faces->@*,
        0+ $facets->@*,
        YN(is_pure($facets)),
        YN(is_matroid($faces)),
        YN(is_graph($faces)),
        NA($betti{$str}->@{'pdim', 'codim', 'degree'}),
        PAR($betti{$str}->{'mvector'}),
        NA($EDdeg{$str}->@{'affine', 'projective'}),
        NA($MLdeg{$str}->{'ml'}),
    );
}

say $tt =~ s/\s+$//gmr;
