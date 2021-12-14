#!/usr/bin/env perl

use Modern::Perl 2018;
use Getopt::Long qw(:config bundling);

use CInet::Base;
use CInet::ManySAT;

use Array::Set qw(set_diff set_union set_intersect);
use Algorithm::Combinatorics qw(subsets permutations);
use List::UtilsBy qw(uniq_by);

GetOptions(
    "-M" => \my $matroids,
) or die 'failed to parse options';

my $n = shift // die 'need ground set size';
my $N = [ 1 .. $n ];
my $cube = Cube($N);

sub getvar {
    $cube->pack([ [], shift ])
}

sub getset {
    $cube->unpack(0, shift)->[1]
}

sub fmtset {
    '{' . join(', ', sort shift->@*) . '}'
}

sub str {
    '{ ' . join(', ', sort map fmtset($_), shift->@*) . ' }'
}

sub act {
    my ($FF, $p) = @_;
    [ map { [ map { $p->[$_-1] } @$_ ] } @$FF ]
}

sub facets {
    my $FF = shift;
    my @max;
    # TODO: Slow double loop over @$FF.
    F: for my $F (@$FF) {
        my $f = fmtset $F;
        for my $G (@$FF) {
            my $g = fmtset $G;
            next if $f eq $g;
            next F if set_diff($F, $G)->@* == 0;
        }
        push @max, $F;
    }
    \@max
}

# TODO:
sub connected {
    my $FF = facets shift;
    my @CC = @$FF;
    my $again = 1;
    # Join facets which intersect for as long as possible.
    while ($again) {
        $again = 0;
        my @CCnew;
        F: for my $F (@CC) {
            my $f = fmtset $F;

            # Already merged into the new array?
            for my $G (@CCnew) {
                my $g = fmtset $G;
                next if $f eq $g;
                next F if set_diff($F, $G)->@* == 0;
            }

            # Try to merge into $F.
            my @merge;
            for my $G (@CC) {
                my $g = fmtset $G;
                next if $f eq $g;
                push @merge, $G if set_intersect($F, $G)->@*;
            }

            if (@merge) {
                push @CCnew, set_union($F, @merge);
                $again++;
            }
            else {
                push @CCnew, $F;
            }
        }
        @CC = @CCnew;
    }

    @CC = uniq_by { fmtset $_ } @CC;
    @CC == 1
}

my $solver = CInet::ManySAT->new;
for my $A (subsets($N)) {
    my $a = getvar($A);
    for my $B (subsets($A)) {
        my $b = getvar($B);
        next if $a == $b;
        $solver->add([ -$a, +$b ]); # A => B
    }
}

# Only matroids requested?
if ($matroids) {
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
}

# Simplicial complex should be covering the whole set $N
$solver->add([ getvar([ $_ ]) ]) for @$N;

my @group = permutations($N);
my $all = $solver->all;
my %seen;
while (defined(my $cpx = $all->next)) {
    my $FF = [ map getset($_), grep $_ > 0, @$cpx ];
    next if $seen{str $FF};
    #next if not connected $FF;
    $seen{str act($FF, $_)}++ for @group;
    say str facets $FF;
}
