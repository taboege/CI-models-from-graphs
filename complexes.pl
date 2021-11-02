#!/usr/bin/env perl

use Modern::Perl 2018;

use CInet::Base;
use CInet::ManySAT;

use Array::Set qw(set_diff set_union set_intersect);
use Algorithm::Combinatorics qw(subsets permutations);
use List::UtilsBy qw(uniq_by);

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

sub reduce {
    my @group = permutations($N);
    my (@rep, %seen);
    for my $FF (@_) {
        next if $seen{str $FF};
        push @rep, $FF;
        $seen{str act($FF, $_)}++ for @group;
    }
    @rep
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
# Simplicial complex should be covering the whole set $N
$solver->add([ getvar([ $_ ]) ]) for @$N;

my @cpxes;
my $all = $solver->all;
while (defined(my $cpx = $all->next)) {
    push @cpxes, [ map getset($_), grep $_ > 0, @$cpx ];
}

# Reduce modulo symmetry.
@cpxes = reduce @cpxes;
# Filter for connected complexes.
@cpxes = grep { connected $_ } @cpxes;
# Output only their facets.
say str facets $_ for @cpxes;
