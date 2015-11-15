#!/usr/bin/perl

$ret = 0;
$! = 1;

$f0 = "$ARGV[0].s";
$f1 = "$ARGV[1].s";

system("sort $ARGV[0] > $f0");
system("sort $ARGV[1] > $f1");

open(A, $f0) or die "can't open ", $f0;
open(B, $f1) or die "can't open ", $f1;

$prec = $ARGV[2];

while(<A>)
{

    $lineA = $_;
    $lineB = <B>;

    @tokA = split /,/, $lineA;
    @tokB = split /,/, $lineB;

    foreach $a (@tokA)
    {
        $b = shift @tokB;

        if (!compare($a, $b))
        {
            $ret = 1;
            print "< ", $lineA;
            print "> ", $lineB;
            print "\n";
            next;
        }
    }
}

while(<B>)
{
    $ret = 1;
    $lineB = $_;
    print "> ", $lineB;
    print "\n";
}

exit $ret;

##############

sub compare
{
    my $a = $_[0];
    my $b = $_[1];
    
    if ($a eq "0" || $a eq "inf" || $a == 0) # "0" or "inf" or a string that does not have a numeric conversion
    {
    	# exact string comparison
    	return $a eq $b;
    }
    else
    {
    	# numeric comparison to given precision
    	return abs($a - $b)/$a <= $prec;
    }
    
    #return $a == $b; # numeric comparison. All strings compare equal
}

