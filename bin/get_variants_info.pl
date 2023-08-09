#!/usr/bin/env perl
#
# Shiping Zhang
#

BEGIN {

  use Getopt::Std;
  use DBI;
  use FileHandle;

  STDOUT->autoflush();
}

END {
}

use strict;

use vars qw(@args $opt_o);
use vars qw($file $tmpin $tmpout);
use vars qw(%versions %vars %variants @variants $variant @varset);
use vars qw($pos $ref $alt $status $locus $type %dz $dz);
use vars qw($from1 $from2 $to1 $to2 $range $cr);
use vars qw($total $count $freq $row %mitomap %gnomad);
use vars qw($yes %yes);
use vars qw($cmd $jar $name $qual $haplo $haplov);
use vars qw($ntchg $aachg $conserv %conservs %mitotips);
use vars qw(%loci $start $end $strand);
use vars qw($num $aa1 $aa2 %mutations);
use vars qw($db $dbh $sql $sth);
use vars qw(@nts $path $msge);

$path = $0;
$path =~ s'/[^/]*$'';

$0 =~ s/.*\///;

@args = @ARGV;

getopts('o:') or usage();

usage() if @ARGV < 2;

$tmpin = "/tmp/$$.hsd";
$tmpout = "/tmp/$$.out";

# $jar = '/scr1/users/liuc9/tools/haplogrep3/haplogrep3.jar';
$jar = shift;

if(! -s $jar)
{
    die "$jar doesn't exist!\n";
}

$cmd = "java -jar $jar classify --extend-report --in $tmpin --out $tmpout --tree phylotree-rcrs\@17.2 >& /dev/null";

$file = "/tmp/$0.log";
open(LOG, ">$file") or die "Can't write $file\n";

$db = shift;
$file = shift;

print LOG $cmd, "\n";

$cr = 576;

$dbh = DBI->connect("dbi:SQLite:dbname=$db","","") or die "Cannot connect: $DBI::errstr"; ;

$sql = 'select sequence from refseq';
$sth = $dbh->prepare($sql) or die $dbh->errstr;
$sth->execute() or die $dbh->errstr;
$row = $sth->fetchrow_arrayref;
$sth->finish();
@nts = split('', ' ' . $row->[0]);

$sql = 'select pos, ref, alt, prediction_perc from mitotip';
$sth = $dbh->prepare($sql) or die $dbh->errstr;
$sth->execute() or die $dbh->errstr;
while($row = $sth->fetchrow_arrayref)
{
    ($pos, $ref, $alt, $ntchg) = @$row;
    $mitotips{"$ref$pos$alt"} = sprintf("%.02f%%", $ntchg * 100);
}
$sth->finish();

$sql = 'select pos, rate from conservation_rate';
$sth = $dbh->prepare($sql) or die $dbh->errstr;
$sth->execute() or die $dbh->errstr;
while($row = $sth->fetchrow_arrayref)
{
    ($pos, $conserv) = @$row;
    $conservs{$pos} = $conserv;
}
$sth->finish();

#$sql = 'select locus, strand, pos, ref, alt, num, codon1, codon2, aa1, aa2 from aachanges';
$sql = 'select locus, pos, ref, alt, num, aa1, aa2 from aachanges';
$sth = $dbh->prepare($sql) or die $dbh->errstr;
$sth->execute() or die $dbh->errstr;
while($row = $sth->fetchrow_arrayref)
{
    ($locus, $pos, $ref, $alt, $num, $aa1, $aa2) = @$row;
    push @{$mutations{"$ref$pos$alt"}}, [$locus, $aa1, $num, $aa2];
}
$sth->finish();

$sql = "select count from summary where name='genbank'";
$sth = $dbh->prepare($sql) or die $dbh->errstr;
$sth->execute();
$row = $sth->fetchrow_arrayref;
$total = $row->[0];
$sth->finish();

$sql = 'select distinct position, refna, regna, dz, status from (select position, refna, regna, dz, status from mmutation union select position, refna, regna, dz, status from rtmutation) t';
$sth = $dbh->prepare($sql) or die $dbh->errstr;
$sth->execute();
while($row = $sth->fetchrow_arrayref)
{
    ($pos, $ref, $alt, $dz, $status) = @$row;
    $dz{"$ref$pos$alt"} = [$dz, $status];
}
$sth->finish();

$sql = 'select tpos, tnt, qnt, fl_count from variants_count';
$sth = $dbh->prepare($sql) or die $dbh->errstr;
$sth->execute() or die $dbh->errstr;
while($row = $sth->fetchrow_arrayref)
{
    ($pos, $ref, $alt, $count) = @$row;
    $mitomap{"$ref$pos$alt"} = [$count, sprintf("%.03f", $count*100/$total)];
}
$sth->finish();

$sql = 'select position, refna, regna, ac_hom, af_hom from gnomad';
$sth = $dbh->prepare($sql) or die $dbh->errstr;
$sth->execute();
while($row = $sth->fetchrow_arrayref)
{
    ($pos, $ref, $alt, $count, $freq) = @$row;
    $gnomad{"$ref$pos$alt"} = [$count, sprintf("%.03f", $freq*100)];
}
$sth->finish();

$sql = "select common_name, starting, ending, strand from locus where type='m'";
$sth = $dbh->prepare($sql) or die $dbh->errstr;
$sth->execute() or die $dbh->errstr;
while($row = $sth->fetchrow_arrayref)
{
    ($locus, $start, $end, $strand) = @$row;
    $loci{$locus} = [$start, $end, $strand];
}
$sth->finish();

$sql = "select common_name, type from locus where starting < ending and starting <= ? and ending >= ? or starting > ending and (starting <= ? or ending >= ?) order by CASE type WHEN 't' THEN 3 WHEN 'm' THEN 2 WHEN 'r' THEN 1 ELSE 0 END, upper(common_name)";
$sth = $dbh->prepare($sql) or die $dbh->errstr;

if($opt_o)
{
    open(OUT, ">$opt_o") or die "Can't write $opt_o\n";
}
else
{
    *OUT = *STDOUT;
}

print OUT <<END;
Sample\tPosition\tRef\tAlt\tLocus\tHaplogroup\tVerbose_haplogroup\tDisease\tStatus\tntchange\taachange\tConservation\tMitomap Frequency\tGnomad Frequency
END

open(IN, $file) or die "Can't open $file\n";
while(<IN>)
{
    #print;
    next if /^\s*$/;
    next if /^\s*#/;
    @args = ();
    s/#.*$//;
    s/^\s+//;
    s/\s+$//;
    chomp;
    foreach $variant (split(/\s+/))
    {
        push @args, $variant if $variant;
    }
    $name = shift @args;
    get_info();
}
close(IN);
close(LOG);

if($opt_o)
{
    close(OUT);
}

unlink($tmpin);
unlink($tmpout);

exit(0);

sub get_info
{
    @variants = ();
    @varset = ();
    %variants = ();

    # start with a real position

    $from1 = $cr;
    $to1 = -1;
    $from2 = 16569;
    $to2 = -1;

    foreach (@args)
    {
        $variant = $_;
        $variant =~ y/a-z/A-Z/;
        if($variant =~ /^(\d+)\.1([A-Z]+)$/)  # 709.1G
        {
            $pos = $1;
            $ref = $nts[$pos];
            $alt = $2;
        }
        elsif($variant =~ /^(\d+)([A-Z]+)>([A-Z]+)$/)  # 146T>C
        {
            $pos = $1;
            $ref = $2;
            $alt = $3;
            if($ref ne $nts[$pos])
            {
                $msge = "$ref != $nts[$pos] at position $pos for sample $name and variant $variant\n";
                print LOG $msge;
                die $msge;
            }
            $variant =~ s/[A-Z]+>//;
        }
        elsif($variant =~ /^([A-Z]+)(\d+)([A-Z]+)$/)  # T146C
        {
            $pos = $2;
            $ref = $1;
            $alt = $3;
            if($ref ne $nts[$pos])
            {
                $msge = "$ref != $nts[$pos] at position $pos for sample $name and variant $variant\n";
                print LOG $msge;
                die $msge;
            }
            $variant =~ s/^[A-Z]+//;
        }
        elsif($variant =~ /^(\d+)([A-Zd:]+)$/)   # 146C, 263d
        {
            $pos = $1;
            $ref = $nts[$pos];
            $alt = $2;
        }
        else
        {
            $msge = "$variant is not a valid variant.\n";
            print LOG $msge;
            die $msge;
        }

        if($alt =~ /[A-Z][A-Z]/)  # insertion
        {
            if($alt =~ /\d$ref/)
            {
                $variant =~ s/$ref/.1/;
            }
            else
            {
                $variant = '';
            }
        }
        push @variants, [$pos, $ref, $alt];
        $vars{"$ref$pos$alt"} = $variant;
        if($variant)
        {
            $variants{$variant} = $_;
            if($pos > $cr)
            {
                if($from2 > $pos)
                {
                    $from2 = $pos;
                }
                if($to2 < $pos)
                {
                    $to2 = $pos;
                }
            }
            else
            {
                if($from1 > $pos)
                {
                    $from1 = $pos;
                }
                if($to1 < $pos)
                {
                    $to1 = $pos;
                }
            }
            push @varset, $variant;
        }
    }

    %yes = ();
    $haplo = '';
    $qual = '';
    if(@varset)
    {
        $pos = join(' ', @varset);
        $range = '';
        if($to1 > 0)
        {
           $range = "$from1-$to1;";
        }
        if($to2 > 0)
        {
           $range .= "$from2-$to2;";
        }


        open(TMP, ">$tmpin") or die "Can't write $tmpin!\n";
        print TMP "ID\tRange\tHaplogroup\tPolymorphisms\n";
        print TMP "Sample\t$range\t?\t$pos\n";
        close(TMP);

        if(system($cmd))
        {
            die "$cmd failed!\n$!\n$?\n";
        }

        open(RST, $tmpout);
        $_ = <RST>;
        # print $_;
        $_ = <RST>;
        # print $_;
        close(RST);

        if($_)
        {
            (undef, $haplov, undef, $qual, undef, undef, $yes) = split(/\t/);
            $haplov =~ s/"//g;
            $qual =~ s/"//g;
            if($haplov =~ /^([A-Za-z\d+])[^A-Za-z\d]/)
            {
                $haplo = $1;
            }
            elsif($haplov =~ /^(...)./)
            {
                $haplo = $1;
            }
            else
            {
                $haplo = $haplov;
            }
            $yes =~ s/"//g;
            foreach (split(/\s+/, $yes))
            {
                $yes{$_} = 1;
            }
        }
    }

    $conserv = '';
    foreach (sort bypos @variants)
    {
        ($pos, $ref, $alt) = @$_;
        $variant = "$ref$pos$alt";
        $sth->execute($pos, $pos, $pos, $pos) or die $dbh->errstr;
        $locus = '';
        $aachg = '';
        while($row = $sth->fetchrow_arrayref)
        {
            #if($row->[1] =~ /CR:/)
            #{
            #    $row->[0] = $row->[1];
            #}
            #elsif($row->[2] eq 't')
            #{
            #    $row->[0] =~ s/^MT-T/tRNA./;
            #}
            #else
            #{
            #    $row->[0] =~ s/^MT-//;
            #}
            $locus .= ',' if $locus;
            $locus .= $row->[0];
            $type = $row->[1];

            if($type eq 'm')
            {
                $start = length($ref);
                $end = length($alt);
                if($start > 1 || $end > 1)
                {
                    if($start > $end)
                    {
                        if(($start - $end)%3)
                        {
                            $aachg = 'frameshift';
                        }
                    }
                    elsif($start < $end)
                    {
                        if(($end - $start)%3)
                        {
                            $aachg = 'frameshift';
                        }
                    }
                }
                else
                {
                    foreach (@{$mutations{"$ref$pos$alt"}})
                    {
                        $aachg .= ',' if $aachg;
                        $aachg .= "$_->[0]:$_->[1]$_->[2]$_->[3]";
                    }
                }
            }
        }
#Sample\tPosition\tRef\tAlt\tLocus\tHaplogroup\tVerbose_haplogroup\tDisease\tStatus\tntchange\taachange\tConservation\tMitomap Frequency\tGnomad Frequency
        $sth->finish();
        print OUT "$name\t$pos\t$ref\t$alt\t$locus";
        if($yes{$vars{$variant}})
        {
            print OUT "\t$haplo\t$haplov";
        }
        else
        {
            print OUT "\t\t";
        }
        if($dz{$variant})
        {
            ($dz, $status) = @{$dz{$variant}};
        }
        else
        {
            $dz = $status = '';
        }
        print OUT "\t", $dz, "\t", $status;
        if($alt eq 'd' || $alt eq ':')
        {
            $ntchg = 'deletion';
        }
        elsif(length($alt) > length($ref))
        {
            $ntchg = 'insertion';
        }
        elsif($ref =~ /[AG]/ && $alt =~ /[AG]/ ||
              $ref =~ /[TC]/ && $alt =~ /[TC]/)
        {
            $ntchg = 'transition';
        }
        else
        {
            $ntchg = 'transversion';
        }
        print OUT "\t", $ntchg;
        if($type eq 'n')
        {
            $aachg = 'non-coding';
        }
        elsif($type eq 'r')
        {
            $aachg = 'rRNA';
        }
        elsif($type eq 't')
        {
            $aachg = 'MitoTIP ' . $mitotips{"$ref$pos$alt"};
        }
        print OUT "\t", $aachg;
        $conserv = $conservs{$pos};
        print OUT "\t", $conserv;
        if($mitomap{$variant})
        {
            ($count, $freq) = @{$mitomap{$variant}};
        }
        else
        {
            $freq = '';
        }
        print OUT "\t", $freq;
        if($gnomad{$variant})
        {
            ($count, $freq) = @{$gnomad{$variant}};
        }
        else
        {
            $freq = '';
        }
        print OUT "\t", $freq;
        print OUT "\n";
    }
}

sub aachange
{
    my ($off, $aa, $cdn);

    $aa = shift;
    ($start, $end, $strand) = @{$loci{$aa}};
    if($strand eq 'H')
    {
        $off = $pos - $start;
        $off -= $off%3;
        $off += $start;
        $cdn = join('', @nts[$off .. $off+2]);
    }
    else
    {
    }
}

sub bypos
{
    if($a->[0] != $b->[0])
    {
        $a->[0] <=> $b->[0];
    }
    elsif($a->[1] ne $b->[1])
    {
        $a->[1] cmp $b->[1];
    }
    else
    {
        $a->[2] cmp $b->[2];
    }
}

sub usage
{
  my $usage = shift;
  $usage = "\n" . $usage . "\n" if $usage;
  $usage .= <<END;

  $0 - get variants info

  Usage: $0 [options] <db_file> <variants_file>

    -o <file>: specify output file (default screen)

  Error messages are logged in /tmp/$0.log
END
    $usage .= "\n";

  die $usage;

}

