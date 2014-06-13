#!/soft/bin/perl
# DGK

use strict;
use warnings;

# Objective
# This script should parse the output of the GEM mapper
# It can parse single mappings, split mappings or paired end mappings, and it
# will select by default those mappings that are unique. In the case of paired
# ends it will choose the two closest mappings that belong to the same
# chromosome if possible
#
# To reproduce behaviour of gem2gff_multi.pl use -multi -closest no_mismatches

use Pod::Usage;
use Getopt::Long;

# Declare variables & Get command line options
my $mismatches=2;
my $closest_hit=0;
my $no_filters;
my @infiles;
my $paired=0;
my $split=0;
my $multi=0;
my $sorted=0;
my $limit=1;
my $tmpdir;
my $mismatch_strict;
my $help;
my $man;
my $read_length;
my $coord_parser;

GetOptions('mismatches|m=i' => \$mismatches,
	   'closest=i' => \$closest_hit,
	   'input|i=s' => \@infiles,
	   'paired|p' => \$paired,
	   'split|s' => \$split,
	   'multi' => \$multi,
	   'sorted' => \$sorted,
	   'mismatch_strict' => \$mismatch_strict,
	   'tmpdir=s' => \$tmpdir,
	   'limit=s' => \$limit,
	   'no_filters' => \$no_filters,
	   'help' => \$help,
	   'man' => \$man);

# Print help if needed
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

# Check we have enough information
if (@infiles) {
    if ($paired &&
	(@infiles != 2)) {
	die "Separate mapping files for both reads required for paired data\n";
    }
} else {
    pod2usage("No input supplied");
}

# Decide on the kind of input and select the subroutine to use
if ($paired && $split) {
    pod2usage("Can't handle paired end split data...Yet");
}
if ($paired && $no_filters) {
    pod2usage("ERROR: -paired and -no_filters options are mutually exclusive");
}
if ($no_filters && !$multi) {
    $multi=1;
    print STDERR "WARNING: Setting -multi option to 1 for no_filters\n";
}


# Define the parsers
my %parsers=%{get_parsers($mismatches,
			  $closest_hit,
			  $sorted)};

# Do the parsing
if ($paired) {
    print STDERR "Parsing paired reads\n";
    if ($multi) {
	print STDERR "WARNING: ignoring -multi option for paired end reads\n";
    }
    $parsers{'paired'}->($infiles[0],
			 $infiles[1]);
} elsif ($split) {
    print STDERR "Parsing split reads\n";
    # Get the splitmapping parsers
    $coord_parser=get_split_coord_parser();
    $parsers{'split'}->($infiles[0]);
} else {
    print STDERR "Parsing single reads\n";
    $parsers{'single'}->($infiles[0]);
}


exit;

# This sub should take a gem standard (not split) hit and parse it
sub process_gem_hit {
    my $match=shift;
    my $length=shift;
    my $matchsum=shift;

    if ($match =~/^-$/) {
	return();
    }

    # Get chromosome
    my ($chr,$start)=split(':',$match);

    # Get strand
    my $strand=substr($start,0,1,'');
    if ($strand eq 'F') {
	$strand='+';
    } elsif ($strand eq 'R') {
	$strand='-';
    } else {
	warn "Strange strand for $match\n";
    }

    # Get qualities
    my ($coord,$qualities);
    if ($start=~/@/) {
	($coord,$qualities)= split(/\@/,$start);
    } else {
	$coord=$start;
    }
    $coord=~/(\d+)([^\d].+)*/;
    $start=$1;
    my $mismatch_string = $2;
    my $end=$start + $length - 1;
    my $mismatch_no=0;
    if ($mismatch_string) {
	my $string=$mismatch_string;
	# Deal with insertions
	$string=~s/<.\d+>/I/g;
	$mismatch_no=$string=~s/[^\d]\d+/-/g;
    }

    unless ($qualities) {
	$qualities='None';
    }
    return([$chr,$start,$end,$strand,$mismatch_no,$qualities,$matchsum]);
}

sub parse_gem_line {
    my $line=shift;
    chomp($line);
    my %line;

    my @line=split(/\t/,$line);
    if (@line == 4) {
	# No qualities
    } elsif (@line == 5) {
	$line{'qual'}=splice(@line,2,1);
    } else {
	warn "Error in GEM format: $line\n";
    }
    $line{'id'}=$line[0];
    $line{'seq'}=$line[1];
    $line{'matches'}=$line[2];
    $line{'hits'}=$line[3];
    unless ($line{'hits'}) {
	warn "Problem parsing line: $line\n";
    }
    $line{'length'}=length($line{'seq'});

    return(\%line);
}

sub get_fh {
    my $filename=shift;
    my $write=shift;

    chomp($filename);
    my $openstring;
    my $fh;
    
    if ($write) {
	if ($filename=~/.gz$/) {
	    warn "Piping $filename through gzip\n";
	    # Use gzip -7 as the increase in compression between 7 an 9 is quite
	    # small and the extra time invested quite large
	    $openstring="| gzip -7 -c > $filename";
	} else {
	    $openstring=">".$filename;
	}
    } else {
	if ($filename=~/.gz$/) {
	    warn "Opening gzipped file $filename\n";
	    $openstring="gunzip -c $filename |";
	} else {
	    $openstring=$filename;
	}
    }
    
    open($fh,$openstring) ||
	die "Unable to open $filename: $!,$?\n";

    return($fh);
}

sub shell_sort_file {
    my $file=shift;
    if ($file=~/.gz$/) {
	die "I'm sorry, but I currently cannot sort gzipped files\n";
    }
    my $command="sort -o $file $file";
    if ($tmpdir) {
	$command.=" -T $tmpdir";
    }
    my $beforecount=`wc $file`;
    system($command);
    my $aftercount=`wc $file`;

    unless ($beforecount eq $aftercount) {
	die "Problem sorting\n";
    }
}

sub get_closest_hit {
    my $read1=shift;
    my $read2=shift;

    my $ref_coord=int(($read1->[1] + $read1->[2]) / 2);

    my $closest;
    my $dist= 1e10;
    foreach my $read (@{$read2}) {
	if ($read1->[0] ne $read->[0]) {
	    # Exclude interchromosomal
	    next;
	} elsif ($read1->[3] eq $read->[3]) {
	    # Both reads have the same strand, which should not happed in
	    # paired unless there is something strange going on
	    next;
	} elsif (($read->[3] eq '+' && $read->[2] > $read1->[1]) ||
		 ($read->[3] eq '-' && $read->[1] < $read1->[2])) {
	    # Check if the orientation is correct
	    next;
	}
	my $prob_coord=int(($read->[1] + $read->[2]) / 2);
	my $prob_dist=abs($ref_coord - $prob_coord);
	if ($prob_dist < $dist) {
	    $dist=$prob_dist;
	    $closest=$read;
	}
    }
    return($closest);
}

# This is the main body of the parsers
sub line_reader {
    my $fh=shift;

    my $reader=sub {
	if (my $line=<$fh>) {
	    my %read=%{parse_gem_line($line)};
	    
	    my @map=split(':',$read{'matches'});

	    # Get best match position
	    my $best_match=0;
	    for (my $i=0;$i<=$mismatches;$i++) {
		if ($map[$i] eq '!') {
		    $best_match=-1;
		    last;
		} elsif ($map[$i] != 0) {
		    $best_match=$i;		    
		    last;
		} 
	    }

	    # Skip hits with too many matches
	    if ($best_match < 0) {
		next;
	    }

	    # From the @map array remove those cases that have more mismatches
	    # than the mismatch threshold
	    if ($mismatch_strict) {
		splice(@map,$mismatches + 1);
	    }

	    # From the @map array get the number of matches we have to take
	    # First remove the first part of the @map array up to the first
	    # match

	    splice(@map,0,$best_match);
	    # After remove the second part from the position after the first
	    # match plus the closest hit number we want
	    my $offset=$closest_hit + 1;
	    unless ($offset > @map) {
		splice(@map,$offset);
	    }
	    my $matches=0;
	    foreach my $hits (@map) {
		$matches+=$hits;
	    }

	    if ($matches < 1) {
		# If there are no matches skip the rest
		$read{'matches'}=[];
		return(\%read);
	    }

	    my @matches=split(',',$read{'hits'});
	    if ($matches > @matches) {
		if (@matches) {
		    $matches=@matches;
		} else {
		    $read{'matches'}=[];
		    return(\%read);
		}
	    }
	    splice(@matches,$matches);
	    my @processed;
	    foreach my $match (@matches) {
		my $processed=process_gem_hit($match,
					      $read{'length'},
					      $read{'matches'});
		if ($processed) {
		    push @processed, $processed;
		}
	    }
	    $read{'matches'}=[@processed];
	    return(\%read);
	} else {
	    return();
	}
    };
    return($reader);
}

sub get_split_coord_parser {
    my %parser;

    $parser{'FF'}=sub{
	my $up_coord_string=shift;
	my $down_coord_string=shift;
	my $pos=shift;

	# Check input
	# First the up string
	my $ranges1=check_single_start($up_coord_string,
					      $pos);

	# Check the down string
	my $ranges2=check_multi_start($down_coord_string,
					     $pos,
					     1);

	return($ranges1,$ranges2);
    };
    
    $parser{'FR'}=sub{
	my $up_coord_string=shift;
	my $down_coord_string=shift;
	my $pos=shift;

	# Check input
	# First the up string
	my $ranges1=check_single_start($up_coord_string,
					      $pos);

	# Check the down string
	my $ranges2=check_single_start($down_coord_string,
					      $pos,
					      1);

	return($ranges1,$ranges2);
    };

    $parser{'RR'}=sub{
	my $up_coord_string=shift;
	my $down_coord_string=shift;
	my $pos=shift;

	# Check input
	# First the up string
	my $ranges1=check_multi_start($up_coord_string,
				      $pos);

	# Check the down string
	my $ranges2=check_single_start($down_coord_string,
					      $pos,
					      1);

	return($ranges1,$ranges2);
    };
    
    $parser{'RF'}=sub{
	my $up_coord_string=shift;
	my $down_coord_string=shift;
	my $pos=shift;

	# Check input
	# First the up string
	my $ranges1=check_multi_start($up_coord_string,
				      $pos);

	# Check the down string
	my $ranges2=check_multi_start($down_coord_string,
				      $pos,
				      1);

	return($ranges1,$ranges2);
    };

    return(\%parser);
}

sub check_multi_start {
    my $coord_string=shift;
    my $pos=shift;
    my $down=shift;
    my @ranges;

    my %starts;
    $pos=~s/(\[|])//g;
    $coord_string=~s/(\[|])//g;
    
    my @pos=sort {$a <=> $b} (split(/[-;]/,$pos));
    my ($pos_start,$pos_end)=($pos[0],$pos[-1]);
    my @splits=($pos_start .. $pos_end);
    if ($down) {
	foreach my $split (@splits) {
	    $split=$read_length - $split;
	}
	@splits=sort {$b <=> $a} @splits;
    } else {
	@splits=sort {$a <=> $b} @splits;
    }

    my @coords=split(/[-;]/,$coord_string);
    my @coords_sort=sort {$a <=> $b} @coords;
    my ($hit_start,$hit_end)=($coords_sort[0],$coords_sort[-1]);
    my @starts=($hit_start .. $hit_end);

    if ($coords[1] &&
	($coords[0] > $coords[1])) {
	# This means we are in the munus strand, so we invert the array
	@starts=reverse(@starts);
    }

    unless (@splits == @starts) {
	warn "Wrong number of starts $pos, $coord_string\n";
	print STDERR join("\t",
			  @splits),"\n";
	print STDERR join("\t",
			  @starts),"\n";
	return();
    }

    for (my $i=0;$i<@splits;$i++) {
	push @ranges, [$starts[$i],$starts[$i] + $splits[$i] - 1];
    }

    return(\@ranges);
}

sub check_single_start {
    my $coord_string=shift;
    my $pos=shift;
    my $down=shift;
    my @ranges;

    my %starts;
    $pos=~s/(\[|])//g;
    $coord_string=~s/(\[|])//g;

    my @pos=sort {$a <=> $b} (split(/[-;]/,$pos));
    my ($pos_start,$pos_end)=($pos[0],$pos[-1]);
    my @splits=($pos_start .. $pos_end);
    if ($down) {
	foreach my $split (@splits) {
	    $split= $read_length - $split;
	}
	@splits=sort {$b <=> $a} @splits;
    } else {
	@splits=sort {$a <=> $b} @splits;
    }
    my @starts=split(/[-;]/,$coord_string);
    @starts=($starts[0] .. $starts[-1]);

    unless ((@splits == @starts) ||
	    (@starts == 1)){
	warn "Wrong number of starts $pos, $coord_string\n";
	return();
    }

    for (my $i=0;$i<@splits;$i++) {
	my $start;
	if (@starts > 1) {
	    $start=$starts[$i];
	} else {
	    $start=$starts[0]
	}
	push @ranges, [$start,$start + $splits[$i] - 1];
    }

    return(\@ranges);
}

sub get_parsers {
    my $mismatches=shift;
    my $closest=shift;
    my $sorted=shift;
    my %parser;
    my $qualities;
    my %stats;

    $parser{'single'}=sub {
	my $infile=shift;
	my $outfilename=$infile;
	if ($multi) {
	    $outfilename.='.multi.gtf';
	} else {
	    $outfilename.='.unique.gtf';
	}
	$outfilename=~s/.*\///;
	$outfilename=~s/(.gem)*.\d+.map//g;

	my $infh=get_fh($infile);
	my $outfh=get_fh($outfilename,1);

	$infile=~s/.*\///;
	while (my $line=<$infh>) {
	    my %line=%{parse_gem_line($line)};

	    # Determine if we have qualities
	    if ($line{'qual'}) {
		$qualities=1;
	    }
	    my $read_length=length($line{'seq'});
	    my @map=split(':',$line{'matches'});

	    unless($no_filters) {
		if (@map < $mismatches - 1) {
		    $mismatches=@map - 1;
		    warn "File was mapped with $mismatches mismatches, reducing mismatch threshold accordingly \n";
		}

		# Get best match position
		my $best_match=0;
		for (my $i=0;$i<=$mismatches;$i++) {
		    if ($map[$i] eq '!') {
			$best_match=-1;
			last;
		    } elsif ($map[$i] != 0) {
			$best_match=$i;		    
			last;
		    } 
		}
		
		# Skip hits with too many matches
		if ($best_match < 0) {
		    $stats{'too_many'}++;
		    next;
		}
		
		# From the @map array get the number of matches we have to take
		# First remove the first part of the @map array up to the first
		# match
		splice(@map,0,$best_match);
				
		# After remove the second part from the position after the first
		# match plus the closest hit number we want
		my $offset=$closest_hit + 1;
		unless ($offset > @map) {
		    splice(@map,$offset);
		}
	    }

	    my $matches=0;
	    foreach my $hits (@map) {
		$matches+=$hits;
	    }

	    if ($matches < 1) {
		# If there are no matches skip the rest
		$stats{'no_hits'}++;
		next;
	    } elsif ($matches > 1) {
		# Unless we are going for multi matches, if the best match is
		# not unique discard it.
		$stats{'multi_hits'}++;
		unless ($multi) {
		    next;
		}
	    } else {
		$stats{'unique_hits'}++;
	    }

	    # If the script reaches this point it has found at least one good
	    # hit
	    my @matches=split(',',$line{'hits'});
	    if ($matches > @matches) {
		if (@matches) {
		    $stats{'missing_close'}++;
		    $matches=@matches;
		} else {
		    $stats{'too_many'}++;
		    next;
		}
	    }
	    splice(@matches,$matches);
	    my @processed;
	    foreach my $match (@matches) {
		my $processed=process_gem_hit($match,
					      $line{'length'},
					      $line{'matches'});
		if ($processed) {
		    push @processed, $processed;
		}
	    }

	    $line{'matches'}=[@processed];

	    foreach my $match (@{$line{'matches'}}) { 
		print_gff($outfh,
			  $infile,
			  $line{'id'},
			  $match,
			  'single');
	    }
	}
	close($infh);
	close($outfh);
	foreach my $key (keys %stats) {
	    print STDERR join("\t",
			      $key,
			      $stats{$key}),"\n";
	}
    };

    $parser{'paired'}=sub {
	my $infile1=shift;
	my $infile2=shift;
	my $type='paired';

	unless($infile1 && $infile2) {
	    die "Two files required for paired ends\n";
	}

	# Sort both files
	unless ($sorted) {
	    print STDERR 'Sorting files...';
	    shell_sort_file($infile1);
	    shell_sort_file($infile2);
	    print STDERR "done\n";
	}
 
	my $outfilename=$infile1.'.'.$infile2.'.paired.gtf';
	$outfilename=~s/gem.\d+.map//g;
	$outfilename=~s/\//./g;
	$outfilename=~s/\.+/./g;
	print STDERR "Printing to $outfilename\n";
	my %stats;
	my %distances;
	
	# Open the required files
	my $outfh=get_fh($outfilename,1);
	my $reads1fh=get_fh($infile1);
	my $reads2fh=get_fh($infile2);
	
	my $reader1=line_reader($reads1fh);
	my $reader2=line_reader($reads2fh);
	
	print STDERR 'Mixing & matching...';
	while ((my $read1=&$reader1) &&
	       (my $read2=&$reader2)){
	    $stats{'total'}++;
	    if ((@{$read1->{'matches'}}) &&
		(@{$read2->{'matches'}})) {
		# We have matches for both ends
		if ((@{$read1->{'matches'}} == 1) && 
		    (@{$read2->{'matches'}} == 1)) {
		    $stats{'double_unique'}++;
		    print_gff($outfh,
			      $infile1,
			      $read1->{'id'},
			      $read1->{'matches'}->[0],
			      $type);
		    print_gff($outfh,
			      $infile2,
			      $read2->{'id'},
			      $read2->{'matches'}->[0],
			      $type);
		} elsif (@{$read1->{'matches'}} == 1) {
		    $stats{'read1_unique'}++;
		    # select the closest hit
		    my $closest=get_closest_hit($read1->{'matches'}->[0],
						$read2->{'matches'});
		    if ($closest) {
			print_gff($outfh,
				  $infile1,
				  $read1->{'id'},
				  $read1->{'matches'}->[0],
				  $type);
			print_gff($outfh,
				  $infile2,
				  $read2->{'id'},
				  $closest,
				  $type);
		    } else {
			print_gff($outfh,
				  $infile1,
				  $read1->{'id'},
				  $read1->{'matches'}->[0],
				  'unmatched');
		    }
		} elsif (@{$read2->{'matches'}} == 1) {
		    $stats{'read2_unique'}++;
		    # select the closest hit
		    my $closest=get_closest_hit($read2->{'matches'}->[0],
						$read1->{'matches'});
		    if ($closest) {
			print_gff($outfh,
				  $infile1,
				  $read1->{'id'},
				  $closest,
				  $type);
			print_gff($outfh,
				  $infile2,
				  $read2->{'id'},
				  $read2->{'matches'}->[0],
				  $type);
		    } else {
			print_gff($outfh,
				  $infile2,
				  $read2->{'id'},
				  $read2->{'matches'}->[0],
				  'unmatched');
		    }
		} else {
		    $stats{'both_multi'}++;
		}
	    } elsif (@{$read1->{'matches'}}) {
		# We have matches only for the first end
		if (@{$read1->{'matches'}} == 1) {
		    $stats{'read1_unique'}++;
		    print_gff($outfh,
			      $infile1,
			      $read1->{'id'},
			      $read1->{'matches'}->[0],
			      'unmatched');
		} else {
		    $stats{'read1_multi'}++;
		}
		
	    } elsif (@{$read2->{'matches'}}) {
		# We have matches only for the second
		if (@{$read2->{'matches'}} == 1) {
		    $stats{'read2_unique'}++;
		    print_gff($outfh,
			      $infile2,
			      $read2->{'id'},
			      $read2->{'matches'}->[0],
			      $type);
		} else {
		    $stats{'read2_multi'}++;
		}
	    } else {
		# We have no match
		$stats{'no_hit'}++;
	    }
	}
	print STDERR "done\n";
	
	# Close the files
	close($reads1fh);
	close($reads2fh);
	close($outfh);
	
	foreach my $key (keys %stats) {
	    print STDERR join("\t",
			      $key,
			      $stats{$key}),"\n";
	}
    };
    
    $parser{'split'}=sub {
	my $infile=shift;
	my $infh=get_fh($infile);

	my %strands;
	while (my $line=<$infh>) {
	    chomp($line);
	    my @line=split(/\t+/,$line);
	    if (@line == 4) {
		# No qualities
	    } elsif (@line == 5) {
		splice(@line,2,1);
	    } else {
		warn "Error in GEM format: $line\n";
	    }

	    unless ($read_length) {
		$read_length=length($line[1]);
	    }

	    if ($line[2]=~/-/) {
		$stats{'too_many'}++;
		next;
	    }

	    my @map=split(':',$line[2]);

	    unless ($no_filters) {
		if (@map < $mismatches - 1) {
		    $mismatches=@map - 1;
		    warn "File was mapped with $mismatches mismatches, reducing mismatch threshold accordingly \n";
		}

		# Get best match position
		my $best_match=0;
		for (my $i=0;$i<=$mismatches;$i++) {
		    if (($map[$i] eq '!') ||
			($map[$i] eq '-')) {
			$best_match=-1;
			last;
		    } elsif ($map[$i] != 0) {
			$best_match=$i;		    
			last;
		    } 
		}
		
		# Skip hits with too many matches
		if ($best_match < 0) {
		    $stats{'too_many'}++;
		    next;
		}

		# From the @map array get the number of matches we have to take
		# First remove the first part of the @map array up to the first
		# match
		splice(@map,0,$best_match);
		
		# After remove the second part from the position after the first
		# match plus the closest hit number we want
		my $offset=$closest_hit + 1;
		unless ($offset > @map) {
		    splice(@map,$offset);
		}
	    }

	    my $matches=0;

	    foreach my $hits (@map) {
		$matches+=$hits;
	    }

	    if ($matches < 1) {
		# If there are no matches skip the rest
		$stats{'no_hits'}++;
		next;
	    } elsif ($matches > 1) {
		# Unless we are going for multi matches, if the best match is
		# not unique discard it.
		$stats{'multi_hits'}++;
		unless ($multi) {
		    next;
		}
	    } else {
		$stats{'unique_hits'}++;
	    }

	    # If the script reaches this point it has found at least one good
	    # hit
	    my @matches=split(',',$line[3]);
	    if ($matches > @matches) {
		if (@matches) {
		    $stats{'missing_close'}++;
		    $matches=@matches;
		} else {
		    $stats{'too_many'}++;
		    next;
		}
	    }
	    splice(@matches,$matches);

	    # Parse the hits in this case we are taking only the best one
	    foreach my $match (@matches) {
		my ($pos,$coord_string)=split('=',$match);
		my %hit=%{get_coords($coord_string,
				     $pos)};

		$strands{$hit{'up_strand'}.$hit{'down_strand'}}++;
	    
		# Print two Gtf entries per split read using the same id to be
		# able to group them after
		my $read_id=join('; ',
				 'ID "'.$line[0].'"',
				 'Identity "'.$line[2].'"',
				 'Split "'.$pos.'"');

		# set the strands to + or -
		my $up_strand=$hit{'up_strand'};
		if ($up_strand eq 'F') {
		    $hit{'up_strand'}='+';
		} elsif ($up_strand eq 'R') {
		    $hit{'up_strand'}='-';
		} else {
		    warn "Unknown strand $read_id\n";
		}
		my $down_strand=$hit{'down_strand'};
		if ($down_strand eq 'F') {
		    $hit{'down_strand'}='+';
		} elsif ($down_strand eq 'R') {
		    $hit{'down_strand'}='-';
		} else {
		    warn "Unknown strand $read_id\n";
		}
		
		$infile=~s/.*\///;
		for (my $i=0;$i<@{$hit{'up_start'}};$i++) {
		    unless ($hit{'up_start'}->[$i] &&
			    $hit{'down_end'}->[$i]) {
			warn "Problem parsing line: $i $line\n";
			print STDERR join("\t",@{$hit{'up_start'}->[$i]}),"\n";
			next;
		    }
		    print join("\t",
			       $hit{'up_chr'},
			       $infile,
			       'splitread',
			       @{$hit{'up_start'}->[$i]},
			       '.',
			       $hit{'up_strand'},
			       '.',
			       $read_id),"\n";
		    print join("\t",
			       $hit{'down_chr'},
			       $infile,
			       'splitread',
			       @{$hit{'down_end'}->[$i]},
			       '.',
			       $hit{'down_strand'},
			       '.',
			       $read_id),"\n";
		}
	    }
	}
	    foreach my $key (keys %stats) {
		print STDERR join("\t",
				  $key,
				  $stats{$key}),"\n";
	    }
    };
    
    return(\%parser);
}

sub print_gff {
    my $outfh=shift;
    my $readfile=shift;
    my $read_id=shift;
    my $read=shift;
    my $type=shift;

    my $string= join(' ',
		     'read_id','"'.$read_id.'";',
		     'mismatches','"'.$read->[4].'";',
		     'qualities','"'.$read->[5].'";',
		     'matches','"'.$read->[6].'";');
    print $outfh join("\t",
		      $read->[0],
		      $readfile,
		      $type.'_read',
		      $read->[1],
		      $read->[2],
		      '.',
		      $read->[3],
		      '.',
		      $string),"\n";
}

sub get_coords {
    my $string=shift;
    my $pos=shift;
    my %hit;

    # First we need to parse everything
    # Get the strand of the upstream and of the downstream
    my ($up_string,$down_string)=split(/~/,$string);
    my ($up_coords,$down_coords);
    ($hit{'up_chr'},$up_coords)=split(':',$up_string);
    ($hit{'down_chr'},$down_coords)=split(':',$down_string);
    $hit{'up_strand'}=substr($up_coords,0,1,'');
    $hit{'down_strand'}=substr($down_coords,0,1,'');

    # Depending of the way the strands are combined we must parse in different
    # ways
    my $parser_key=$hit{'up_strand'}.$hit{'down_strand'};
    if (exists $coord_parser->{$parser_key}){
	($hit{'up_start'},
	 $hit{'down_end'})=$coord_parser->{$parser_key}->($up_coords,
							 $down_coords,
							 $pos);
    } else {
	warn "Unknown strand combination $parser_key\n";
    }

    return(\%hit);
}


__END__

=head1 NAME
    
    ?
    
=head1 SYNOPSIS
    
    gem2gff_parser.pl -i <file/s>
    
  Options:
    -help            brief help message
    -man             full documentation

    Mandatory
    -input/i         Indicates the input file or files (2 for paired end reads)

    Optional:
    -mismatches      Maximum allowed number of mismatches. Default 2
    -closest         Number of mismatches more than the best hit allowed for
                     keeping a hit. Default 0
    -paired|p        Input are two files with paired reads. Default Not paired
    -multi           Keep all hits within mismatch threshold. Default keep only
                     unique hits.This option is ignored when using paired ends.
                     Check -man for more info.
    -sorted          Input files are already sorted. Default sort input files.
                     Only used for paired input
    -tmpdir          Temporary dir to use for file sorting. Defaults to the
                     default dir for the system sort
    -mismatch_strict Only meaningful if -closest > 0. Check -man for more info.
    -split           Input is a file resulting from splitmapping.
    -no_filters      All mappings will be returned

=head1 OPTIONS
    
=over 8
    
=item B<-help>
    
    Print a brief help message and exits.
    
=item B<-man>
    
    Prints the manual page and exits.

=item B<-multi>

    Keep all hits within mismatch threshold, not only unique.
    This option is ignored when using paired ends because as far as this script
    is concerned, the point of using them is to recover unique maps from
    multimaps.

=item B<-mismatch_strict>

    Only meaningful if -closest > 0. In this case those matches within
    [closest] mismatches from the best hit will be kept by default even if they 
    have more than [mismatches] mismatches. The reason for this is that a hit
    with one or two more mismatches may be equally valid if htere is a SNP, so
    a unique hit with no hit within a threshold of 2 mismatches is more unique
    than one for which there is another hit with one more mismatch. This flag
    will discard nay match with more mismatches than [mismatch] regardless of it
    being within the [closest] threshold

=item B<-no_filters>

    All mappings will be returned. Ignores the options -mismatchces and closest
    and automatically sets multi. Does not work with paired.

=back
    
=head1 DESCRIPTION
    
    This program is not documented yet, sorry

=cut
