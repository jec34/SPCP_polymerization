#!/usr/bin/perl

###############################################################################
#
# polym_init.pl
# This file is part of the Polymatic distribution.
#
# Description: Performs a polymerization initialization for use with the
# Polymatic code. Finds all linker atoms as defined by the 'link' command in
# the input script and adds the appropriate artificial charges as given by the
# 'charge' command in the input script. Reads in and writes out a LAMMPS data
# file.
#
# Author: Lauren J. Abbott
# Version: 1.0
# Date: February 15, 2013
#
# Syntax:
#  ./polym.pl -i data.lmps
#             -t types.txt
#             -s polym.in
#             -o new.lmps
#
# Parameters:
#  1. data.lmps, LAMMPS data file of initial system (-i)
#  2. types.txt, data types text file (-t)
#  3. polym.in, input script specifying polymerization options (-l)
#  4. new.lmps, updated LAMMPS data file after polymerization step (-o)
#
###############################################################################
#
# Polymatic: a general simulated polymerization algorithm
# Copyright (C) 2013 Lauren J. Abbott
#
# Polymatic is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# Polymatic is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License (COPYING)
# along with Polymatic. If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

use strict;
use Math::Trig();
use POSIX();

# Variables
my ($fileData, $fileTypes, $fileInput, $fileOut);
my ($header, $lengthA, $lengthB, $lengthC, $xLo, $xHi, $yLo, $yHi, $zLo, $zHi);

my ($numAtoms, $numBonds, $numAngles, $numDiheds, $numImprops, $numMols);
my (@atomMol, @atomType, @atomQ, @atomPos); # atomID => value
my (@bonds, @angles, @diheds, @improps, @molecules); # bondID => atoms
my (@atomBonds, @atomAngles, @atomDiheds, @atomImprops); # atomID => atoms
my (@atomBondNums, @atomAngleNums, @atomDihedNums, @atomImpropNums);

my ($numAtomTypes, $numBondTypes, $numAngleTypes);
my ($numDihedTypes, $numImpropTypes);
my (@atomTypes, @bondTypes, @angleTypes); # id => string
my (@dihedTypes, @impropTypes);
my (%atomTypesKey, %bondTypesKey, %angleTypesKey); # string => id
my (%dihedTypesKey, %impropTypesKey);

my (@masses, @pairCoeffs, @bondCoeffs);
my (@angleCoeffs, @bondBondCoeffs, @bondAngCoeffs);
my (@dihedCoeffs, @midBondTorsCoeffs, @endBondTorsCoeffs);
my (@angTorsCoeffs, @angAngTorsCoeffs, @bondBond13Coeffs);
my (@impropCoeffs, @angAngCoeffs);

my ($polymCutoff, $polymLink1, $polymLink2, $polymLinkNew1, $polymLinkNew2);
my ($polymCharge1, $polymCharge2, $polymClosest1, $polymClosest2);
my ($polymBondFlag, $polymAlignFlag, $polymChargeFlag, $polymIntraFlag);
my (@polymTypes, @polymConnect, @polymBonds, @polymPlanes, @polymVectors);

#
# Main Program
#

# Read command-line arguments
readArgs();

# Read input files
readLammps($fileData);
readTypes($fileTypes);
readPolymInput($fileInput);

# Add charges
addCharges();

# Write updated file
writeLammps($fileOut);

###############################################################################
# Subroutines

# errExit( $error )
# Exit program and print error
sub errExit
{
	printf "Error: %s\n", $_[0];
	exit 2;
}

# readArgs( )
# Read command line arguements for program
sub readArgs
{
	# Variables
	my (@args, $flag);
	@args = @ARGV;

	while(scalar(@args) > 0)
	{
		$flag = shift(@args);

		# Input file
		# LAMMPS file of system to run polymerization step
		if ($flag eq "-i")
		{
			$fileData = shift(@args);
			errExit("LAMMPS data file '$fileData' does not exist.")
				if (! -e $fileData);
		}

		# Types file (optional)
		# Text file specifying atom type numerical values and strings
		elsif ($flag eq "-t")
		{
			$fileTypes = shift(@args);
			errExit("Data types file '$fileTypes' does not exist.")
				if (! -e $fileTypes);
		}

		# Input script
		# Parameters for the polymerization step
		elsif ($flag eq "-s")
		{
			$fileInput = shift(@args);
			errExit("Input script '$fileInput' does not exist.")
				if (! -e $fileInput);
		}

		# Output file
		# LAMMPS file with updated connectivity after new bond(s)
		elsif ($flag eq "-o")
		{
			$fileOut = shift(@args);
		}

		# Help/syntax
		elsif ($flag eq "-h")
		{
			printf "Syntax: ./polym_init.pl -i data.lmps -t types.txt ".
				"-s polym.in -o new.lmps\n";

			exit 2;
		}

		# Error
		else
		{
			errExit("Did not recognize command-line flag.\n".
				"Syntax: ./polym_init.pl -i data.lmps -t types.txt ".
				"-s polym.in -o new.lmps\n");
		}
	}

	# Check values are defined
	errExit("Output file name not properly defined.") if (!defined($fileOut));
	errExit("Data file was not properly defined.") if (!defined($fileData));
	errExit("Types file was not properly defined.") if (!defined($fileTypes));
	errExit("Input script was not properly defined.") if (!defined($fileInput));
}

# addCharges( )
# Add charges to linker atoms
sub addCharges
{
	# Variables
	my ($type1, $type2);
	my (@link1, @link2);

	# Skip if no charges defined
	return if (!$polymChargeFlag);

	# Get atom types of linking atoms
	$type1 = $atomTypesKey{$polymLink1};
	$type2 = $atomTypesKey{$polymLink2};
	errExit("Atom type '$polymLink1' was not found.") if (!defined($type1));
	errExit("Atom type '$polymLink2' was not found.") if (!defined($type2));

	# Get groups of linking atoms
	@link1 = group(\@atomType, $type1);
	@link2 = group(\@atomType, $type2);
	errExit("No atoms of type '$polymLink1' found.") if (scalar(@link1) == 0);
	errExit("No atoms of type '$polymLink2' found.") if (scalar(@link2) == 0);

	# Add charges to linking atoms
	for (my $i = 0; $i < scalar(@link1); $i++) {
		$atomQ[$link1[$i]] += $polymCharge1;
	}
	for (my $i = 0; $i < scalar(@link2); $i++) {
		$atomQ[$link2[$i]] += $polymCharge2;
	}
}

# removeCharges( )
# Remove charges from linker atoms
sub removeCharges
{
	# Variables
        my ($type1, $type2);
        my (@link1, @link2);

        # Skip if no charges defined
        return if (!$polymChargeFlag);

        # Get atom types of linking atoms
        $type1 = $atomTypesKey{$polymLink1};
        $type2 = $atomTypesKey{$polymLink2};
        errExit("Atom type '$polymLink1' was not found.") if (!defined($type1));
        errExit("Atom type '$polymLink2' was not found.") if (!defined($type2));

        # Get groups of linking atoms
        @link1 = group(\@atomType, $type1);
        @link2 = group(\@atomType, $type2);
        errExit("No atoms of type '$polymLink1' found.") if (scalar(@link1) == 0);
        errExit("No atoms of type '$polymLink2' found.") if (scalar(@link2) == 0);

        # Remove charges from linking atoms
        for (my $i = 0; $i < scalar(@link1); $i++) {
                $atomQ[$link1[$i]] -= $polymCharge1;
        }
        for (my $i = 0; $i < scalar(@link2); $i++) {
                $atomQ[$link2[$i]] -= $polymCharge2;
        }
}

# group( @array, $value )
# Return all indices of @array that have a value that matches $value
# i.e., all $i such that $array[$i] = $value
sub group
{
	# Variables
	my (@array) = @{$_[0]};
	my $value = $_[1];
	my @matches;

	for (my $i = 0; $i < scalar(@array); $i++) {
		push(@matches, $i) if ($array[$i] == $value);
	}

	return @matches;
}

# readPolymInput( $file )
# Read input script for Polymatic polymerization step
sub readPolymInput
{
	# Variables
	my $file = $_[0];
	my ($command, $secConnect, $secTypes);
	my @params;

	# Initialize flags
	$polymBondFlag = 0;
	$polymAlignFlag = 0;
	$polymChargeFlag = 0;
	$polymIntraFlag = 0;
	$secConnect = 0;
	$secTypes = 0;

	open POLYMIN, "< $file" or die "Error opening file '$file': $!";

	while (my $line = <POLYMIN>)
	{
		chomp($line);
		$line =~ s/^\s+//;

		# Split line by spaces
		@params = split(' ', $line);
		$command = $params[0];

		# Section flags
		if ($command eq 'connect') {
			$secConnect = 1;
			$secTypes = 0;
			next;
		} elsif ($command eq 'types') {
			$secTypes = 1;
			$secConnect = 0;
			next;
		}

		# Connect records
		if ($secConnect) {
			$polymConnect[$command] = [ split(',', $params[1]) ];
		}

		# Type records
		elsif ($secTypes) {
			$polymTypes[$command] = $params[1];
		}

		# Cutoff
		elsif ($command eq 'cutoff') {
			$polymCutoff = $params[1];
		}

		# Linking atoms
		elsif ($command eq 'link') {
			($polymLink1, $polymLinkNew1) = split(',', $params[1]);
			($polymLink2, $polymLinkNew2) = split(',', $params[2]);
		}

		# Intramolecular bonding
		elsif ($command eq 'intra') {
			$polymIntraFlag = 1 if ($params[1] eq 'true');
		}

		# Linking atom charges
		elsif ($command eq 'charge') {
			$polymCharge1 = $params[1];
			$polymCharge2 = $params[2];
			$polymChargeFlag = 1;
		}

		# Additional bonds
		elsif ($command eq 'bond') {
			push( @polymBonds, [$params[1], $params[2]] );
			$polymBondFlag = 1;
		}

		# Plane alignment checks
		elsif ($command eq 'plane') {
			push ( @polymPlanes, [$params[1], $params[2], $params[3]] );
			$polymAlignFlag = 1;
		}

		# Vector alignment checks
		elsif ($command eq 'vector') {
			push ( @polymVectors, [$params[1], $params[2], $params[3]] );
			$polymAlignFlag = 1;
		}
	}

	close POLYMIN;

	# Check for required parameters
	errExit("Cutoff radius not properly defined.") if (!defined($polymCutoff));
	errExit("Reactive atoms not properly defined.")
		if (!defined($polymLink1) || !defined($polymLinkNew1)
			|| !defined($polymLink2) || !defined($polymLinkNew2));
}

# readTypes( $file )
# Read LAMMPS data types from a text types file
sub readTypes
{
	# Variables
	my $file = $_[0];
	my ($num, $string);

	# Section flags
	my $secAtom = 0;
	my $secBond = 0;
	my $secAngle = 0;
	my $secDihed = 0;
	my $secImprop = 0;

	open INTYPES, "< $file" or die "Error opening file '$file': $!";

	while (my $line = <INTYPES>)
	{
		chomp($line);
		$line =~ s/^\s+//;

		# Skip if commented line
		next if (substr($line,0,1) eq "#");

		# New section
		if (substr($line,0,10) eq "atom types" ||
			substr($line,0,10) eq "bond types" ||
			substr($line,0,11) eq "angle types" ||
			substr($line,0,14) eq "dihedral types" ||
			substr($line,0,14) eq "improper types")
		{
			$secAtom = 0;
			$secBond = 0;
			$secAngle = 0;
			$secDihed = 0;
			$secImprop = 0;

			$secAtom = 1 if (substr($line,0,10) eq "atom types");
			$secBond = 1 if (substr($line,0,10) eq "bond types");
			$secAngle = 1 if (substr($line,0,11) eq "angle types");
			$secDihed = 1 if (substr($line,0,14) eq "dihedral types");
			$secImprop = 1 if (substr($line,0,14) eq "improper types");

			next;
		}

		# Split line by space
		($num, $string) = split(' ', $line);

		# Atom section
		if ($secAtom) {
			$atomTypes[$num] = $string;
			$atomTypesKey{$string} = $num;
		}

		# Bond section
		elsif ($secBond) {
			$bondTypes[$num] = $string;
			$bondTypesKey{$string} = $num;
		}

		# Angle section
		elsif ($secAngle) {
			$angleTypes[$num] = $string;
			$angleTypesKey{$string} = $num;
		}

		# Dihedral section
		elsif ($secDihed) {
			$dihedTypes[$num] = $string;
			$dihedTypesKey{$string} = $num;
		}

		# Improper section
		elsif ($secImprop) {
			$impropTypes[$num] = $string;
			$impropTypesKey{$string} = $num;
		}
	}

	close INTYPES;

	# Check counts against data file
	errExit("Number of atom types in `$file' is not consistent.")
		if ($numAtomTypes != scalar(@atomTypes) - 1);
	errExit("Number of bond types in `$file' is not consistent.")
		if ($numBondTypes != scalar(@bondTypes) - 1);
	errExit("Number of angle types in `$file' is not consistent.")
		if ($numAngleTypes != scalar(@angleTypes) - 1);
	errExit("Number of dihedral types in `$file' is not consistent.")
		if ($numDihedTypes != scalar(@dihedTypes) - 1);
	errExit("Number of improper types in `$file' is not consistent.")
		if ($numImpropTypes != scalar(@impropTypes) - 1);
}

# readLammps( $file )
# Read LAMMPS data file
sub readLammps
{
	# Variables
	my $file = $_[0];
	my ($num, $mol, $type, $q, $x, $y, $z, $atom1, $atom2, $atom3, $atom4);
	my ($secMasses, $secPairCoeff, $secBondCoeff, $secAngleCoeff);
	my ($secBondBondCoeff, $secBondAngleCoeff, $secDihedCoeff);
	my ($secMidBondTorsCoeff, $secEndBondTorsCoeff, $secAngTorsCoeff);
	my ($secAngAngTorsCoeff, $secBondBond13Coeff, $secImpropCoeff);
	my ($secAngleAngleCoeff, $secAtoms, $secBonds, $secAngles);
	my ($secDiheds, $secImprops);
	my @temp;
	my $emptyLine = 2;
	my $bodyStart = 13;

	# Initialize
	$header = "";
	$lengthA = 0;
	$lengthB = 0;
	$lengthC = 0;
	$xLo = 0;
	$xHi = 0;
	$yLo = 0;
	$yHi = 0;
	$zLo = 0;
	$zHi = 0;
	$numAtoms = 0;
	$numBonds = 0;
	$numAngles = 0;
	$numDiheds = 0;
	$numImprops = 0;
	$numMols = 0;
	@atomMol = ();
	@atomType = ();
	@atomQ = ();
	@atomPos = ();
	@bonds = ();
	@angles = ();
	@diheds = ();
	@improps = ();

	# Open and read from file
	open INLMPS, "< $file" or die "Error opening file '$file': $!";

	my $i = 0;
	while (my $line = <INLMPS>)
	{
		chomp($line);
		$line =~ s/^\s+//;
		$i++;

		# Header
		if ($i == 1) {
			$header = $line;
			next;
		}

		# Counts
		elsif ($i == 3) {
			@temp = split(' ', $line);
			errExit("Atoms count not on proper line.")
				if ($temp[1] ne "atoms");
			$numAtoms = $temp[0];
			$bodyStart++ if ($numAtoms > 0);
			next;
		} elsif ($i == 4) {
			@temp = split(' ', $line);
			errExit("Bonds count not on proper line.")
				if ($temp[1] ne "bonds");
			$numBonds = $temp[0];
			$bodyStart++ if ($numBonds > 0);
			next;
		} elsif ($i == 5) {
			@temp = split(' ', $line);
			errExit("Angles count not on proper line.")
				if ($temp[1] ne "angles");
			$numAngles = $temp[0];
			$bodyStart++ if ($numAngles > 0);
			next;
		} elsif ($i == 6) {
			@temp = split(' ', $line);
			errExit("Dihedrals count not on proper line.")
				if ($temp[1] ne "dihedrals");
			$numDiheds = $temp[0];
			$bodyStart++ if ($numDiheds > 0);
			next;
		} elsif ($i == 7) {
			@temp = split(' ', $line);
			errExit("Impropers count not on proper line.")
				if ($temp[1] ne "impropers");
			$numImprops = $temp[0];
			$bodyStart++ if ($numImprops > 0);
			next;
		}

		# Data type counts and box size
		elsif ($i < $bodyStart)
		{
			@temp = split(' ', $line);

			# Save counts
			if ($temp[2] eq "types")
			{
				if ($temp[1] eq "atom") {
					$numAtomTypes = $temp[0];
				}
				elsif ($temp[1] eq "bond") {
					$numBondTypes = $temp[0];
				}
				elsif ($temp[1] eq "angle") {
					$numAngleTypes = $temp[0];
				}
				elsif ($temp[1] eq "dihedral") {
					$numDihedTypes = $temp[0];
				}
				elsif ($temp[1] eq "improper") {
					$numImpropTypes = $temp[0];
				}
			}

			elsif ($temp[2] eq "xlo") {
				$xLo = $temp[0];
				$xHi = $temp[1];
				$lengthA = $xHi - $xLo;
			} elsif ($temp[2] eq "ylo") {
				$yLo = $temp[0];
				$yHi = $temp[1];
				$lengthB = $xHi - $xLo;
			} elsif ($temp[2] eq "zlo") {
				$zLo = $temp[0];
				$zHi = $temp[1];
				$lengthC = $xHi - $xLo;
			}

			next;
		}

		# Flag sections
		if (substr($line,0,6) eq "Masses" ||
			substr($line,0,11) eq "Pair Coeffs" ||
			substr($line,0,11) eq "Bond Coeffs" ||
			substr($line,0,12) eq "Angle Coeffs" ||
			substr($line,0,15) eq "BondBond Coeffs" ||
			substr($line,0,16) eq "BondAngle Coeffs" ||
			substr($line,0,15) eq "Dihedral Coeffs" ||
			substr($line,0,24) eq "MiddleBondTorsion Coeffs" ||
			substr($line,0,21) eq "EndBondTorsion Coeffs" ||
			substr($line,0,19) eq "AngleTorsion Coeffs" ||
			substr($line,0,24) eq "AngleAngleTorsion Coeffs" ||
			substr($line,0,17) eq "BondBond13 Coeffs" ||
			substr($line,0,15) eq "Improper Coeffs" ||
			substr($line,0,17) eq "AngleAngle Coeffs" ||
			substr($line,0,5) eq "Atoms" ||
			substr($line,0,10) eq "Velocities" ||
			substr($line,0,5) eq "Bonds" ||
			substr($line,0,6) eq "Angles" ||
			substr($line,0,9) eq "Dihedrals" ||
			substr($line,0,9) eq "Impropers" )
		{
			$secMasses = 0;
			$secPairCoeff = 0;
			$secBondCoeff = 0;
			$secAngleCoeff = 0;
			$secBondBondCoeff = 0;
			$secBondAngleCoeff = 0;
			$secDihedCoeff = 0;
			$secMidBondTorsCoeff = 0;
			$secEndBondTorsCoeff = 0;
			$secAngTorsCoeff = 0;
			$secAngAngTorsCoeff = 0;
			$secBondBond13Coeff = 0;
			$secImpropCoeff = 0;
			$secAngleAngleCoeff = 0;
			$secAtoms = 0;
			$secBonds = 0;
			$secAngles = 0;
			$secDiheds = 0;
			$secImprops = 0;

			$secMasses = 1
				if (substr($line,0,6) eq "Masses");
			$secPairCoeff = 1
				if (substr($line,0,11) eq "Pair Coeffs");
			$secBondCoeff = 1
				if (substr($line,0,11) eq "Bond Coeffs");
			$secAngleCoeff = 1
				if (substr($line,0,12) eq "Angle Coeffs");
			$secBondBondCoeff = 1
				if (substr($line,0,15) eq "BondBond Coeffs");
			$secBondAngleCoeff = 1
				if (substr($line,0,16) eq "BondAngle Coeffs");
			$secDihedCoeff = 1
				if (substr($line,0,15) eq "Dihedral Coeffs");
			$secMidBondTorsCoeff = 1
				if (substr($line,0,24) eq "MiddleBondTorsion Coeffs");
			$secEndBondTorsCoeff = 1
				if (substr($line,0,21) eq "EndBondTorsion Coeffs");
			$secAngTorsCoeff = 1
				if (substr($line,0,19) eq "AngleTorsion Coeffs");
			$secAngAngTorsCoeff = 1
				if (substr($line,0,24) eq "AngleAngleTorsion Coeffs");
			$secBondBond13Coeff = 1
				if (substr($line,0,17) eq "BondBond13 Coeffs");
			$secImpropCoeff = 1
				if (substr($line,0,15) eq "Improper Coeffs");
			$secAngleAngleCoeff = 1
				if (substr($line,0,17) eq "AngleAngle Coeffs");
			$secAtoms = 1
				if (substr($line,0,5) eq "Atoms");
			$secBonds = 1
				if (substr($line,0,5) eq "Bonds");
			$secAngles = 1
				if (substr($line,0,6) eq "Angles");
			$secDiheds = 1
				if (substr($line,0,9) eq "Dihedrals");
			$secImprops = 1
				if (substr($line,0,9) eq "Impropers");

			next;
		}

		@temp = split(' ', $line);

		# Mass section
		if ($secMasses && length($line) > $emptyLine) {
			$num = shift(@temp);
			$masses[$num] = [@temp];
		}

		# Nonbond Coeffs section
		elsif ($secPairCoeff && length($line) > $emptyLine) {
			$num = shift(@temp);
			$pairCoeffs[$num] = [@temp];
		}

		# Bond Coeffs section
		elsif ($secBondCoeff && length($line) > $emptyLine) {
			$num = shift(@temp);
			$bondCoeffs[$num] = [@temp];
		}

		# Angle Coeffs section
		elsif ($secAngleCoeff && length($line) > $emptyLine) {
			$num = shift(@temp);
			$angleCoeffs[$num] = [@temp];
		}

		# Bond Bond Coeffs section
		elsif ($secBondBondCoeff && length($line) > $emptyLine) {
			$num = shift(@temp);
			$bondBondCoeffs[$num] = [@temp];
		}

		# Bond Angle Coeffs section
		elsif ($secBondAngleCoeff && length($line) > $emptyLine) {
			$num = shift(@temp);
			$bondAngCoeffs[$num] = [@temp];
		}

		# Dihedral Coeffs section
		elsif ($secDihedCoeff && length($line) > $emptyLine) {
			$num = shift(@temp);
			$dihedCoeffs[$num] = [@temp];
		}

		# Middle Bond Torsion Coeffs section
		elsif ($secMidBondTorsCoeff && length($line) > $emptyLine) {
			$num = shift(@temp);
			$midBondTorsCoeffs[$num] = [@temp];
		}

		# End Bond Torsion Coeffs section
		elsif ($secEndBondTorsCoeff && length($line) > $emptyLine) {
			$num = shift(@temp);
			$endBondTorsCoeffs[$num] = [@temp];
		}

		# Angle Torsion Coeffs section
		elsif ($secAngTorsCoeff && length($line) > $emptyLine) {
			$num = shift(@temp);
			$angTorsCoeffs[$num] = [@temp];
		}

		# Angle Angle Torsion Coeffs section
		elsif ($secAngAngTorsCoeff && length($line) > $emptyLine) {
			$num = shift(@temp);
			$angAngTorsCoeffs[$num] = [@temp];
		}

		# Bond Bond 1-3 Torsion Coeffs section
		elsif ($secBondBond13Coeff && length($line) > $emptyLine) {
			$num = shift(@temp);
			$bondBond13Coeffs[$num] = [@temp];
		}

		# Improper Coeffs section
		elsif ($secImpropCoeff && length($line) > $emptyLine) {
			$num = shift(@temp);
			$impropCoeffs[$num] = [@temp];
		}

		# Angle Angle Coeffs section
		elsif ($secAngleAngleCoeff && length($line) > $emptyLine) {
			$num = shift(@temp);
			$angAngCoeffs[$num] = [@temp];
		}

		# Atoms section
		elsif ($secAtoms && length($line) > $emptyLine)
		{
			($num, $mol, $type, $q, $x, $y, $z) = @temp;
			$atomMol[$num] = $mol;
			$atomType[$num] = $type;
			$atomQ[$num] = $q;
			$atomPos[$num] = [$x,$y,$z];

			$numMols = $mol if ($mol > $numMols);
			push ( @{$molecules[$mol]}, $num );
		}

		# Bonds section
		elsif ($secBonds && length($line) > $emptyLine)
		{
			($num, $type, $atom1, $atom2) = @temp;
			$bonds[$num] = [$type, $atom1, $atom2];

			push ( @{$atomBonds[$atom1]}, $atom2 );
			push ( @{$atomBonds[$atom2]}, $atom1 );
			push ( @{$atomBondNums[$atom1]}, $num );
			push ( @{$atomBondNums[$atom2]}, $num );
		}

		# Angles section
		elsif ($secAngles && length($line) > $emptyLine)
		{
			($num, $type, $atom1, $atom2, $atom3) = @temp;
			$angles[$num] = [$type, $atom1, $atom2, $atom3];

			push ( @{$atomAngles[$atom1]}, $atom2, $atom3 );
			push ( @{$atomAngles[$atom3]}, $atom2, $atom1 );
			push ( @{$atomAngleNums[$atom1]}, $num );
			push ( @{$atomAngleNums[$atom2]}, $num );
			push ( @{$atomAngleNums[$atom3]}, $num );
		}

		# Dihedrals section
		elsif ($secDiheds && length($line) > $emptyLine)
		{
			($num, $type, $atom1, $atom2, $atom3, $atom4) = @temp;
			$diheds[$num] = [$type, $atom1, $atom2, $atom3, $atom4];

			push ( @{$atomDiheds[$atom1]}, $atom2, $atom3, $atom4 );
			push ( @{$atomDiheds[$atom4]}, $atom3, $atom2, $atom1 );
			push ( @{$atomDihedNums[$atom1]}, $num );
			push ( @{$atomDihedNums[$atom2]}, $num );
			push ( @{$atomDihedNums[$atom3]}, $num );
			push ( @{$atomDihedNums[$atom4]}, $num );
		}

		# Impropers section
		elsif ($secImprops && length($line) > $emptyLine)
		{
			($num, $type, $atom1, $atom2, $atom3, $atom4) = @temp;
			$improps[$num] = [$type, $atom1, $atom2, $atom3, $atom4];

			push ( @{$atomImprops[$atom2]}, $atom1, $atom3, $atom4 );
			push ( @{$atomImpropNums[$atom1]}, $num );
			push ( @{$atomImpropNums[$atom2]}, $num );
			push ( @{$atomImpropNums[$atom3]}, $num );
			push ( @{$atomImpropNums[$atom4]}, $num );
		}
	}

	close INLMPS;

	# Check for errors
	errExit("Length of atom data does not match atom count.")
		if (scalar(@atomMol) - 1 != $numAtoms ||
			scalar(@atomType) - 1 != $numAtoms ||
			scalar(@atomQ) - 1 != $numAtoms ||
			scalar(@atomPos) - 1 != $numAtoms);

	errExit("Length of bond data does not match bond count.")
		if (@bonds && scalar(@bonds) - 1 != $numBonds);

	errExit("Length of angle data does not match angle count.")
		if (@angles && scalar(@angles) - 1 != $numAngles);

	errExit("Length of dihedral data does not match dihedral count.")
		if (@diheds && scalar(@diheds) - 1 != $numDiheds);

	errExit("Length of improper data does not match improper count.")
		if (@improps && scalar(@improps) - 1 != $numImprops);

	errExit("Length of atom types data does not match atom types count.")
		if (scalar(@masses) - 1 != $numAtomTypes ||
			scalar(@pairCoeffs) - 1 != $numAtomTypes);

	errExit("Length of bond types data does not match bond types count.")
		if (@bondCoeffs && scalar(@bondCoeffs) - 1 != $numBondTypes);

	errExit("Length of angle types data does not match angle types count.")
		if (@angleCoeffs && scalar(@angleCoeffs) - 1 != $numAngleTypes);

	errExit("Length of dihedral types data does not match dihedral types count.")
		if (@dihedCoeffs && scalar(@dihedCoeffs) - 1 != $numDihedTypes);

	errExit("Length of improper types data does not match improper types count.")
		if (@impropCoeffs && scalar(@impropCoeffs) - 1 != $numImpropTypes);
}

# writeLammps( $file )
# Write LAMMPS data file
sub writeLammps
{
	# Variables
	my $file = $_[0];
	my @values;

	# Check for necessary information
	errExit("LAMMPS file cannot be written, box dimensions not defined.")
		if (!defined($xLo) || !defined($xHi) ||
			!defined($yLo) || !defined($yHi) ||
			!defined($zLo) || !defined($zHi));

	errExit("Cannot write LAMMPS file, atoms not defined properly.")
		if ($numAtoms == 0 ||
			scalar(@atomMol)-1 != $numAtoms ||
			scalar(@atomType)-1 != $numAtoms ||
			scalar(@atomQ)-1 != $numAtoms ||
			scalar(@atomPos)-1 != $numAtoms);

	errExit("Cannot write LAMMPS file, bonds not defined properly.")
		if ($numBonds > 0 && scalar(@bonds)-1 != $numBonds);

	errExit("Cannot write LAMMPS file, angles not defined properly.")
		if ($numAngles > 0 && scalar(@angles)-1 != $numAngles);

	errExit("Cannot write LAMMPS file, dihedrals not defined properly.")
		if ($numDiheds > 0 && scalar(@diheds)-1 != $numDiheds);

	errExit("Cannot write LAMMPS file, impropers not defined properly.")
		if ($numImprops > 0 && scalar(@improps)-1 != $numImprops);

	errExit("Cannot write LAMMPS file, atom types not defined properly.")
		if ($numAtomTypes == 0 || scalar(@masses)-1 != $numAtomTypes ||
			scalar(@pairCoeffs)-1 != $numAtomTypes);

	errExit("Cannot write LAMMPS file, bond types not defined properly.")
		if ($numBondTypes > 0 && scalar(@bondCoeffs)-1 != $numBondTypes);

	errExit("Cannot write LAMMPS file, angle types not defined properly.")
		if ($numAngleTypes > 0 && scalar(@angleCoeffs)-1 != $numAngleTypes);

	errExit("Cannot write LAMMPS file, dihedral types not defined properly.")
		if ($numDihedTypes > 0 && scalar(@dihedCoeffs)-1 != $numDihedTypes);

	errExit("Cannot write LAMMPS file, bond types not defined properly.")
		if ($numImpropTypes > 0 &&
			scalar(@impropCoeffs)-1 != $numImpropTypes);

	# Open and write to file
	open FILE, "> $file" or die "Error opening file '$file': $!";

	# Heading
	printf FILE "$header\n";
	printf FILE "\n";

	printf FILE "%d atoms\n", $numAtoms;
	printf FILE "%d bonds\n", $numBonds;
	printf FILE "%d angles\n", $numAngles;
	printf FILE "%d dihedrals\n", $numDiheds;
	printf FILE "%d impropers\n", $numImprops;
	printf FILE "\n";

	printf FILE "%d atom types\n", $numAtomTypes
		if ($numAtomTypes > 0);
	printf FILE "%d bond types\n", $numBondTypes
		if ($numBondTypes > 0);
	printf FILE "%d angle types\n", $numAngleTypes
		if ($numAngleTypes > 0);
	printf FILE "%d dihedral types\n", $numDihedTypes
		if ($numDihedTypes > 0);
	printf FILE "%d improper types\n", $numImpropTypes
		if ($numImpropTypes > 0);
	printf FILE "\n";

	printf FILE "%10.6f %10.6f xlo xhi \n", $xLo, $xHi;
	printf FILE "%10.6f %10.6f ylo yhi \n", $yLo, $yHi;
	printf FILE "%10.6f %10.6f zlo zhi \n", $zLo, $zHi;
	printf FILE "\n";

	# Atom Coeffs
	if ($numBondTypes > 0)
	{
		# Masses
		printf FILE "Masses\n\n";
		for (my $i = 1; $i <= $numAtomTypes; $i++)
		{
			errExit("Did not find mass for atom type $i.")
				if (!$masses[$i]);
			@values = @{$masses[$i]};

			printf FILE "  %-7d %11.6f\n",
				$i, $values[0];
		}
		printf FILE "\n";

		# Pair Coeffs
		printf FILE "Pair Coeffs\n\n";
		for (my $i = 1; $i <= $numAtomTypes; $i++)
		{
			errExit("Did not find pair coeffs for atom type $i.")
				if (!$pairCoeffs[$i]);
			@values = @{$pairCoeffs[$i]};

			printf FILE "  %-7d %11.6f %11.6f\n",
				$i, $values[0], $values[1];
		}
		printf FILE "\n";
	}

	# Bond Coeffs
	if ($numBondTypes > 0)
	{
		# Bond Coeffs
		printf FILE "Bond Coeffs\n\n";
		for (my $i = 1; $i <= $numBondTypes; $i++)
		{
			errExit("Did not find bond coeffs for bond type $i.")
				if (!$bondCoeffs[$i]);
			@values = @{$bondCoeffs[$i]};

			printf FILE "  %-7d %11.6f %11.6f %11.6f %11.6f\n",
				$i, $values[0], $values[1], $values[2], $values[3];
		}
		printf FILE "\n";
	}

	# Angle Coeffs
	if ($numAngleTypes > 0)
	{
		# Angle Coeffs
		printf FILE "Angle Coeffs\n\n";
		for (my $i = 1; $i <= $numAngleTypes; $i++)
		{
			errExit("Did not find angle coeffs for angle type $i.")
				if (!$angleCoeffs[$i]);
			@values = @{$angleCoeffs[$i]};

			printf FILE "  %-7d %11.6f %11.6f %11.6f %11.6f\n",
				$i, $values[0], $values[1], $values[2], $values[3];
		}
		printf FILE "\n";

		# BondBond Coeffs
=a
		printf FILE "BondBond Coeffs\n\n";
		for (my $i = 1; $i <= $numAngleTypes; $i++)
		{
			errExit("Did not find bond bond coeffs for angle type $i.")
				if (!$bondBondCoeffs[$i]);
			@values = @{$bondBondCoeffs[$i]};

			printf FILE "  %-7d %11.6f %11.6f %11.6f\n",
				$i, $values[0], $values[1], $values[2];
		}
		printf FILE "\n";

		# BondAngle Coeffs
		printf FILE "BondAngle Coeffs\n\n";
		for (my $i = 1; $i <= $numAngleTypes; $i++)
		{
			errExit("Did not find bond angle coeffs for angle type $i.")
				if (!$bondAngCoeffs[$i]);
			@values = @{$bondAngCoeffs[$i]};

			printf FILE "  %-7d %11.6f %11.6f %11.6f %11.6f\n",
				$i, $values[0], $values[1], $values[2], $values[3];
		}
		printf FILE "\n";
=cut
	}

	if ($numDihedTypes > 0)
	{
		# Dihed Coeffs
		printf FILE "Dihedral Coeffs\n\n";
		for (my $i = 1; $i <= $numDihedTypes; $i++)
		{
			errExit("Did not find dihedral coeffs for dihedral type $i.")
				if (!$dihedCoeffs[$i]);
			@values = @{$dihedCoeffs[$i]};

			printf FILE
				"  %-7d %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f\n",
				$i, $values[0], $values[1], $values[2], $values[3], $values[4],
				$values[5];
		}
		printf FILE "\n";

		# MiddleBondTorsion Coeffs
=a
		printf FILE "MiddleBondTorsion Coeffs\n\n";
		for (my $i = 1; $i <= $numDihedTypes; $i++)
		{
			errExit("Did not find mbt coeffs for dihedral type $i.")
				if (!$midBondTorsCoeffs[$i]);
			@values = @{$midBondTorsCoeffs[$i]};

			printf FILE "  %-7d %11.6f %11.6f %11.6f %11.6f\n",
				$i, $values[0], $values[1], $values[2], $values[3];
		}
		printf FILE "\n";

		# EndBondTorsion Coeffs
		printf FILE "EndBondTorsion Coeffs\n\n";
		for (my $i = 1; $i <= $numDihedTypes; $i++)
		{
			errExit("Did not find ebt coeffs for dihedral type $i.")
				if (!$endBondTorsCoeffs[$i]);
			@values = @{$endBondTorsCoeffs[$i]};

			printf FILE "  %-7d %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f ".
				"%11.6f %11.6f\n",
				$i, $values[0], $values[1], $values[2], $values[3], $values[4],
				$values[5], $values[6], $values[7];
		}
		printf FILE "\n";

		# AngleTorsion Coeffs
		printf FILE "AngleTorsion Coeffs\n\n";
		for (my $i = 1; $i <= $numDihedTypes; $i++)
		{
			errExit(
				"Did not find angle torsion coeffs for dihedral type $i.")
				if (!$angTorsCoeffs[$i]);
			@values = @{$angTorsCoeffs[$i]};

			printf FILE "  %-7d %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f ".
				"%11.6f %11.6f\n",
				$i, $values[0], $values[1], $values[2], $values[3], $values[4],
				$values[5], $values[6], $values[7];
		}
		printf FILE "\n";

		# AngleAngleTorsion Coeffs
		printf FILE "AngleAngleTorsion Coeffs\n\n";
		for (my $i = 1; $i <= $numDihedTypes; $i++)
		{
			errExit("Did not find angle angle torsion coeffs for dihedral ".
				"type $i.")
				if (!$angAngTorsCoeffs[$i]);
			@values = @{$angAngTorsCoeffs[$i]};

			printf FILE "  %-7d %11.6f %11.6f %11.6f\n",
				$i, $values[0], $values[1], $values[2];
		}
		printf FILE "\n";

		# BondBond13 Coeffs
		printf FILE "BondBond13 Coeffs\n\n";
		for (my $i = 1; $i <= $numDihedTypes; $i++)
		{
			errExit("Did not find bond bond 13 coeffs for dihedral type $i.")
				if (!$bondBond13Coeffs[$i]);
			@values = @{$bondBond13Coeffs[$i]};

			printf FILE "  %-7d %11.6f %11.6f %11.6f\n",
				$i, $values[0], $values[1], $values[2];
		}
		printf FILE "\n";
=cut
	}


	if ($numImpropTypes > 0)
	{
		# Improper Coeffs
		printf FILE "Improper Coeffs\n\n";
		for (my $i = 1; $i <= $numImpropTypes; $i++)
		{
			errExit("Did not find improper coeffs for improper type $i.")
				if (!$impropCoeffs[$i]);
			@values = @{$impropCoeffs[$i]};

			printf FILE "  %-7d %11.6f %11.6f\n",
				$i, $values[0], $values[1];
		}
		printf FILE "\n";

		# AngleAngle Coeffs
=a
		printf FILE "AngleAngle Coeffs\n\n";
		for (my $i = 1; $i <= $numImpropTypes; $i++)
		{
			errExit("Did not find angle angle coeffs for improper type $i.")
				if (!$angAngCoeffs[$i]);
			@values = @{$angAngCoeffs[$i]};

			printf FILE "  %-7d %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f\n",
				$i, $values[0], $values[1], $values[2], $values[3], $values[4],
				$values[5];
		}
		printf FILE "\n";
=cut
	}


	if ($numAtoms > 0)
	{
		# Atoms
		printf FILE "Atoms\n\n";
		for (my $i = 1; $i <= $numAtoms; $i++)
		{
			errExit("Did not find atom '$i'.")
				if (!$atomPos[$i] || !$atomMol[$i] || !$atomType[$i] ||
					!defined($atomQ[$i]));
			@values = @{$atomPos[$i]};

			printf FILE "  %-7d %-4d %-4d %10.6f %13.6f %13.6f %13.6f\n",
				$i, $atomMol[$i], $atomType[$i], $atomQ[$i],
				$values[0], $values[1], $values[2];
		}
		printf FILE "\n";
	}

	if ($numBonds > 0)
	{
		# Bonds
		printf FILE "Bonds\n\n";
		for (my $i = 1; $i <= $numBonds; $i++)
		{
			errExit("Did not find bond '$i'.") if (!$bonds[$i]);
			@values = @{$bonds[$i]};

			printf FILE "  %-7d %-4d %-7d %-7d\n",
				$i, $values[0], $values[1], $values[2];
		}
		printf FILE "\n";
	}

	if ($numAngles > 0)
	{
		# Angles
		printf FILE "Angles\n\n";
		for (my $i = 1; $i <= $numAngles; $i++)
		{
			errExit("Did not find angle '$i'.") if (!$angles[$i]);
			@values = @{$angles[$i]};

			printf FILE "  %-7d %-4d %-7d %-7d %-7d\n",
				$i, $values[0], $values[1], $values[2], $values[3];
		}
		printf FILE "\n";
	}

	if ($numDiheds > 0)
	{
		# Dihedrals
		printf FILE "Dihedrals\n\n";
		for (my $i = 1; $i <= $numDiheds; $i++)
		{
			errExit("Did not find dihedral '$i'.") if (!$diheds[$i]);
			@values = @{$diheds[$i]};

			printf FILE "  %-7d %-4d %-7d %-7d %-7d %-7d\n",
				$i, $values[0], $values[1], $values[2], $values[3], $values[4];
		}
		printf FILE "\n";
	}

	if ($numImprops > 0)
	{
		# Impropers
		printf FILE "Impropers\n\n";
		for (my $i = 1; $i <= $numImprops; $i++)
		{
			errExit("Did not find improper '$i'.") if (!$improps[$i]);
			@values = @{$improps[$i]};

			printf FILE "  %-7d %-4d %-7d %-7d %-7d %-7d\n",
				$i, $values[0], $values[1], $values[2], $values[3], $values[4];
		}
		printf FILE "\n";
	}
	# Close file
	close FILE;
}
