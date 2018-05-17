package Build::ChemNamelist;

#-------------------------------------------------------------------------------------
# generates species lists for chemistry namelist settings
#-------------------------------------------------------------------------------------
# Date         Contributor      Modification
#-------------------------------------------------------------------------------------
# 26 Jan 2011  Francis Vitt     Created 
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

use strict;
use Exporter;
use FindBin qw($Bin);
use lib "$Bin/perl5lib";
use Build::ChemPreprocess qw(get_species_list);

our @ISA = qw(Exporter);
our @EXPORT = qw(set_dep_lists set_aero_modes_info chem_has_species);
our $VERSION = 1.00;

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub chem_has_species
{
    my ( $cfg, $species ) = @_;
    my $chem_proc_src = $cfg->get('chem_proc_src');
    my $chem_src_dir = $cfg->get('chem_src_dir');
    my @species_list;
    if ($chem_proc_src) {
        @species_list = get_species_list($chem_proc_src);
    } else {
        @species_list = get_species_list($chem_src_dir);
    }
    my %hash;
    @hash{@species_list}=();

    return ( exists $hash{$species} );

}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub set_dep_lists
{
    my ( $cfgdir, $chem_proc_src, $chem_src_dir, $nl, $print_lvl ) = @_;

    my ( $gas_wetdep_list, $aer_wetdep_list, $aer_drydep_list, $aer_sol_facti, $aer_sol_factb, 
         $aer_scav_coef, $gas_drydep_list ) ;

    my @species_list ;
    if ($chem_proc_src) {
	if (defined $ENV{CASEBUILD}) {
            #needed to expand $CASEBUILD in $chem_proc_src for CESM scripts
	    my $root = $ENV{CASEBUILD};
            $chem_proc_src =~ s/\$CASEBUILD/$root/;
	}
	@species_list = get_species_list($chem_proc_src);
    } else {
	if (defined $ENV{CODEROOT}) {
            #needed to expand $CODEROOT in $chem_src_dir for CESM scripts
	    my $root = $ENV{CODEROOT};
            $chem_src_dir =~ s/\$CODEROOT/$root/;
	}
	@species_list = get_species_list($chem_src_dir);
    }
    if ($print_lvl>=2) {print "Chemistry species : @species_list \n" ;}

    $gas_wetdep_list = get_gas_wetdep_list( $cfgdir, $print_lvl, @species_list );

    $aer_wetdep_list = get_aer_wetdep_list( $cfgdir, $print_lvl, @species_list );

    $gas_drydep_list = get_gas_drydep_list( $cfgdir, $print_lvl, @species_list );

    $aer_drydep_list = get_aer_drydep_list( $cfgdir, $print_lvl, @species_list );

    # set solubility factors for aerosols
    if (length($aer_wetdep_list)>2){ 
        my @values = split(',', $aer_wetdep_list);
        my $first = 1; my $pre = "";
        foreach my $val (@values){
            $val =~ tr/'//d;
            my $sol_facti;
            my $sol_factb;
            if ($val =~ /DST/) {
                $sol_facti = $nl->get_value('dust_sol_facti');
                if (!defined $sol_facti) { $sol_facti = 0.3; }
                $sol_factb = $nl->get_value('dust_sol_factb');
                if (!defined $sol_factb) { $sol_factb = 0.3; }
            } elsif ($val =~ /SSLT/) {
                $sol_facti = $nl->get_value('sslt_sol_facti');
                if (!defined $sol_facti) { $sol_facti = 0.3; }
                $sol_factb = $nl->get_value('sslt_sol_factb');
                if (!defined $sol_factb) { $sol_factb = 0.3; }
            } else {
                $sol_facti = $nl->get_value(${val}.'_sol_facti');
                if (!defined $sol_facti) { $sol_facti = 0.3; }
                $sol_factb = $nl->get_value(${val}.'_sol_factb');
                if (!defined $sol_factb) { $sol_factb = 0.3; }
            }
            $aer_sol_facti .= $pre . $sol_facti ;
            $aer_sol_factb .= $pre . $sol_factb ;

            my $scav_coef = $nl->get_value(${val}.'_scav_coef');
            if (!defined $scav_coef) {
                if ($val =~ /DST03/ or $val =~ /DST04/) { $scav_coef = 0.3; }
                else { $scav_coef = 0.1; }
            }
            $aer_scav_coef .= $pre . $scav_coef ;

            if ($first) { $pre = ","; $first = 0; }
        }
    }

    return ( $gas_wetdep_list, $aer_wetdep_list, $aer_sol_facti, $aer_sol_factb, $aer_scav_coef, 
             $aer_drydep_list, $gas_drydep_list );
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub set_aero_modes_info
{
    my ( $cfg, $data_src, $print_lvl, $mode_types, $modal_species, $modal_groups, 
         $mode_spec_type, $mode_spec, $mode_spec_cw, $mode_spec_src ) = @_;

    my $chem_proc_src = $cfg->get('chem_proc_src');
    my $chem_src_dir = $cfg->get('chem_src_dir');
    my @species_list;
    if ($chem_proc_src) {
        @species_list = get_species_list($chem_proc_src);
    } else {
        @species_list = get_species_list($chem_src_dir);
    }

    my %mymodal_species = %$modal_species;
    my $nmodes = scalar(@$mode_types);
    my $found_modal_species = 0;

    for (my $i = 1; $i <= $nmodes; $i++) {
        my @type ;
        my @spec ;
        my @spec_cw ;
        my @src ;
        my $modal_species = 0;
        my @species = @{ $$modal_groups{ @$mode_types[$i-1] } } ;
        foreach my $spc (@species) {
            if ($data_src eq 'A') {
              foreach my $tracer (@species_list) {
                if ($tracer =~ /^${spc}.*_a${i}$/) {
                    $found_modal_species = 1;
                    push @spec, $tracer ;
                    push @src , $data_src ;
                    my $tracer_cw = $tracer; $tracer_cw =~ s/_a/_c/g;
                    push @spec_cw,  $tracer_cw;
                    push @type, $mymodal_species{$spc};
                }
              }
            } else { # for prescribed modal aerosols do not check against the species list ....
                my $tracer = ${spc}.'_a'.${i};
                $found_modal_species = 1;
                push @spec, $tracer ;
                push @src , $data_src ;
                my $tracer_cw = $tracer; $tracer_cw =~ s/_a/_c/g;
                push @spec_cw,  $tracer_cw;
                push @type, $mymodal_species{$spc};
            }
        }
        if ($found_modal_species) {
            push @$mode_spec, [ @spec ];
            push @$mode_spec_cw, [ @spec_cw ];
            push @$mode_spec_type, [ @type ];
            push @$mode_spec_src, [ @src ];
        }
    }

    if ($print_lvl>=2) { print_modal_info( "mode_spec", @$mode_spec ); }
    if ($print_lvl>=2) { print_modal_info( "mode_spec_cw", @$mode_spec_cw ); }
    if ($print_lvl>=2) { print_modal_info( "mode_spec_type", @$mode_spec_type ); }
    if ($print_lvl>=2) { print_modal_info( "mode_spec_src", @$mode_spec_src ); }

    return ;
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub print_modal_info
{
    my ( $string,  @array ) = @_;
    print "-------------------------------------------------------------------\n";
    print "$string :\n";
    foreach my $row (@array) {
	print "   ";
	foreach my $item (@$row) {
	    print $item, " ";
	}
	print "\n";
    }
    print "-------------------------------------------------------------------\n";
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub get_gas_drydep_list
{
    my ($cfg_dir,$print_lvl,@species_list) = @_;

    my $master_file = "$cfg_dir/namelist_files/master_gas_drydep_list.xml";

    my $list = get_dep_list($master_file,$print_lvl,@species_list);

    if ($print_lvl>=2) {print " dry dep list : $list  \n" ;}

    return ($list);
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub get_aer_drydep_list
{
    my ($cfg_dir,$print_lvl,@species_list) = @_;

    my $master_file = "$cfg_dir/namelist_files/master_aer_drydep_list.xml";

    my $list = get_dep_list($master_file,$print_lvl,@species_list);

    if ($print_lvl>=2) {print " aer drydep list : $list  \n" ;}
    return ($list);
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub get_aer_wetdep_list
{
    my ($cfg_dir,$print_lvl,@species_list) = @_;

    my $master_file = "$cfg_dir/namelist_files/master_aer_wetdep_list.xml";

    my $list = get_dep_list($master_file,$print_lvl,@species_list);

    if ($print_lvl>=2) {print " aer wet dep list : $list  \n" ;}
    return ($list);
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub get_gas_wetdep_list
{
    my ($cfg_dir,$print_lvl,@species_list) = @_;

    my $master_file = "$cfg_dir/namelist_files/master_gas_wetdep_list.xml";

    my $list = get_dep_list($master_file,$print_lvl,@species_list);

    if ($print_lvl>=2) {print " gas wet dep list : $list  \n" ;}

    return ($list);
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
sub get_dep_list
{
    my ($master_file,$print_lvl,@species_list) = @_;

    if ($print_lvl>=2){ print "Using chemistry master list file $master_file \n"; }

    my @master_list = read_master_list_file($master_file);

    my $list = '';
    my $first = 1; my $pre = "";
    foreach my $name (sort @species_list) {
	foreach my $item (@master_list) {
	    if ($name eq $item) { 
		$list .= $pre .  quote_string($name) ;
                if ($first) { $pre = ","; $first = 0; }
	    }
	}
    }

    if ( length($list)<1 ) {$list = quote_string(' ') ;}

    return ($list);
}

#-------------------------------------------------------------------------------
sub read_master_list_file
{
    my ($master_file) = @_;

    require XML::Lite;

    my @master_list ;
    my $xml = XML::Lite->new($master_file);
    my $root = $xml->root_element();
    my @children = $root->get_children();
    foreach my $child (@children) {
	my $content = $child->get_content();
	my @list = split( ('\s+|\s*,+\s*') ,$content);
	foreach my $item (@list) {
	    if ( length( $item) > 0 ){
		push ( @master_list, $item );
	    }
	}
    }

    return (@master_list);
}
#-------------------------------------------------------------------------------

sub quote_string {
    my $str = shift;
    $str =~ s/^\s+//;
    $str =~ s/\s+$//;
    unless ($str =~ /^['"]/) {        #"'
        $str = "\'$str\'";
    }
    return $str;
}
1; # to appease require 
