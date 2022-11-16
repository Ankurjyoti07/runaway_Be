#!/usr/bin/env perl
use strict; # recommended
use 5.16.0; # recommended
use binary_grid2; # required
use binary_grid::C; # backend : C or Perl
use rob_misc qw/ncpus/;
use Data::Dumper;
use IO::Handle;
use File::Copy;

############################################################
#
# Example script to demonstrate how to use the
# binary_grid2 module.
#
# For full documentation, please see binary_grid2.pdf
# in the doc/ directory of binary_c`
#
############################################################

# my $outdir ='/home/alexkemp/Desktop/PhDResearchStuff/binary_c_stuff/dmlossdmh_standard3_Hnovalowres_deriv0/';
# my $outdir ='/home/alexkemp/data/data/binary_c/osevent_dmlossdmh_ts1p0mod1000/';
# my $outdir ='/home/alexkemp/data/data/binary_c/FIXEDmaybe_standard3_Hnovahighres/';
# my $outdir ='/home/alexkemp/Desktop/PhDResearchStuff/binary_c_stuff/fixednovatimesteps_dmlossdmhcrit_ts0p001/';
# my $outdir = '/home/alexkemp/data/data/binary_c/DUMMY/';
my $outdir ='/home/alexkemp/data/data/binary_c/std4_metallicity0p0001_RlowHnova2lin/';

mkdir($outdir);
unlink("$outdir/"."grid.txt");
copy('/home/alexkemp/progs/stars/binary_c/src/perl/scripts2/kemplate.pl',$outdir);
# number of computational threads to launch
my $nthreads = rob_misc::ncpus();


############################################################
# Binary_c should output data that we can understand. There
# are two ways to do this:
#
# 1) Put output statements, using PRINTF, in binary_c's
#    log_every_timestep() function. This requires a rebuild
#    of libbinary_c.so and a reinstall of the binary_grid module
#    every time you change the PRINTF statement.
#
#
# 2) Put a list of hashes in the C_auto_logging grid option
#
#  $population->set(
#    C_auto_logging => {
#        'MY_STELLAR_DATA' =>
#            [
#             'model.time',
#             'star[0].mass',
#             'model.probability',
#             'model.dt'
#            ]
#    });
#
#  where MY_STELLAR_DATA is the key of the hash {...} and is also 
#  the header matched in the parse_data function (below). The list [...] 
#  contains the names of the variables to be output, which are all
#  assumed to be in stardata.
#
#  This option does not require a rebuild of libbinary_c.so or a
#  reinstall of binary_grid.
#
#
# 3) Put output statements, using PRINTF, into the C_logging_code
#    grid option
#
#  $population->set( C_logging_code => ' PRINTF("...\n"); ');
#
#    You have access to the stardata variable, so you can
#    output everything that is available to log_every_timestep();
#
#  This option does not require a rebuild of libbinary_c.so or a
#  reinstall of binary_grid.
#
############################################################
# make a new stellar population



# make a new stellar population
my $population = binary_grid2->new(
    # options can be given now ... 
    metallicity => 0.0001, # mass fraction of "metals"
    max_evolution_time => 15000, # Myr 
    nthreads=>$nthreads, # number of threads
    );

# ... or options can be set manually later.
$population->set(
    vb=>1, # turn on verbose logging (can be 0,1,2...)
    return_array_refs=>1, # quicker data parsing mode
    log_args=>1,
    sort_args=>1,
    save_args=>1,
    log_args_dir=>'/tmp',
    individual_novae=>1,
    nova_retention_method => 'NOVA_RETENTION_ALGORITHM_WANGWU',
    nova_retention_fraction => 0.0,
    WD_accretion_rate_novae_upper_limit_hydrogen_donor => 'DONOR_RATE_ALGORITHM_WANGWU',
    WD_accretion_rate_novae_upper_limit_helium_donor => 'DONOR_RATE_ALGORITHM_WANGWU',
    WD_accretion_rate_novae_upper_limit_other_donor => 'DONOR_RATE_ALGORITHM_WANGWU',
    WD_accretion_rate_new_giant_envelope_lower_limit_hydrogen_donor => 'DONOR_RATE_ALGORITHM_WANGWU',
    WD_accretion_rate_new_giant_envelope_lower_limit_helium_donor => 'DONOR_RATE_ALGORITHM_WANGWU',
    WD_accretion_rate_new_giant_envelope_lower_limit_other_donor => 'DONOR_RATE_ALGORITHM_WANGWU',
    nova_timestep_accelerator_num => -1,
    mass_accretion_for_eld => -1,
    WDWD_merger_algorithm => 'WDWD_MERGER_ALGORITHM_HYBRID_PERETS2019_SATO2016',
    alpha_ce=>1,
    lambda_ce=>'LAMBDA_CE_WANG_2016',
    AGB_core_algorithm=>2,
    AGB_radius_algorithm=>2,
    AGB_luminosity_algorithm=>2,
    AGB_3dup_algorithm=>2,
    accretion_limit_eddington_steady_multiplier=>-1,#wind or RLOF, default should be -1?
    accretion_limit_eddington_WD_to_remnant_multiplier=>-1,
    #unstable RLOF (from a WD) onto a WD, NS or BH, default should be -1 (ie. ignore the limit)
    accretion_limit_eddington_LMMS_multiplier=>1,
    #unstable RLOF from a convective, low mass main sequence, default should be 1 
    chandrasekhar_mass=>1.38,
    RLOF_method=>'RLOF_METHOD_CLAEYS',
    accretion_limit_thermal_multiplier=>10, #clayes et al 2014, eq 14
    wind_angular_momentum_loss=>'WIND_ANGMOM_LOSS_SPHERICALLY_SYMMETRIC', #sepherically symmetric.
    RLOF_angular_momentum_transfer_model=>'RLOF_ANGULAR_MOMENTUM_TRANSFER_MODEL_CONSERVATIVE',
    nonconservative_angmom_gamma=>'RLOF_NONCONSERVATIVE_GAMMA_ISOTROPIC',
    hachisu_disk_wind=>'True',
    wRLOF_method=>'WRLOF_Q_DEPENDENT',
    max_neutron_star_mass=>2.5, #2.5
    BH_prescription=>'BH_FRYER12_RAPID',
    nova_eta_shift=>0,
    #eta_ruiter2013=>0.7,
    # #log_dt_secs=>1,
    # C_auto_logging => {
    #     'MY_STELLAR_DATA' =>
    #         [
    #          'model.time',
    #          'star[0].mass',
    #          'model.probability',
    #          'model.dt'
    #         ]
    # },
 #   C_logging_code => 
 #       'PRINTF("stardata = %p\n",stardata);'
 #   ,
 
## or enter more complicated code yourself:
#       
#    C_logging_code => '
#             PRINTF("MY_STELLAR_DATA %g %g %g %g\n",
#                 stardata->model.time,
#                 stardata->star[0].mass,
#                 stardata->model.probability,
#                 stardata->model.dt);
#                       '
    );

# scan command line arguments for extra options
$population->parse_args();     

# duplicity is 0 for single stars, 1 for binary stars
# and 2 for a mixed population sampled at equal times
my $duplicity = 1;

if($duplicity == 0)
{
    # make a grid of $nstars single binary stars, log-spaced,
    # with masses between $mmin and $mmax
    my $nstars = 100;
    my $mmin = 0.1;
    my $mmax = 80.0;
    $population->add_grid_variable(
        'name'       => 'lnm1',
        'longname'   =>'Primary mass',
        'range'      =>[log($mmin),log($mmax)],
        'resolution' => $nstars, # just a counter for the grid
        'spacingfunc'=>"const(log($mmin),log($mmax),$nstars)",
        'precode'    =>'$m1=exp($lnm1);',
        'probdist'   =>'Kroupa2001($m1)*$m1',
        'dphasevol'  =>'$dlnm1',
        );
}
elsif($duplicity == 1)
{
    # make a population of binary stars
    my $resolution = {
        m1 => 80,
        q =>50,
        P =>60
    };
    # my $resolution = {
    #     m1 => 20,
    #     q =>20,
    #     P =>20
    # };
    
    my $mmin =0.8;
    my $mmax = 20.0;
    $population->{_grid_options}{binary} = 1;
    # flat in ln(m)
    # $population->add_grid_variable
    #     (
    #     'name'       => 'lnm1', 
    #     'longname'   =>'Primary mass', 
    #     'range'      =>[log($mmin),log($mmax)],
    #     'resolution' => $resolution->{m1},
    #     'spacingfunc'=>"const(log($mmin),log($mmax),$resolution->{m1})",
    #     'precode'    =>'$m1=exp($lnm1);',
    #     'probdist'   =>'Kroupa2001($m1)*$m1',
    #     'dphasevol'  =>'$dlnm1',
    #     );

    # flat in m:
    $population->add_grid_variable
        (
        'name'       => 'm1', 
        'longname'   =>'Primary mass', 
        'range'      =>[$mmin,$mmax],
        'resolution' => $resolution->{m1},
        'spacingfunc'=>"const($mmin,$mmax,$resolution->{m1})",
        'probdist'   =>'Kroupa2001($m1)*$m1',
        'dphasevol'  =>'$dm1',
        );

#     # q=M1/M2 distribution flat in q between 0.1/M1 and 1.0
#     $population->add_grid_variable
#         (
#         'condition'  =>'$self->{_grid_options}{binary}==1',
#         'name'       =>'q',
#         'longname'   =>'Mass ratio',
#         'range'      =>['0.1/$m1',1.0],
#          'resolution'=>$resolution->{q},
#         'spacingfunc'=>"const(0.1/\$m1,1.0,$resolution->{q})",
#         'probdist'   =>"flatsections\(\$q,\[
# \{min=>0.1/\$m1,max=>0.8,height=>1.0\},
# \{min=>0.8,max=>1.0,height=>1.0\},
# \]\)",
#         precode     =>'$m2=$q*$m1;',
#         dphasevol   =>'$dq',
#         );
    # q=M1/M2 distribution flat in lnq between 0.1/M1 and 1.0
    $population->add_grid_variable
        (
        'condition'  =>'$self->{_grid_options}{binary}==1',
        'name'       =>'lnq',
        'longname'   =>'log Mass ratio',
        'range'      =>['log(0.1/$m1)',log(1.0)],
         'resolution'=>$resolution->{q},
        'spacingfunc'=>"const(log(0.1/\$m1),log(1.0),$resolution->{q})",
        'probdist'   =>"flatsections\(\$lnq,\[
\{min=>0.1/\$m1,max=>0.8,height=>1.0\},
\{min=>0.8,max=>1.0,height=>1.0\},
\]\)",
        precode     =>'$m2=exp($lnq)*$m1;',
        dphasevol   =>'$dlnq',
        );

     # orbital period Duquennoy and Mayor 1991 distribution
#      my $Prange = [-2.0,9.0];
#      $population->add_grid_variable
#          (
#           'name'       =>'logper',
#           'longname'   =>'log(Orbital_Period)',
#           'range'      =>$Prange,
#           'resolution' =>$resolution->{P},
#           'spacingfunc'=>"const($Prange->[0],$Prange->[1],$resolution->{P})",
#           'precode'=>"\$per = 10.0 ** \$logper;
# my \$eccentricity = 0.0;
# \$sep=calc_sep_from_period(\$m1,\$m2,\$per) if(defined \$m1 && defined \$m2);
# ",
#           'probdist'=>"gaussian(\$logper,4.8,2.3,$Prange->[0],$Prange->[1])",
#           'dphasevol'=>'$dln10per'
#          );
    # #flat_in_log_a
        my $seprange = ['log(3)','log(1e6)'];
        $population->add_grid_variable(
        # name of the variable
        'name'=>'lnsep',
        'longname'=>'ln(Orbital_Separation)',
        'range'=>$seprange,
        'resolution'=>$resolution->{P},
        'spacingfunc'=>"const($seprange->[0],$seprange->[1],$resolution->{P})",
        # precode has to calculation the period (required for binary_c)
        'precode'=>"my \$eccentricity = 0.0;my\$sep=exp(\$lnsep);my\$per=calc_period_from_sep(\$m1,\$m2,\$sep);",
        # flat in log-separation (dN/da~1/a) distribution (Opik law)
        'probdist'=>'const(log(3.0),log(1e6))',
        # phase volume contribution
        'dphasevol'=>'$dlnsep'
        );
}
else
{
    # Make a mixed population of single and binary stars
    # with equal sampling in equal times.
    my $dt = 100.0;
    $population->set(
        time_adaptive_mass_grid_step => $dt,
        );
    my $sampling_factor = 0.5; # deliberate oversampling (<0.5)
    my $time_adaptive_options = {
        max_evolution_time =>
            $population->{_bse_options}{max_evolution_time},
        stellar_lifetime_table_nm=>1000,
        time_adaptive_mass_grid_log10_time=>0,
        time_adaptive_mass_grid_step=>$dt,
        time_adaptive_mass_grid_nlow_mass_stars=>10,
        nthreads           =>$nthreads,
        thread_sleep       =>1,
        mmin               =>0.1,
        mmax               =>80.0,
        mass_grid_log10_time=>0,
        mass_grid_step      =>$dt*$sampling_factor,
        extra_flash_resolution=>0, # 1 = broken?
        mass_grid_nlow_mass_stars=>10,
        debugging_output_directory=>'/tmp',
        max_delta_m         =>1.0,
        savegrid            =>1,
        vb                  =>0,
        metallicity=>$population->{_bse_options}{metallicity},
        agbzoom             =>0,
    };
    distribution_functions::bastard_distribution(
        $population, {
            mmin         =>0.1,
            mmax         =>80.0,
            m2min        =>0.1,
            nm2          =>10,
            nper         =>10,
            qmin         =>0.0,
            qmax         =>1.0,
            necc         =>undef,
            useecc       =>undef,
            agbzoom      =>0,
            time_adaptive=>$time_adaptive_options,
        });
}

# link population to custom data parser function
$population->set(
    parse_bse_function_pointer => \&main::parse_data
    );

my %init = $population->initial_abundance_hash('Karakas2002',0.02);
my %isotope_hash = $population->isotope_hash();
my @isotope_list = $population->isotope_list();
my %nuclear_mass_hash = $population->nuclear_mass_hash();
my @nuclear_mass_list = $population->nuclear_mass_list();
my @sources = $population->source_list();
my @ensemble = $population->ensemble_list();

if(0){
print Data::Dumper->Dump([
    #\%init,
    #\%isotope_hash,
    #\@isotope_list,
    #\%nuclear_mass_hash,
    \@nuclear_mass_list,
    #\@sources,
    #\@ensemble
                         ]);
}

# uncomment this to show version information
#print $population->evcode_version_string();

# uncomment this to show the evcode's args list
#print join("\n",@{$population->evcode_args_list()});

# evolution the stellar population (this takes some time)


my $pdcall=0;
output($population);
$population->evolve();
# output log of setup

print STDERR $outdir;#print outdir

# done : exit
exit;

############################################################
# subroutines 
############################################################

sub parse_data
{
    my $initialline = $population->tbse_line();
    shift @$initialline;

    open (my $out, '>>',"$outdir/grid.txt");
    my $m1init= $initialline->[2];
    print {$out} $m1init." ";
    $m1init=~ s/\./_/g;
    my $m2init= $initialline->[3];
    print {$out} $m2init." ";
    $m2init=~ s/\./_/g;
    my $ainit= $initialline->[33];
    print {$out} $ainit." ";
    print {$out} "\n";
    $ainit=~ s/\./_/g;
    close $out;

    
    my $i=0;
    my ($population, $results) = @_;
    $pdcall=$pdcall+1;
    open(my $out, '>', "$outdir/M1_".$m1init."M2_".$m2init."a_".$ainit.".txt");
    close $out;
    #say "parsed_data called ".$pdcall." times";
    my $keep_me=0;#decides whether to keep the file or not.

    while(1)
    {
        # subsequent calls to tbse_line contain
        # (references to) arrays of data 
        my $la = $population->tbse_line();
        # first element is the "header" line
        my $header = shift @$la;

        # break out of the look if this is 'fin'
        last if ($header eq 'fin');

        # check if $header matches one of your
        # expected data lines, if so, act

        if($header eq 'KEMPSTARDATA')
        {
            # matched MY_STELLAR_DATA header
            #
            # ... so do stuff with $la.You have to
            # enable a line in binary_c which outputs
            # the mass, the probability and the timestep
            # (see log_every_timestep.c for examples)
            # and starts with MY_STELLAR_DATA.
            #
            my $modeltime = $la->[0];
            my $dt = $la->[1];
            my $mass1 = $la->[2];
            my $mass2 = $la->[3];
            my $coremass1 = $la->[4];
            my $coremass2 = $la->[5];
            my $stellartype1 = $la->[6];
            my $stellartype2 = $la->[7];
            my $radius1 = $la->[8];
            my $radius2 = $la->[9];
            my $roche_radius1 = $la->[10];
            my $roche_radius2 = $la->[11];
            my $luminosity1 = $la->[12];
            my $luminosity2 = $la->[13];
            my $teff1 = $la->[14];
            my $teff2 = $la->[15];
            my $mdot1 = $la->[16];
            my $mdot2 = $la->[17];
            my $num_thermal_pulses1 = $la->[18];
            my $num_thermal_pulses2 = $la->[19];
            my $dmH1 = $la->[20];
            my $dmH2 = $la->[21];
            my $dmHcrit1 = $la->[22];
            my $dmHcrit2 = $la->[23];
            my $dmHe1 = $la->[24];
            my $dmHe2 = $la->[25];
            my $dmHecrit1 = $la->[26];
            my $dmHecrit2 = $la->[27];
            my $nova_eta1 = $la->[28];
            my $nova_eta2 = $la->[29];
            my $num_novae1 = $la->[30];
            my $num_novae2 = $la->[31];
            my $WD_accretion_regime = $la->[32];
            my $sep = $la->[33];
            my $per = $la->[34];
            my $e = $la->[35];
            my $ce = $la->[36];
            my $interesting_flag = $la->[37];
            my $pre_nova_acc_rate = $la->[38];
            my $nova_overshoot_factor= $la->[39];
            my $dangmom= $la->[40];
            my $dangmom_postnova= $la->[41];
            my $dangmom_gw= $la->[42];
            my $dangmom_windloss= $la->[43];
            my $dangmom_windgain= $la->[44];
            my $dangmom_rlofloss= $la->[45];
            my $dangmom_rlofgain= $la->[46];
            my $dangmom_cbdisc= $la->[47];
            my $dangmom_tides= $la->[48];
            my $dangmom_nonconservativeloss= $la->[49];
            my $dangmom_nova= $la->[50];
            my $dangmom_artificial= $la->[51];
            my $angmom= $la->[52];

            
            if($stellartype1==15 && $stellartype2==15){
                #double degenerate Type 1a? flawed... but it doesn't spam.
                #single degenerate chanel is contained within 'interesting_flag'.
                $keep_me=1;
            }

            if($stellartype1>9 && $stellartype1 < 13 &&
             $stellartype2 >9 && $stellartype2 < 13){
                ###keep all double-WD sys
                #$keep_me=1
            }

            if($stellartype1==13 || $stellartype2==13){
                #keep all neutron stars
                $keep_me=1;
            }

            
            if($interesting_flag>0){
                # interesting flag raised, could be 
                # 1.Hnova, 2.Henova, 3.CO -> ONe conv,
                # 4.AIC->NS, 5.MCh Type Ia, 6.ELD, 
                # 7. double-WD#doesn't actually work, 8. Mchand_coalesence Type Ia 
                # 9. Violent_merger Type Ia, 10. HeCO Hybrid merger subluminous type Ia 
                # 11. HeCO Hybrid merger 'normal' Type Ia, 12. Subluminous HeCO hybrid AND violent type Ia
                # 13. Normal HeCO hybrid AND violent merger. 14. Merger induced ECSNe (occurs when merger product = ONeWD and > Mchand)
                $keep_me=1;
            }




            # bin mass to nearest 1.0 Msun
            #$mass = $population->rebin($mass, 1.0);
            #$results->{prob}->{$probability} += $probability;
            # add up the mass distribution in a histogram
            #$results->{mass_distribution}->{$mass} += $probability * $timestep;
            #$results->{mass}->{$mass1} = $mass1;

            #my @mass[i]=Dumper($mass);
            if ($i % 1 ==0)
            {
                #say $i;
                open (my $out, '>>',"$outdir/M1_".$m1init."M2_".$m2init."a_".$ainit.".txt");
                if($out)
                {
                    print {$out} $modeltime . " ";       #0
                    print {$out} $dt . " ";              #1                   
                    print {$out} $mass1 . " ";           #2         
                    print {$out} $mass2 . " ";           #3        
                    print {$out} $coremass1 . " ";       #4
                    print {$out} $coremass2 . " ";       #5
                    print {$out} $stellartype1 . " ";    #6
                    print {$out} $stellartype2 . " ";    #7
                    print {$out} $radius1 . " ";         #8
                    print {$out} $radius2 . " ";         #9
                    print {$out} $roche_radius1 . " ";   #10
                    print {$out} $roche_radius2 . " ";   #11
                    print {$out} $luminosity1 . " ";     #12
                    print {$out} $luminosity2 . " ";     #13
                    print {$out} $teff1 . " ";           #14
                    print {$out} $teff2 . " ";           #15
                    print {$out} $mdot1 . " ";           #16   
                    print {$out} $mdot2 . " ";           #17
                    print {$out} $num_thermal_pulses1 . " ";#18   
                    print {$out} $num_thermal_pulses2 . " ";#19
                    print {$out} $dmH1 . " ";            #20
                    print {$out} $dmH2 . " ";            #21
                    print {$out} $dmHcrit1 . " ";        #22
                    print {$out} $dmHcrit2 . " ";        #23
                    print {$out} $dmHe1 . " ";           #24
                    print {$out} $dmHe2 . " ";           #25
                    print {$out} $dmHecrit1 . " ";       #26
                    print {$out} $dmHecrit2 . " ";       #27
                    print {$out} $nova_eta1 . " ";       #28
                    print {$out} $nova_eta2 . " ";       #29
                    print {$out} $num_novae1 . " ";      #30
                    print {$out} $num_novae2 . " ";      #31
                    print {$out} $WD_accretion_regime . " ";#32
                    print {$out} $sep . " ";             #33
                    print {$out} $per . " ";             #34
                    print {$out} $e . " ";               #35
                    print {$out} $ce . " ";              #36
                    print {$out} $interesting_flag . " ";#37
                    print {$out} $pre_nova_acc_rate . " ";#38
                    print {$out} $nova_overshoot_factor . " ";#39
                    print {$out} $dangmom . " ";#40
                    print {$out} $dangmom_postnova . " ";#41
                    print {$out} $dangmom_gw . " ";#42
                    print {$out} $dangmom_windloss . " ";#43
                    print {$out} $dangmom_windgain . " ";#44
                    print {$out} $dangmom_rlofloss . " ";#45
                    print {$out} $dangmom_rlofgain . " ";#46
                    print {$out} $dangmom_cbdisc . " ";#47
                    print {$out} $dangmom_tides . " ";#48
                    print {$out} $dangmom_nonconservativeloss . " ";#49
                    print {$out} $dangmom_nova . " ";#50
                    print {$out} $dangmom_artificial . " ";#51
                    print {$out} $angmom . " ";#52
                    print {$out} "\n";
                    #$out->autoflush;
                }
                close $out;

            }
            $i=$i+1;
        }
    }
    if ($keep_me==0){
        unlink("$outdir/M1_".$m1init."M2_".$m2init."a_".$ainit.".txt");
    }
}

############################################################

sub output
{
    my $population = shift;
    print "OUTPUT to $outdir\n";


    open(my $out, '>', "$outdir/log");
    if($out)
    {
    print {$out} Dumper(\%{$population->{_grid_options}});
    print {$out} Dumper(\%{$population->{_bse_options}});
    }
    close $out;
}
