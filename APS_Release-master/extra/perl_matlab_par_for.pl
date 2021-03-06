#!/usr/bin/perl -w
## ------------------ Disclaimer  ------------------
# 
# BG Group plc or any of its respective subsidiaries, affiliates and 
# associated companies (or by any of their respective officers, employees 
# or agents) makes no representation or warranty, express or implied, in 
# respect to the quality, accuracy or usefulness of this repository. The code
# is this repository is supplied with the explicit understanding and 
# agreement of recipient that any action taken or expenditure made by 
# recipient based on its examination, evaluation, interpretation or use is 
# at its own risk and responsibility.
# 
# No representation or warranty, express or implied, is or will be made in 
# relation to the accuracy or completeness of the information in this 
# repository and no responsibility or liability is or will be accepted by 
# BG Group plc or any of its respective subsidiaries, affiliates and 
# associated companies (or by any of their respective officers, employees 
# or agents) in relation to it.
## ------------------ License  ------------------ 
# GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
## github
# https://github.com/AnalysePrestackSeismic/
##
# this is a perl script to run compiled matlab code in a parrallel for loop on a node
##
##############################################################################################
##  load required packages
## use eval to test the existance of modules
##
#### config is a core perl module, so assume it exists
use Config;
## check to see if perl version was built with threads
$Config{useithreads} or die('Recompile Perl with threads to run this program.');
##
##
$usagemessage = "perl_matlab_par_for.pl  <command to run> <as many args as you like, each is a seperate item to execute>";
##  If no arguments are supplied, prompt user with how to use script.###
if ($#ARGV < 0) {
    die($usagemessage);
}
## 
## set defaults ##############################################################################
$noofthreads = 4;
##
## check the command line for switches #######################################################
$para_commnda = $ARGV[0];
#print "a: ".$para_commnda."\n";
$para_commndb = $ARGV[1];
#print "b: ".$para_commndb."\n";
$para_commndc = $ARGV[2];
#print "c: ".$para_commndc."\n";
$para_commndd = $ARGV[3];
#print "d: ".$para_commndd."\n";
$para_commnde = $ARGV[4];
#print "e: ".$para_commnde."\n";
$para_commndf = $ARGV[5];
#print "f: ".$para_commndf."\n";
$para_commndg = $ARGV[6];
#print "g: ".$para_commndg."\n";
##
@fileargs = split(' ', $ARGV[7]);
##
$#item_to_loop = -1;
#print @fileargs ;
#print "total number of files :".$#fileargs."\n";
##
for ($argi=0; $argi <= $#fileargs; $argi++){
    $item_to_loop[$argi] = $fileargs[$argi];
    #print "filearg ".$argi." ".$fileargs[$argi]."\n";
}
print "last filearg ".$#fileargs." ".$fileargs[$#fileargs]."\n";
##
## fork manager is a cpan module, so need to test if it has been installed
## if not installed, install in a csh with
## setenv FTP_PASSIVE 1
## perl -MCPAN -e "install 'Parallel::ForkManager'"
## 
$usethread = 0;
##
eval("use Parallel::ForkManager");
if ($@ =~ m/Can\'t locate/) {
    #print "forkmanager package did not exist\n".$@."\n see in code of this program how to install\n";
} else {
    #print "package did exist\n";
    $usethread = 1;
}    
##
## test to see if threads enabled, if not then uses serial loop
$i = 0;
if ($usethread == 1){
    ##
    #print "in threaded loop now =================\n";
    ## start number of threads to use
    my $manager = new Parallel::ForkManager( $noofthreads );
    $#runcomms = -1;
    ##
    ##test the digital signatures of each file
    foreach (@item_to_loop) {
        print $item_to_loop[$i]."\n";
        $runcomms[$i] = $para_commnda." ".$para_commndb." ".$para_commndc." ".$para_commndd." ".$para_commnde." ".$para_commndf." ".$para_commndg." ".$item_to_loop[$i];
        #
        #print $runcomms[$i]."\n";
        $i++;
    }
    ###
    ##now run the commands
    umask(000);
    ##
    foreach my $command (@runcomms) {
      $manager->start and next;
      sleep 2;
      #print "running $command \n";
      system( $command );
      $manager->finish;
    }
    $manager->wait_all_children;
    #
    print "all threads completed\n";
    ##
} else {
    foreach (@item_to_loop) {
        #print $item_to_loop[$i]."\n";
        $commcj = $para_commnda." ".$para_commndb." ".$para_commndc." ".$para_commndd." ".$para_commnde." ".$para_commndf." ".$para_commndg." ".$item_to_loop[$i];
        @out = `$commcj`;
        chomp(@out);
        $i++;
    }
}
##
##
exit 1;


