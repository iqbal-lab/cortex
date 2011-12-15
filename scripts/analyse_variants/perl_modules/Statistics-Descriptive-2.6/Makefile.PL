use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'		=> 'Statistics::Descriptive',
    'VERSION_FROM'	=> 'Descriptive.pm', # finds $VERSION
    ($] >= 5.005 ?
      ('ABSTRACT'   	=> 'Supplies statistical methods for perl5',
       'AUTHOR'     	=> 'Colin Kuskie (ckuskie@cadence.com)') : ()),
    'dist'		=> { COMPRESS => "gzip -9f" }
);
