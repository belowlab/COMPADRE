package LogConfig::configure;
use strict;
use File::Spec;
use Log::Log4perl;
use Exporter 'import';

our @EXPORT_OK = qw(configure_logger get_logger_level);

# We need to convert the verbosity integer to a loglevel value that the 
# logger uses
sub get_logger_level {
  my ($verbosity) = @_;
  # Log level such as WARN, DEBUG, INFO
  my $loglevel;

  if ($verbosity == 0) {
    $loglevel = "WARN";
  }
  elsif ($verbosity == 1) {
    $loglevel = "INFO";
  }
  elsif ($verbosity == 2) {
    $loglevel = "DEBUG";
  }
  elsif ($verbosity < 2) {
    $loglevel = "TRACE";
  }
  else {
    die "Unknown loglevel $verbosity";
  }
  return($loglevel);
}

sub configure_logger {
  my ($logger_filepath, $verbosity) = @_;

  Log::Log4perl::Logger::create_custom_level(
      "PROGINFO",       
      "FATAL",
      10000
  );

  # set the configuration for the logger
  my $conf = qq(
      log4perl.rootLogger              = $verbosity, Screen, File
      
      log4perl.appender.Screen         = Log::Log4perl::Appender::Screen
      log4perl.appender.Screen.stderr  = 0
      log4perl.appender.Screen.layout  = Log::Log4perl::Layout::PatternLayout
      log4perl.appender.Screen.layout.ConversionPattern = %m%n
      
      log4perl.appender.File           = Log::Log4perl::Appender::File
      log4perl.appender.File.filename  = $logger_filepath
      log4perl.appender.File.mode      = append
      log4perl.appender.File.layout    = Log::Log4perl::Layout::PatternLayout
      log4perl.appender.File.layout.ConversionPattern = %d %p %m%n
  );

  # Initialize Log4perl
  # We use 'init_once' so it doesn't crash if you accidentally load it twice
  Log::Log4perl->init_once(\$conf);

  # Redirect standard warn and die signals to the logger
  $SIG{__WARN__} = sub {
      local $Log::Log4perl::caller_depth = $Log::Log4perl::caller_depth + 1;
      Log::Log4perl::get_logger()->warn(@_);
  };
  
  $SIG{__DIE__} = sub {
      if($^S) {
          # We're in an eval {} and don't want log this message but catch it later
          return;
      }
      local $Log::Log4perl::caller_depth = $Log::Log4perl::caller_depth + 1;
      Log::Log4perl::get_logger()->fatal(@_);
      die @_; # Now die really
  };
}

1;
