package Socket::socket_helper;

use strict;
use IO::Socket::INET; 
use JSON::PP;
use Time::HiRes qw(sleep);

# we need to export the send_to_compadre_helper function for it to be used in other scripts
use Exporter qw(import);
our @EXPORT_OK = qw(send_to_compadre_helper);
# To make the error handling logic easier to parse, we are going to create 
# a separate function. This function can be expanded to make it more clear 
# what errors occurred
sub check_errors {    

    my ($error_msg) = @_;

    if ($error_msg =~ /^ERROR:\s*(.*)$/s) {
        my $error_details = $1;
        die "Python COMPADRE helper error: $error_details\n";
    }
    elsif ($error_msg =~ /^(ERSA_ERROR|PADRE_ERROR|SOCKET_ERROR|POP_CLASSIFIER_ERROR):\s*(.*)$/s) {
        my $error_type = $1;
        my $error_details = $2;
        die "Python $error_type: $error_details\n";
    }
}

sub send_to_compadre_helper {
    my ($data, $port) = @_;
    $port //= 6000;  
    my $host = $ENV{COMPADRE_HOST} // 'localhost';


    # We are going to implement an exponential backoff retry strategy for 
    # the socket just incase there are weird network issues. The 
    # IO::Socket::INET->new returns undef if it timeouts waiting for the 
    # connect. We can use this to check if we have successfully connected 
    # to the socket or if we can't connect to the socket for some reason
    my $socket;
    my $max_attempts = 5;
    my $retry_wait = 0.1;


    for (my $attempt = 1; $attempt <= $max_attempts; $attempt++) {
      $socket = IO::Socket::INET->new(
          PeerAddr => $host,
          PeerPort => $port,
          Proto    => 'tcp',
          Timeout  => 5,  # Increase timeout for debugging
      );

      # This checks if we have successfully connected to the socket 
      # and then exits if we have
      last if $socket;
      
      # if we make it here then we have going through all of the retry 
      # attempts and still haven't connected to the server so we should 
      # terminate the program
      if ($attempt == $max_attempts) {
        die "failed to connect to the Python server at $host:$port: $!\n. Exiting program now...";
      }

      # Now that we have checked if we connected successfully and checked if 
      # we have reached the number of retries we can perform the exponential  
      # delay where we sleep for the set time and then multiple the sleep time 
      # to grow the delay
      sleep($retry_wait);
      $retry_wait *= 2;
    }

    # Enable keep-alive on the socket
    $socket->sockopt(SO_KEEPALIVE, 1) if $socket->can('sockopt');
    $socket->blocking(1);


    my $bytes_sent = $socket->print($data);
    
    if (!defined $bytes_sent) {
        close($socket);
        die "Failed to send data to Python server: $!\n";
    }
    
    $socket->flush();
    
    # Read response with error checking
    my $response_json = <$socket>;
    my $response; # This value represents the parsed output from the python server
    if (defined $response_json) {
      # we now need to parse the response since it is a json object. Use block eval as a catch incase the JSON parsing fails for some reason.
      my $parsed_response = eval {decode_json($response_json) };
      
      # make sure the json parsing didn't fail but return a message if it 
      # did and close the socket
      if ($@) {
        close($socket);
        die "Failed to parse the JSON response from the Python server: $@\nResponse was: $response_json\n";
      }

      # we can now check the error key of the message
      if ($parsed_response->{status} eq "error") { 
        # Make sure the socket gets closed on error
        close($socket);
        # Check if response indicates a broad error where the message starts 
        # with    
        my $error_msg = $parsed_response->{message};
        # We are going to check the different errors and then kill the program
        check_errors($error_msg);
      }
      elsif ($parsed_response->{status} eq "success") {
        # here we need to gather the response
        $response = $parsed_response->{result};
      } 
    }
      # Technically there is a case here where the status != "error" or 
      # "success" but those are the only values I coded in the response JSON
    else {
        $response = "No response";
        warn "No response received from Python server\n";
    }
    
    close($socket);
    return $response;
}
# Return 1 to indicate the module was successfully loaded.
1;
