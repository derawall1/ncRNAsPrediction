use strict;
use LWP::Simple qw( $ua );

# Make a request command (uncomment lines below if you want something else)
my $request = HTTP::Request->new( GET => 'https://probis.nih.gov/rest/align?structure_id1=4C0O.ABC&structure_id2=2CPE.A');
#my $request = HTTP::Request->new( GET => 'http://probis.cmm.ki.si/rest/align?structure_id1=1all.A&bsite1=CYC.175.A.7&structure_id2=3dbj.C&bsite2=CYC.201.C.7');
#my $request = HTTP::Request->new( GET => 'http://probis.cmm.ki.si/rest/align?structure_id1=1all.A&structure_id2=3dbj.C&return=pdb');
#my $request = HTTP::Request->new( GET => 'http://probis.cmm.ki.si/rest/scan?structure_id=5cyt.R');
#my $request = HTTP::Request->new( GET => 'http://probis.cmm.ki.si/rest/scan?structure_id=5cyt.R&bsite=HEM.105.R.5');
#my $request = HTTP::Request->new( GET => 'http://probis.cmm.ki.si/rest/scan?structure_id=5cyt.R&z_score=2.0&return=json');

# Decide about the content type you want to get in return (default is XML) (applies to get_alignments and get_representative; other two commands return "text/plain")
$request->header(Accept => "application/json");
#$request->content_type( 'application/xml' );

# Send the HTTP request
my $response = $ua->request( $request );

# Check to see if there is an error
unless( $response->is_success ) {
print "\n Error: ", $response->status_line, "\n";
}

# Output response
print "ProBiS returned:\n", $response->content;

konc@cmm.ki.si