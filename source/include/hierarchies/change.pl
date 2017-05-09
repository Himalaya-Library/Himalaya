#!/usr/bin/perl -w
my %upcut = (
    "Dmsqst1" => "3",
    "Dmst12" => "3",
    "Mgl" => "4",
    "Mst1" => "3"
    );
#@indices = ("1Full.inc","12Full.inc","2Full.inc");
#foreach my $index (@indices){
#$fileName = "sigS$index";
#open(my $file, "<", "h9q2/$fileName");
open(my $file, "test.t");
my $sig = do{
    local $/;
    <$file>;
};
close $file;
@denoms = ();
foreach my $key(keys %upcut){
    $cut = $upcut{$key};
    
    for(my $i = $cut + 1; $i <= 10; $i++){
	if($i < 10){
	    $sig =~ s/pow$i\($key\)/upcut\(pow$i\($key\),cut\)/g;
	}
	else{
            $sig =~ s/pow$i\($key\)/upcut\(power$i\($key\),cut\)/g; 
	}
    }
}
print "$sig\n";
#open(my $file2, ">", $fileName);
#print $file2 $sig;
#close($file2);
#}

    

