# stew           <img src="https://github.com/advaitb/stew/blob/main/stew.gif" width="70" height="70">
Create a nice stew of your reads!


**stew** uses HyperLogLog arrays for diversity sampling.


### Installation:

* Source:    
	```
	git clone https://github.com/advaitb/stew.git
	cd stew    
	cmake -DCMAKE_BUILD_TYPE=Release  
	make 
	```   

### Parameters:
```
Usage: stew [Subcommand] [options] [input.*|input1.*|input2.*] [out.*...]
Subcommands:
	S - Single end read mode
	P - Paired end read mode

Main options:
	-t (--threads) - Number of threads [Default: 1]
	-p (--platters) - Number of platters (arrays) of HLL structures [Default: 10, Max: 50]
	-c (--cups) - Number of cups (bits) in each HLL platter (array) [Default: 8, Min: 4, Max: 16]
	-k (--kmers) - Kmer size [Default: 23, Max: 100]
	-x (--select) - Selectivity for similarity [Default: 0.5, Min: 0 (least selective), Max: 1 (most selective)]
	-m (--momentum) - Momentum applied to boost score (Useful in bigger datasets) [Default: 0.000001, Max 0.001]
	-h (--help) - Print usage
	-v (--version) - Print version

Subcommand postitional options:
	S
		S [input.*] [out.*]
	P
		P [input1.*] [input2.*] [out1.*] [out2.*]
```
 
