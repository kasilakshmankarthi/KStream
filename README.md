# KStream
Modified Streamv1 and Streamv2 to suit QDT model runs

Added support to run individual subtests

### Streamv1 benchmark
Usage:

./binaries/KStreamv1.\<arraysize\>.\<platform\>.\<compiler\>.elf [-T \<copy/scale/add/triad/all\>] \
                    [-R \<repetitions\>] \
                    [-A \<enable(1)/disable angel signals(0)\>]

KStreamv1 subtests strings: </br>
 "all" </br>
  "copy" </br>
  "scale" </br>
  "add" </br>
  "triad" </br>

#####Example:
To run all subtests in KStreamv1 (for 16K array with 10 repeats) </br>
  ./binaries/KStreamv1.16K.x86_64.gcc.6.2.0.linux.elf -T all -R 10 -A 0

To run copy subtest in KStreamv2 (for 16K array with 10 repeats) </br>
  ./binaries/KStreamv1.16K.x86_64.gcc.6.2.0.linux.elf -T copy -R 10 -A 0

### Streamv2 benchmark
Usage:

./binaries/KStreamv2.\<arraysize\>.\<platform\>.\<compiler\>.elf [-T \<fill/copy/daxpy/sum/all\>]  [-M \<size(K/M)\>]\
                [-P \<NPAD (offset. Match R)\>] [-R \<repetitions\>] \
                [-A \<enable(1)/disable angel signals(0)\>] \n

KStreamv2 subtests strings </br>
  "all" </br>
  "fill" </br>
  "copy" </br>
  "daxpy" </br>
  "sum" </br>
  "fillzero" </br>

Currently M and P are suported as command line arguments to use static arrays and not malloc

#####Example:
To run all subtests in KStreamv2 (for 16K array with 10 repeats) </br>
  ./binaries/KStreamv2.16K.x86_64.gcc.6.2.0.linux.elf -T all -R 10 -A 0

To run copy subtest in KStreamv2 (for 16K array with 10 repeats) </br>
  ./binaries/KStreamv2.16K.x86_64.gcc.6.2.0.linux.elf -T copy -R 10 -A 0


