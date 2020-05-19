# Configuration Files

`eDNA-dada2` scripts are designed to be called with a single argument, the path to a configuration file.  The configuration file is a convenient way to organize and record input parameters to a process that can become part of the meta data for a reproducible project.

A [configuration file](https://github.com/BigelowLab/edna-dada2/blob/master/config/dada2_example.yml) is a text file that is specially formatted to be machine readable. While there are a variety of formats to chose among, we have chosen [yaml](https://learnxinyminutes.com/docs/yaml) as our preferred format.  

#### Required elements

Required elements include a version id (please - no spaces!), and email address, and input and output path descriptions.  If one or more of the required elements are missing then the script will stop and report an error.  Obviously you will want to tailor these for your processing run.  

```
version: dada2_example
email: "workerbee@beehive.org"
input_path: "/home/btupper/edna/data/examples/dada2/MiSeq_SOP"
output_path: "/home/btupper/edna/data/examples/dada2/MiSeq_SOP_results"
```

A useful habit is to name the configuration file with the name of the version run, in this case name it `dada2_example.yaml`.  If you are stumped about a the best version naming system, try starting with something like `v00.000`, `v00.001`, `v00.002`, ...  It looks boring and mundane, but this approach allows you up to 100 major versions `v00-v99`, each with `000-999` minor versions.  You can also replace the major version with something project or customer minded, like `Penobscot_2019.000` or `Jones2020.000`.  Of course, you can adapt to fit your own needs but try to be systematic from the get go. 


#### Global elements

Global elements define attributes that are defined for the entire R session.  Each has a default value as shown below so you can skip these, but best practice is to include them.

```
verbose: info
multithread: auto
```

`verbose` is used by the messaging/logging system.  Leaving it at `info` is probably the best choice unless you are debugging.  Options include 'trace', 'debug', 'info' (default), 'warn', 'error', and 'fatal'.

`multithread` can be either a number of cpus (cores), the word 'auto' or a boolean 'yes' or 'no'.   If you are running you process in a PBS queue then your best choice is 'auto', in which case the software will select the correct number for you. (Tech note, the software will grab the environmental value of `$NCPUS`.) If you are **not** operating within a PBS queue then you should set this to the correct number, such as 4, 28 or whatever you want to allocate. 

#### Individual step elements

After the required and global elements, the configuration file can be roughly organized around the steps you take in your script.  These elements can be fairly organic in designed, and are meant to facilitate setting parameters/arguments for individual steps.  The example below shows a section of elements for executing the `filterAndTrim()` function from the `dada2` R package. The name is arbitrary, but again employing a systematic naming system can be useful.

```
dada2_filterAndTrim_filtN:
  name: "filtN"
  truncLen: 
  - 275
  - 225
  maxN: 0.0
  maxEE:
  - 2
  - 2
  truncQ: 2
  rm.phix: yes
  compress: yes
```
 
Within the script that calls `dada2::filterAndTrim()` the script will harvest the variables defined in the config file.  The script will translate the call like this...

```
results <- dada2::filterAndTrim(inputs, truncLen = c(275,225), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE)
```

There is some customization to the script required if you want to add/remove parameters to the call - in other words, the configuration system helps but it isn't smart.


### Other Familiar Elements

This is a simple run-down of some of the elements we have used so far.  Not all are needed of required, but they are documented here to provide some ideas on how we have used configurations to date. 

This can be used to set the number of `dada2::plotQualityProfiles()` to add to the output PDF.

```
dada2_plotQualityProfiles:
  nplots: 2
```

The `dada2::addSpecies()` requires a reference fasta which we can specify here.

```
dada2_addSpecies:
  refFasta: "/mnt/storage/labs/countway_nfs/dada2_data/Highland/Cyano_16S/cutadapt/silva_species_assignment_v138_dada2.fasta"
```

The `datautils::taxa_remove()` accepts a named list of items to remove by variable form taxa outputs.  In the example below, we would filter out rows with 'Chloroplast' in the `Order` column, and filter out rows with either 'foo' or 'bar' in the `Phylum` column.

```
taxa_remove:
  Order:
    - "Chloroplast"
  Phylum:
    - "foo"
    - "bar"
```

Provide references for finding/counting primers.

```
primer:
  FWD: "CYGCGGTAATTCCAGCTC"
  REV: "AYGGTATCTRATCRTCTTYG"
```