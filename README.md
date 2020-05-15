# edna-dada2

[Maine eDNA](https://umaine.edu/edna/) 

Divisive Amplicon Denoising Algorithm [dada2](https://benjjneb.github.io/dada2/index.html)

[charlie wiki](https://github.com/BigelowLab/charlie/wiki)

### Requirements

 + Install into your home R library 2 packages from github
    
    - [compressr](https://github.com/BigelowLab/compressr)
    
    - [dadautils](https://github.com/BigelowLab/dadautils)
    
 + Load the dada2 module
 
``` 
module use /mod/bigelow
module load dada2
```

### Organization

 + `config` configuration files in yaml format
 
 + `docs` miscellaneous notes
 
 + `pbs` shell scripts to use when calling `qsub`
 
 + `rscript` R scripts to manage workflows
 
 + `scratch` a place to put files that will be ignored by git
 
 
### Workflows

Whether as a PBS job, cron or manually workflow scripts accept at least one argument (configuration filename)

```
$ Rscript /path/to/Rscript/my_script.R /path/to/config/file.yaml
```

The [config file](config/dada2_example.yml) will contain among other things the identities of input and output paths.  The output path will be populated with subdirectories and a diagnostic and results files. Diagnostic files may include [software audits](docs/example-audit.txt) and [output logs](docs/example-log.txt).  Results files will vary by script.

#### Configuration file

A [configuration file](config/dada2_example.yml) is a text file that is specially formatted to be machine readable. While there are a variety of formats to chose among, we have chosen [yaml](https://learnxinyminutes.com/docs/yaml) as our preferred format.  
See the [wiki page](https://github.com/BigelowLab/edna-dada2/wiki/Configuration-Files) for details.


