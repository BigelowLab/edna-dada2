# edna-dada2

[Maine eDNA](https://umaine.edu/edna/) 

Divisive Amplicon Denoising Algorithm [dada2](https://benjjneb.github.io/dada2/index.html)

[charlie wiki](https://github.com/BigelowLab/charlie/wiki) (because Ben can't remember anything)


### Organization

 + `config` configuration files in yaml format
 
 + `docs` miscellaneous notes
 
 + `pbs` shell scripts to use when calling `qsub`
 
 + `renv` and `renv.lock` [renv](https://cran.r-project.org/package=renv) project management stuff 
 
 + `rscript` R scripts to manage workflows
 
 
### workflows

Whether as a PBS job, cron or manually workflow scripts accept at least one argument (configuration filename)

```
$ Rscript /path/to/Rscript/my_script.R /path/to/config/file.yaml
```

The config file will contain the identities of input and output paths.  The output path will be populated with subdirectories and a diagnostic and results files. Diagnostic files may include [software audits](docs/example-audit.txt) and [output logs](docs/example-log.txt).  Results files will vary by script.

