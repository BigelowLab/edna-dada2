# edna-dada2

[Maine eDNA](https://umaine.edu/edna/) 

Divisive Amplicon Denoising Algorithm [dada2](https://benjjneb.github.io/dada2/index.html)

[charlie wiki](https://github.com/BigelowLab/charlie/wiki) (because Ben can't remember anything)


### Devel Queue

```
qsub -I -q devel -l walltime=8:00:00 -l ncpus=8,mem=8G  -N ben-edna
module use /mod/bigelow
module load dada2
```