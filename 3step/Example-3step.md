---
title: "An example ASV 3-step workflow for 18S"
output: github_document
---



The [dadautils](https://github.com/BigelowLab/edna-dada2) and [charlier](https://github.com/BigelowLab/charlier) 
packages have been designed to permit the user to mix and match steps. Every analysis requires some
level of customization, and this document shows how you might customize the automated workflow
into a semi-automated 3-step workflow: preprocess, user-supervision and postprocess. The need to have the
use provide supervision can arise when working with a novel dataset, a noisy dataset or when it is desirable to
have a "first-look" at data quality before investing forther in the the prusuit of taxonomic classification.

This example shows how you might choose to break your analysis into the three parts.  Keep in mind, this is not 
prescriptive - you have many different options.

### Preprocess

The user provides a [configuration file (YAML format)](https://github.com/BigelowLab/edna-dada2/blob/master/3step/asv_18S_preprocess.yaml) 
and an [R script](https://github.com/BigelowLab/edna-dada2/blob/master/3step/asv_18S_preprocess.R) to a PBS session.

Preprocessing steps might include ...

+ verifying the files are complete (and if forward/reverse pairs),
+ running cutadapt to remove primers,
+ assessing the overlap between forward and reverse reads, 
+ assessing quality of the reads, and 
+ running dada2's filterAndTrim() and learnErrors() functions,
+ recommending read cutoff values for the user to apply to dada filtering.

The configuration file can contain information on what steps to run (or not!) To be a part of your metadata, 
the config should use the `stage` element to clearly indicate that it is 'preprocess'.  At the
end the script should write a copy of this config with the value set to something that indicates
something after preprocessing such as "user-supervision" or "step-2".


```r
CFG <- charlier::read_config("asv_18S_preprocess.yaml")
CFG$stage
```

```
## [1] "preprocess"
```

If you examine the config you will notice that includes metadata for many more workflow steps that you might
in a preprocessing step.  That is OK! Think of the configuration as a menu; every option is listed but
you don't need to order every item. The script you design will cherry pick just the stuff you need.  Note 
that at the end of the proeprocess script we save the original "name-preprocess.yaml" and a copy 
"name-postprocess.yaml" in the output folder. That way you retain the metadata along the way - what you did 
in the preprocess step and what you did in the postprocess step even though they may contain mostly-but-not-completely 
the same info. Super clear! 

### User supervision

At this step the user takes time to review the outputs of the preprocessing step.  If all looks satisfactory
then the user updates the configuration file to reflect any modifications and that the dataset is prepared for 
the post processing step.  You can edit the configuration file to suit your own needs.  One particular change may be
the change from `CFG$dada2_filterAndTrim$truncLen` which started as "auto" in preprocessing but may now hold the path to the 
CSV formatted file that holds the per-file-pair recommended `truncLen` values. The postprocessing script and
`dadautils::filter_and_trim()` function will adapt accordingly.

### Postprocess

The user provides an [updated configuration file (YAML format)](https://github.com/BigelowLab/edna-dada2/blob/master/3step/asv_18S_postprocess.yaml) 
and an [R script](https://github.com/BigelowLab/edna-dada2/blob/master/3step/asv_18S_postprocess.R) to a PBS session.  Don't forget, you should
flag the `stage` in the config as 'postprocess'.


```r
CFG <- charlier::read_config("asv_18S_postprocess.yaml")
CFG$stage
```

```
## [1] "postprocess"
```

Postprocessing steps might include...

+ running dada(),
+ merging pairs,
+ assigning taxonomy,
+ writing diagnostic files.

... but what you actually run is up to you.  

