# Reinstalls the dadautils package in the users local R library
#
# usage: Rscript <script> [--] [--help] [--upgrade] [--opts OPTS] [--path PATH]
#        [--devel DEVEL]
# 
# Update dadautils
# 
# flags:
#   -h, --help     show this help message and exit
#   -u, --upgrade  set to upgrade dependicies
# 
# optional arguments:
#   -x, --opts     RDS file containing argument values
#   -p, --path     path to the package directory [default:
#                  /mnt/storage/data/edna/packages/dadautils]
#   -d, --devel    set to use devtools::load_all() instead of
                 devtools::install() for a temporary installation


library(magrittr)
library(argparser)

ARGS <- arg_parser("Update dadautils") %>%
  add_argument("--path",
               default = "/mnt/storage/data/edna/packages/dadautils",
               help = "path to the package directory") %>%
  add_argument("--upgrade", 
               help = "set to upgrade dependicies", 
               flag = TRUE) %>%
  add_argument("--devel",
               help = "set to use devtools::load_all() instead of devtools::install() for a temporary installation") %>%
  parse_args()   

  if (ARGS$devel){
    devtools::load_all(ARGS$path)
  } else {
    devtools::install(ARGS$path, upgrade = ARGS$upgrade)
}