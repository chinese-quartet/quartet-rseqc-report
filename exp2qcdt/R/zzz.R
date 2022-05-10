.onAttach <- function(lib, pkg, ...){
   packageStartupMessage(welcomeMessage())
}

welcomeMessage <- function(){
   paste0("\n",
         "---------------------\n",
         "Welcome to exp2qcdt version ", utils::packageDescription("exp2qcdt")$Version, "\n",
         # "\n",
         "Type ?exp2qcdt for the main documentation.\n",
         "The github page is: https://github.com/chinese-quartet/quartet-rseqc-report/tree/master/exp2qcdt\n",
         "\n",
         "Need help getting started? Github README: https://github.com/chinese-quartet/quartet-rseqc-report/tree/master/exp2qcdt \n",
         "\n",
         "To suppress this message use:  ", "suppressPackageStartupMessages(library(exp2qcdt))\n",
         "---------------------\n"
   )
}
