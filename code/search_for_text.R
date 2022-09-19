# look for a string in all .R files in a folder

library(dplyr)

folder <-  'code'
{
  files <-list.files(folder, recursive=T, pattern = '*.R') # get files from folder
  # keep those that end on .R

  search_for <-  "calc_LL_for_selection"

for(file_n in files){
  # open file
  zz <- file(paste(folder,file_n, sep='/'), "rt")  # open an output file connection
  # read all lines of the file
  my_text <- readLines(zz) 
  # close connection
  close(zz)
  # my_text[[1]]
  grep(search_for, my_text) -> grpd # grpd : grepped : grEPpEd
  if(length(grpd)>0) {print(file_n)
    print(grpd)}
}
}
