# to use this package you must have and environment variable named MWAS_GWAS_DIR
# pointing to the top-level git repo folder.
# you can add the line "export MWAS_GWAS_DIR=/path/to/mwasgwas" to the file
# ".bash_profile" in your home directory

# required libraries in R
repos = 'http://cran.mtu.edu'
#install.packages('infotheo',repos=repos)
install.packages('randomForest',repos=repos)
install.packages('vegan',repos=repos)
install.packages('optparse',repos=repos)
install.packages('MASS',repos=repos)
install.packages('car',repos=repos)
install.packages('edgeR',repos=repos)
