###
# Duplicate the file from the input and print it to the output
###

awk '{ $2=$1; print; }' $1 > $2