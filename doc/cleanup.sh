for i in `ls -1 main.* | grep --invert-match tex$`; do rm $i; done
rm -rf helpfiles
rm surfing
