# set variables for alignsurface.sh script

subjnum=$1

subjprefix=s

if [ "$subjnum" = "" ]; then
  echo "ERROR: No subject id given"
  exit 1
fi

#subjnums="05 03 04 02 06 07 08" # now set in the script itself
subjid=${subjprefix}${subjnum}

cd ..
rootdir=`pwd`"/"

if [ `echo $rootdir | wc -w ` -ne 1 ]; then
	echo "rootdir ${rootdir} contains spaces - cannot continue"
	exit 1
elif [ ! -e $rootdir ]; then
	echo "rootdir ${rootdir} does not exist - cannot continue"
	exit 1
fi

scriptdir=$rootdir/scripts/
matlabdir=$scriptdir/matlab/

fsrootdir=${rootdir}/fs/
glmdir=${rootdir}/glm+orig/${subjid}/

fsdir=${rootdir}/fs/${subjid}/
surfdir=${fsdir}/${subjid}/surf/
sumadir=${surfdir}/SUMA/
aligndir=${fsdir}/align/
groupfsdir=${rootdir}/fs/groupana/

export AFNI_SLICE_SPACING_IS_GAP=N
export AFNI_DICOM_RESCALE=YES
export FLOATIZE=YES

prefix=${subjid}_
subjcount=`echo $subjnums | wc -w`

TMP=__tmp_

runmatlab () {
	local X="${1};exit"	
	echo "Running matlab with $X"
	
	
	
	pths="matlab /Applications/MATLAB_R2008b.app/bin/matlab /Applications/MATLAB_R2008a/bin/matlab"

	for pth in $pths; do
		which $pth || continue
		${pth} -nojvm -nodisplay -nosplash -r "$X" || exit 1
		return
	done
	exit 1
}

runallsubjects () {
	shift
	for k in $subjnums; do
		cd ${scriptdir}
		./$0 $k $@ || exit 1
	done
}

removepath () {
	local s="$1"
	local length=$((`echo $s | wc -c`-1))
	
	local idx="$length"
	
	while [ $idx -gt 0 ]; do
	  local thechar=`echo $s | cut -c $idx`

	  if [ $thechar = "/" ]; then
	    break;
	  fi
	  local idx=$(($idx-1))
	done
	   
	local result=`echo $s | cut -c $((${idx}+1))-${length}`
	echo $result
}

linkheadbrik() { # to link a file foo/bar+orig.{HEAD,BRIK} to ./bar+orig.{HEAD,BRIK}
	local s="$1"
	local t="${2}/"`removepath "$1"`
	
	echo "linking $s to $t"
	
	for ext in HEAD BRIK; do
	  if [ ! -e "${s}.${ext}" ] ; then
	    echo "File ${s}.${ext} does not exist!"
	    exit 1 
	  fi 
	  if [ -e "${t}.${ext}" ]; then 
		  rm "${t}.${ext}"
	  fi
	  ln "${s}.${ext}" "${t}.${ext}" || ln "${s}.${ext}.gz" "${t}.${ext}.gz" || cp "${s}.${ext}" "${t}.${ext}" || cp "${s}.${ext}.gz" "${t}.${ext}.gz" || exit 1 #if linking is not supported, then copy the file
	done
}

calc() { #simple computation, e.g. use k=`calc "${i}+3"`
	local s="${1}"
	echo "$1" | bc
}

getfield() { #getfield IDX STR takes the IDX-th field in STR, where fields are seperated by spaces
	if [ $# -lt 2 ]; then
		echo "Error: getfield needs 2 arguments"
	else

		local idx=`calc "$1 +1"`
		echo $@ | cut -f${idx} -d' '
	fi
}

getidxs() {
	local s=`echo $@ | wc -w`
	count -digits 1 1 $s
}

infomsg() {
	echo -e "\n\n *** $@ *** \n\n"
}
