#!/bin/bash
#
# this script aligns FreeSurfer (FS) surfaces with anatomical and 
# functional data, for viewing in AFNI. It requires that FreeSurfer 
# and AFNI are set up, and that AFNI +orig files are used.
#
# The only purpose of this script is to ensure that surfaces, 
# anatomical and functional volume data are all aligned, and that 
# these can be viewed easily using AFNI/SUMA
# 
# * Syntax: ./alignsurface.sh SUBJECT ACTION [ICOLD]
# 
#  where 
#  - SUBJECT is the subject id (e.g. "03"), "all" for all subjects.
#  - ACTION is either all (for individual subjects, actions #1 below), 
#    groupall (for group maps, to be done after each subject, actions #2),
#    ico (assuming a full pass with onother ICOLD has been run before),
#    or one single action as listed below
#  - ICOLD is the number of linear divisions for MapIcosehedron, default 100
#
# * List of possible actions
# -) recon: freesurfer surface reconstruction
# 1) tosuma: converts surfaces to AFNI (ASCII) format
# 1) align: alignment between aligned and freesurfer anatomical
# 1) mapico: resamples surface meshes (see "ld" variable below)
# 1) fixhead: fixes the header of the aligned anatomical
# 1) apply: applies alignment parameters
# 1) tlrc: finds alignment to talairach space
# 1) applytlrcvol: applies alignment for anatomical to talairach space
# 1) applytlrcsurf: applies alignment for surfaces to talairach space
# 1) avgsurf: averages pial and white matter surface
# 1) makespec: constructs spec file for viewing in afni
# 1) linkup: links files back to root subject freesurfer directory
# 2) groupavgsurf: compute average of surfaces across subjects
# 2) groupmakespec: constructs group spec file for viewing in afni
#
# Quick start: 
# - set directory structure so that rootdir/s{subjnum} contains
#   s{subjnum}_anatomical+orig (original anatomical) and
#   s{subjnum}_anatomical_al+orig (anatomical aligned to EPI). The 
#   former is used for freesurfer reconstruction, the latter for 
#   alignment of the surfaces to the functional volumes
# - run freesurfer reconstruction: ./alignsurface.sh s{subjnum} recon
#   this will take about a day.
# - run individual subject alignments etc: ./alignsurface.sh s{subjnum} all
# - after processing of all subjects, make group surface averages:
#   ./alignsurface.sh foo groupall
#
# Misc notes
# - I really should learn Python, yes I know.
# - The script may contain bugs etc, use at your own risk and don't sue me
# - Certain assumptions are made of folder locations (see setvars.sh)
# - Ensure your file and directory names do not contain spaces
# - To use a nifti file foo, use 3dcopy foo+orig foo.nii 
# - Note that AFNI is not a fan of oblique datasets. Consider to nuke
#   the transformation matrix with 3drefit -deoblique
#
# NNO Sept 2009, updated June 2010 <n.oosterhof@bangor.ac.uk>

. setvars.sh || exit 1

if [ "$3" = "" ]; then
	ld=100 # number of linear divisions in resampling meshes - see MapIcosahedron
else
	ld="$3"
fi

action="$2"
if [ "$action" = "" ]; then
	action="all"
fi

if [ "${subjnum}" = "all" ]; then
	echo "Running all subjects"
	runallsubjects $@
	exit 0
fi

anatfn=${subjid}_anatomical_fb  #anatomical used for freesurfer recon-all
anatalfn=${subjid}_anatomical_al_fb # assumed to be aligned with functional volume
subjnums="03" # all subject numbers (check subjprefix in setvars.sh)
  
case $action in 
    recon)		
		mkdir $fsdir
		cd $fsdir
		
		export SUBJECTS_DIR=`pwd`
				
		# make sure all paths are set for fs
		source ~/.bash_profile
		
		# convert nifti to fs-mgh format
		rm ${anatfn}.mgh
		3dcopy -overwrite ${anatfn}+orig ${anatfn}.nii || exit 1# remove this line if you use nifti rather than AFNI files
		mri_convert ${anatfn}.nii ${anatfn}.mgh || exit 1
		
		# delete previous freesurfer folder, if it exists
		rm -rf ${subjid}/
		
		# start reconstruction
		recon-all -subject ${subjid} -i ${anatfn}.mgh -all -cw256
	;;
	
	mrimakesurf)
		cd $fsdir
		
		export SUBJECTS_DIR=`pwd`
				
		# make sure all paths are set for fs
		source ~/.bash_profile

		mris_make_surfaces -noaparc -mgz -T1 brain.finalsurfs s03 lh
	;;
	
	tosuma)
	# converts surface files from reconstruction to AFNI SUMA ASCII
		if [ ! -e $surfdir ]; then
		  echo "Directory $surfdir does not exist!"
		  exit 1
		fi
		
		cd $surfdir || exit 1
		rm -rf $sumadir
		
		echo -e "\n\n*** Processing subject $subjid for ld=$ld ***\n\n"
		\@SUMA_Make_Spec_FS -sid ${subjid} || exit 1
		echo "done"
		cd $sumadir
	;;
	
	align)
		# finds the affine transformation from oligned anat to reconstructed fs anat
		rm -rf ${aligndir}
		mkdir ${aligndir}
		cd ${aligndir}
		
		fsvol="${subjid}_SurfVol+orig" #volume from freesurfer
		linkheadbrik ${sumadir}/${fsvol}	.
		linkheadbrik ${fsdir}/${anatalfn}+orig .
		
		matfileinv=${subjid}_SurfVol_al_mat_inv.aff12.1D
		matfile=${subjid}_SurfVol_al_mat.aff12.1D
		
		3dAllineate -base ${anatalfn}+orig -source ${fsvol} -1Dmatrix_save $matfileinv -cmass
		cat_matvec -ONELINE ${matfileinv} -I > ${matfile} || exit 1
		matfile3x4="3x4_${matfile}"
		cat_matvec ${matfile}  > ${matfile3x4} || exit 1
		fsvolal="${subjid}_SurfVol_al+orig"
		3dwarp -matvec_in2out $matfile3x4 -prefix $fsvolal $fsvol || exit 1
	;;
	
	mapico)
	# runs mapicosedron, so that group analysis can be done without interpolation
		cd $sumadir
		echo "Mapping to standard space, ld=${ld}"
		for hemi in l r; do 
		  MapIcosahedron -overwrite -spec ${subjid}_${hemi}h.spec -ld $ld -prefix ico${ld}_ || exit 1
		done
		ls
	;;

	fixhead)
	# removes allineation information in AFNI, as this messes up viewing the surfaces aligned
		cd $aligndir || exit 1
		
		anatin="${subjid}_SurfVol_al+orig"
		anatout="${subjid}_SurfVol_al_nob2s+orig"
		
		rm ${anatout}.????
		3dcopy $anatin $anatout
				
		id="1 0 0 0 0 1 0 0 0 0 1 0"
		
		headers="ALLINEATE_MATVEC_B2S_000000 ALLINEATE_MATVEC_S2B_000000 WARPDRIVE_MATVEC_FOR_000000 WARPDRIVE_MATVEC_INV_000000"
		
		for header in $headers; do
			3drefit -atrfloat $header "${id}" $anatout || exit 1
		done
		
		3drefit -oblique_origin ${anatout}
		
	;;

	apply)
	# applies the allineation information to the surfaces,
		cd $aligndir || exit 1
		matrixfile=${aligndir}/${subjid}_SurfVol_al_mat_inv.aff12.1D
		if [ ! -e $matrixfile ]; then
			echo "Cannot find $matrixfile"
			exit 1
		fi
		
		swap="MATRIX(-1,0,0,0,0,-1,0,0,0,0,1,0)" # RAI <-> LPI
		cat_matvec $swap $matrixfile  $swap > ${matrixfile}.3x4.1D
		
		cd ${sumadir} || exit 1
		for k in `ls ico${ld}_?h.*.asc`; do # apply to all surfaces
			echo $k
			m=`echo $k | sed s/\\.asc/_al\\.asc/g`
			echo "Aligning $k -> $m"
			if [ -e "$m" ]; then 
				rm $m
			fi
			ConvertSurface -i_fs $k -o_fs $m -ixmat_1D ${matrixfile}.3x4.1D
			if [ -e "${fsdir}/${m}" ]; then
				rm "${fsdir}/${m}"
			fi
			mv ${m} ${aligndir}
		done
	;;
	
	tlrc)
		# finds the talairach transformation for the anatomical
		cd $aligndir
		anat="${subjid}_SurfVol_al_nob2s"

		rm ${anat}_tlrc+tlrc.????
		@auto_tlrc -base TT_N27+tlrc. -input ${anat}+orig -suffix _tlrc -no_ss || exit 1
	;;
	
	applytlrcvol)
		# applies talairach transformation for anatomical 
		anatin="${subjid}_SurfVol_al_nob2s_tlrc+tlrc"
		anatout="${subjid}_SurfVol_al_nob2s_tlrc_nomat+tlrc"

		surfinfixes="smoothwm pial inflated"
		cd $aligndir
		
		rm ${anatout}.????
		3dcopy $anatin $anatout 

		swap="MATRIX(-1,0,0,0,0,-1,0,0,0,0,1,0)"
		cat_matvec ${swap} ${anatin}::WARPDRIVE_MATVEC_FOR_000000 ${swap} >al2tlrc.aff12.1D
			
		# remove MATVEC vectors 
		id="1 0 0 0 0 1 0 0 0 0 1 0"
		3drefit -atrfloat WARPDRIVE_MATVEC_FOR_000000 "${id}" -atrfloat WARPDRIVE_MATVEC_INV_000000 "${id}" ${anatout} || exit 1
	;;
	
	applytlrcsurf)		
		cd $aligndir	
		surfinfixes="smoothwm pial inflated"

		for surfinfix in $surfinfixes; do # for each surface name
			cd $aligndir
			
			pat="ico${ld}*${surfinfix}_al*asc" 
			echo "$surfinfix: $pat"
			ls *.asc
		 	surffiles=`ls -1 ico${ld}*${surfinfix}_al.asc` #for each surface file (left and right hemisphere) 
		 	
			for surffile in ${surffiles}; do
				surffile_subj="${surffile}"
			
				surffile_subj_tlrc=`echo ${surffile_subj} | sed s/\\.asc/_tlrc\\.asc/g`
				echo $surffile_subj_tlrc

				rm ${surffile_subj_tlrc}
				ConvertSurface -i $surffile_subj -o $surffile_subj_tlrc -ixmat_1D al2tlrc.aff12.1D || exit 1
			done
		done
	;;
		
		
	avgsurf)
		cd $aligndir
		pwd=`pwd`
		
		infixes="smoothwm_pial_avg smoothwm pial"

		for hemi in l r; do
			for suffix in _al _al_tlrc; do 
				args="freesurfer_asc_averagesurfaces("
				for j in 1 2 3; do
					infix=`getfield $j $infixes`
					file="${pwd}/ico${ld}_${hemi}h.${infix}${suffix}.asc"
					if [ $j -gt 1 ]; then
						args="${args},"
					fi
					args="${args}'${file}'"
				done
				args="${args})"
				
				cd $matlabdir
				runmatlab $args 
			done
		done
	;;		
	
	makespec)
		cd $aligndir || exit 1
		surfinfixes="smoothwm smoothwm_pial_avg pial inflated" #critical that smoothm is first!
		groupprefix=all
		viewlist="orig tlrc" # use "orig tlrc" or "tlrc" to get orig and/or tlrc views
		
		for view in $viewlist; do
			if [ $view = tlrc ]; then
				postfix="_tlrc"
				anatpost="_tlrc_nomat"
			else
				postfix=""
				anatpost=""
			fi
			
			for hemi in l r; do
				surfaces=
				specout="spec_${hemi}h_ico${ld}_${view}.spec"
				
				echo > $specout
				echo "# Created using script $0 on "`date` >> $specout
				echo -e "\nGroup = all\n" >> $specout
				
				for surfinfix in $surfinfixes; do # for each surface name
					surffile="ico${ld}_${hemi}h.${surfinfix}_al${postfix}.asc"
					
					statename="${hemi}h_${surfinfix}"
					skiphemistate="X"
					
					if [ $surfinfix = "smoothwm" ]; then # was smoothwm
						ldp="SAME"
						ldpname=${surffile}
					else
						ldp=$ldpname
					fi
									
					surfaces="${surfaces}\nNewSurface\nSurfaceFormat = ASCII\nSurfaceType = FreeSurfer\n"
					surfaces="${surfaces}FreeSurferSurface = ${surffile}\nLocalDomainParent = ${ldp}\n"
					surfaces="${surfaces}SurfaceState = ${statename}\nEmbedDimension = 3\n\n"
					
					if [ ${hemi} = $skiphemistate ]; then
						continue
					fi
					
					echo "StateDef = $statename" >> $specout
				done
				echo -e $surfaces >> $specout
				echo "Spec for ${hemi}-hemi written to "`pwd`"/$specout"
				
				seesumaout="${hemi}h_seesuma_ico${ld}_${view}.sh"
				
				echo -e "killall afni\nafni -niml &\nsuma -spec ${specout} -sv ${subjid}_SurfVol_al_nob2s${anatpost}+${view}" > ${seesumaout}
				chmod u+x ${seesumaout}
				
				echo "To view ${hemi}-hemi surface, run ./${seesumaout}"
			done
		done
	
	;;

	linkup) #links files from aligndir to fsdir - only tlrc ones
		surfinfixes="smoothwm smoothwm_pial_avg pial inflated"
		
		cd ${fsdir}
		
		linkheadbrik ${aligndir}/${subjid}_SurfVol_al_nob2s_tlrc_nomat+tlrc .
		
		cd ${aligndir}

		for j in `ls spec*spec ?h_seesuma_ico${ld}_tlrc.sh`; do
			echo "Linking $j"
			t=${fsdir}/${j}
			rm ${t}
			ln ${j} ${t}
		done
				
		for surfinfix in $surfinfixes; do 
			for j in `ls ico*_?h.${surfinfix}_al_tlrc.asc`; do
				t=${fsdir}/${j}
				rm ${t}
				ln ${j} ${t}
			done
		done
					
	;;
	
	groupavgsurf)
		if [ ! -e ${groupfsdir} ]; then
			mkdir $groupfsdir
		fi
	
		surfinfixes="pial smoothwm smoothwm_pial_avg inflated"
		
		subjone=`echo $subjids | cut -f1 -d' '`
		for surfinfix in $surfinfixes; do # for each surface name
			for hemi in l r; do 
				cd ${groupfsdir}
				pwd=`pwd`"/"
				fn="ico${ld}_${hemi}h.${surfinfix}_al_tlrc.asc" 
				echo $fn
				
				args="freesurfer_asc_averagesurfaces('${pwd}${fn}'"
				for subjnum in $subjnums; do
					subjid="${subjprefix}${subjnum}"
					subjfn="${groupfsdir}../${subjid}/align/${fn}"
					
					if [ ! -e ${subjfn} ]; then
						echo "File $subjfn not found!"
						exit 1
					fi
					
					args="${args},'${subjfn}'"
				done
				args="${args})"
				echo $args
				cd $matlabdir

				runmatlab $args
			 done
		done
	;;
	
	groupmakespec)
		cd $groupfsdir
		
		anatpostfix="_SurfVol_al_nob2s_tlrc_nomat"
		anatfile="all_avg${anatpostfix}"
		anatfns=
		if [ ! -e ${anatfile}+tlrc.BRIK ]; then
			for subjnum in $subjnums; do
				anatfn="${subjprefix}${subjnum}${anatpostfix}+tlrc"
				anatfnpath="${fsrootdir}/${subjprefix}${subjnum}/align/${anatfn}"
				linkheadbrik ${anatfnpath} . 
				anatfns="${anatfns} $anatfn"
			done
			
			3dMean -prefix ${anatfile} ${anatfns} || exit 1
		fi
		
		surfinfixes="smoothwm smoothwm_pial_avg pial inflated" #critical that smoothm is first!
		groupprefix=all
		for hemi in l r; do
			surfaces=
			specout="${groupprefix}_${hemi}h_ico${ld}.spec"
			
			echo > $specout
			echo "# Created using script $0 on "`date` >> $specout
			echo -e "\nGroup = all\n" >> $specout
			
			for surfinfix in $surfinfixes; do # for each surface name
				surffile="ico${ld}_${hemi}h.${surfinfix}_al_tlrc.asc"
				
				statename="${hemi}h_${surfinfix}"
				skiphemistate="X"
				
				if [ $surfinfix = "smoothwm" ]; then # was smoothwm
					ldp="SAME"
					ldpname=${surffile}
				else
					ldp=$ldpname
				fi
								
				surfaces="${surfaces}\nNewSurface\nSurfaceFormat = ASCII\nSurfaceType = FreeSurfer\n"
				surfaces="${surfaces}FreeSurferSurface = ${surffile}\nLocalDomainParent = ${ldp}\n"
				surfaces="${surfaces}SurfaceState = ${statename}\nEmbedDimension = 3\n\n"
				
				if [ ${hemi} = $skiphemistate ]; then
					continue
				fi
				
				echo "StateDef = $statename" >> $specout
			done
			echo -e $surfaces >> $specout
			echo "Spec for ${hemi}-hemi written to "`pwd`"/$specout"
			
			seesumaout="${hemi}h_seesuma_ico${ld}.sh"
			echo -e "killall afni\nafni -niml &\nsuma -spec ${specout} -sv ${anatfile}+tlrc" > ${seesumaout}
			chmod u+x ${seesumaout}
			
			echo "To view ${hemi}-hemi surface, run ./${seesumaout}"


		done
	;;
		
	all)
		shift 
		shift
		for i in tosuma align mapico fixhead apply tlrc applytlrcvol applytlrcsurf avgsurf makespec linkup; do
#		for i in avgsurf makespec linkup; do

			cd $scriptdir
			echo -e "\n\n*** ${0}: Running $i ***\n\n"
			$0 $subjnum $i $@ || exit 1
		done
		exit 0
	;;
		
	ico)
		shift 
		shift
		for i in mapico apply applytlrcsurf avgsurf makespec linkup; do
			cd $scriptdir
			echo "Running $0 $subjnum $i $@ for $subjnum"
			$0 $subjnum $i $@ || exit 1
		done
		exit 0
	;;
	
	groupall)
		for i in groupavgsurf groupmakespec; do
			cd $scriptdir
			$0 $1 $i || exit 1
		done
		exit 0
	;;
	
	trytlrc)
		tstdir=${fsrootdir}/tst/
		mkdir $tstdir
		cd $tstdir
		
		matfn="mat.3x4.1D"
		head -11 ${fsrootdir}/${subjid}/${subjid}/mri/transforms/talairach.lta | tail -3 > $matfn
		cat $matfn
		pwd
		
	;;
		
	

	*)
		echo "Unknown action $action"
		exit 1
	;;
esac
	
	
