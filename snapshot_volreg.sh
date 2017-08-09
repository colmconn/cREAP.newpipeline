#!/bin/bash

## set -x

programName=`basename $0`

DISPLAY_FILE=$( mktemp )

# control variables

AFNI_DONT_LOGFILE=YES
# unset noclobber

# check if all the needed auxiliary programs exist

(( nerr=0 )) 
errm="** ERROR:"

plist="Xvfb djpeg cjpeg pnmcat pbmtext pnmscale pbmtopgm"
for pppp in $plist ; do
  wwww=`which $pppp`
  if [[ $? != 0 ]] ; then
    (( nerr=nerr+1 ))
    errm="$errm $pppp"
  fi
done

# we can use either pamcomp or pnmcomp, so only need to
# find one of the twain

pcprog=pamcomp
wwww=`which $pcprog`
if [[ $? != 0 ]] ; then
    pcprog=pnmcomp
  wwww=`which $pcprog`
  if [[ $? != 0 ]] ; then
    (( nerr=nerr+1 ))
    errm="$errm (pnmcomp OR pamcomp)"
  fi
fi

if [[ $nerr -gt 0 ]] ; then
  echo "$errm -- not found in path -- @snapshot_volreg fails"
  echo "** WARNING: this script cannot run without installing package netpbm11"
  exit 1
fi

# set the prefix for the anat and epi datasets

adset=$1
abase=${adset##*/}
anat=`basename $abase .gz`
anat=`basename $anat  .nii`
anat=`basename $anat  .HEAD`
anat=`basename $anat  .BRIK`
anat=`basename $anat  +tlrc`
anat=`basename $anat  +orig`
anat=`basename $anat  +acpc`
anat=`basename $anat  +tlrc.`
anat=`basename $anat  +orig.`
anat=`basename $anat  +acpc.`
if [[ $abase == $anat.nii.gz ]] ; then
  asuff=".nii.gz"
elif [[ $abase == $anat.nii ]] ; then
  asuff=.nii
else
  asuff=""
fi

edset=$2
epi=`basename $edset .gz`
epi=`basename $epi   .nii`
epi=`basename $epi   .HEAD`
epi=`basename $epi   .BRIK`
epi=`basename $epi   +orig`
epi=`basename $epi   +acpc`
epi=`basename $epi   +tlrc`
epi=`basename $epi   +tlrc.`
anat=`basename $anat  +tlrc.`
anat=`basename $anat  +tlrc`
anat=`basename $anat  +orig.`
anat=`basename $anat  +acpc.`

# set output image prefix

if [[ $# -gt 2  ]] ; then
    jnam=$3
    jnam=`basename "$jnam" .jpg`
    jnam=`basename "$jnam" .JPG`
    jnam=`basename "$jnam" .jpeg`
    jnam=`basename "$jnam" .JPEG`
else
  jnam=$anat.$epi
fi

# Are we re-using Xvfb?

if [[ $# -gt 3 ]] ; then
    xdisplay=$4
    if xdpyinfo -display $xdisplay > /dev/null 2>&1 ; then
	echo "Found X display number $xdisplay"
    else 
	echo "** ERROR: it doesn't look like Xvfb is running with display $xdisplay"
	exit 1
    fi
fi

exad=`3dinfo -exists $adset`
exep=`3dinfo -exists $edset`

if [[ $exad == 0 || $exep == 0 ]] ; then
    if [[ $exad == 0 ]] ; then
	echo "** ERROR: @snapshot_volreg can't find $adset"
    fi
    if [[ $exep == 0 ]] ;  then
	echo "** ERROR: @snapshot_volreg can't find $edset"
    fi
    exit 1
fi

# set some AFNI GUI things

export AFNI_NOSPLASH=YES
export AFNI_SPLASH_MELT=NO
export AFNI_LEFT_IS_LEFT=YES
export AFNI_IMAGE_LABEL_MODE=5
export AFNI_IMAGE_LABEL_SIZE=2
export AFNI_ENVIRON_WARNINGS=NO
export AFNI_COMPRESSOR=NONE
export OMP_NUM_THREADS=1

# start the X virtual frame buffer on display #xdisplay


if [[ "x$xdisplay" == "x" ]] ; then
    killX=1
    (( ntry=0 ))
    Xnotfound=1
    while  [[ $Xnotfound==1 ]] ; do
	# if( -e /tmp/.X${xdisplay}-lock ) continue
	# # echo " -- trying to start Xvfb :${xdisplay}"
	echo " -- trying to start Xvfb"
	exec 6>$DISPLAY_FILE
	Xvfb -displayfd 6 -screen 0 1024x768x24 >& /dev/null &
	Xpid=$!
	sleep 1
	xdisplay=$( cat $DISPLAY_FILE )
	if xdpyinfo -display :$xdisplay > /dev/null 2>&1 ; then
	    echo "Found X display number $xdisplay"
	    Xnotfound=0
	else 
	    echo "** ERROR: it doesn't look like Xvfb is running with display $xdisplay"
	fi
	if [[ $Xnotfound == 0 ]] ; then break ; fi
	(( ntry=ntry+1 ))
	if [[ $ntry -ge 5 ]] ; then
	    echo "** ERROR: can't start Xvfb -- exiting"
	    exit 1
	fi
    done
fi

## exit
export DISPLAY=:${xdisplay}

ranval=$( echo $RANDOM | shasum - | awk '{print $1}' )

# quasi-random temp filename prefix

zpref=zzerm.X${xdisplay}-${ranval}

# crop the input anat to a standard size

cdset=${zpref}.acrop.nii

3dAutobox -npad 17 -prefix $cdset $adset

# resample the EPI to the anat grid

3dAllineate -input  ${edset}'[0]'     \
            -master ${cdset}          \
            -prefix ${zpref}.epiR.nii \
            -1Dparam_apply IDENTITY   \
            -final cubic
# get the lower and upper thresholds for the edgized EPI

epp=( $( 3dBrickStat -non-zero -percentile 20 60 80 ${zpref}.epiR.nii ) )
eth=${epp[1]}
emx=${epp[3]}

# run AFNI to make the 3 overlay images

anatCM=( `3dCM $cdset` )

anatNN=`3dinfo -nijk $cdset`

astep=`ccalc -int "cbrt($anatNN)/6.111"`

if [[ $astep < 2 ]] ; then astep = 2; fi

#     -com "SET_FUNC_RANGE $emx"                                  \
#     -com "SET_THRESHNEW $eth *"                                 \


afni -noplugins -no_detach                                       \
     -com "SWITCH_UNDERLAY ${cdset}"                             \
     -com "SWITCH_OVERLAY ${zpref}.epiR.nii"                     \
     -com "SET_DICOM_XYZ ${anatCM[*]}"                           \
     -com "SEE_OVERLAY +"                                        \
     -com "SET_XHAIRS OFF"                                       \
     -com "OPEN_WINDOW sagittalimage opacity=7 mont=3x1:$astep"  \
     -com "OPEN_WINDOW axialimage opacity=7 mont=3x1:$astep"     \
     -com "OPEN_WINDOW coronalimage opacity=7 mont=3x1:$astep"   \
     -com "SAVE_JPEG sagittalimage ${zpref}.sag.jpg blowup=2"    \
     -com "SAVE_JPEG coronalimage  ${zpref}.cor.jpg blowup=2"    \
     -com "SAVE_JPEG axialimage    ${zpref}.axi.jpg blowup=2"    \
     -com "QUITT"                                                \
     ${cdset} ${zpref}.epiR.nii

# convert the output JPEGs to PNMs for manipulation

djpeg ${zpref}.sag.jpg > ${zpref}.sag.pnm
djpeg ${zpref}.cor.jpg > ${zpref}.cor.pnm
djpeg ${zpref}.axi.jpg > ${zpref}.axi.pnm

# cat them together, make a text label, overlay it, make output JPEG

pnmcat -tb -jcenter -black ${zpref}.sag.pnm ${zpref}.axi.pnm ${zpref}.cor.pnm           > ${zpref}.pnm
pbmtext -builtin fixed "$jnam"   | pnmcrop | pbmtopgm 1 1 | pnmscale 2                  > ${zpref}.t1.pgm
pbmtext -builtin fixed "<- Left"           | pbmtopgm 1 1                               > ${zpref}.t2.pgm
$pcprog -align=right -valign=bottom ${zpref}.t1.pgm ${zpref}.pnm                        > ${zpref}.t3.pnm
$pcprog -align=right -valign=top    ${zpref}.t2.pgm ${zpref}.t3.pnm | cjpeg -quality 95 > "$jnam.jpg"

# delete the trash data

\rm -f ${zpref}.*

echo "$programName output image = $jnam.jpg"

# stop Xvfb if we started it ourselves

if [[ $killX==1 ]] ; then
    kill -9 $Xpid
    rm -f $DISPLAY_FILE
fi

exit 0
