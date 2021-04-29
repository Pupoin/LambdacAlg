# echo "setup LambdacAlg LambdacAlg-00-00-00 in /publicfs/ucas/user/yuanchy8/workbefs/7.0.6/LambdacAlgtmp"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtLambdacAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtLambdacAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=LambdacAlg -version=LambdacAlg-00-00-00 -path=/publicfs/ucas/user/yuanchy8/workbefs/7.0.6/LambdacAlgtmp  -no_cleanup $* >${cmtLambdacAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=LambdacAlg -version=LambdacAlg-00-00-00 -path=/publicfs/ucas/user/yuanchy8/workbefs/7.0.6/LambdacAlgtmp  -no_cleanup $* >${cmtLambdacAlgtempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtLambdacAlgtempfile}
  unset cmtLambdacAlgtempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtLambdacAlgtempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtLambdacAlgtempfile}
unset cmtLambdacAlgtempfile
return $cmtsetupstatus

