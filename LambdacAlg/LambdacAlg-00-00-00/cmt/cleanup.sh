# echo "cleanup LambdacAlg LambdacAlg-00-00-00 in /publicfs/ucas/user/yuanchy8/workbefs/7.0.6/LambdacAlg"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtLambdacAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtLambdacAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=LambdacAlg -version=LambdacAlg-00-00-00 -path=/publicfs/ucas/user/yuanchy8/workbefs/7.0.6/LambdacAlg  $* >${cmtLambdacAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=LambdacAlg -version=LambdacAlg-00-00-00 -path=/publicfs/ucas/user/yuanchy8/workbefs/7.0.6/LambdacAlg  $* >${cmtLambdacAlgtempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtLambdacAlgtempfile}
  unset cmtLambdacAlgtempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtLambdacAlgtempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtLambdacAlgtempfile}
unset cmtLambdacAlgtempfile
return $cmtcleanupstatus

