# echo "cleanup SimplePIDSvc SimplePIDSvc-00-00-16 in /workfs/bes/yuanchy8/7.0.5"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtSimplePIDSvctempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtSimplePIDSvctempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=SimplePIDSvc -version=SimplePIDSvc-00-00-16 -path=/workfs/bes/yuanchy8/7.0.5  $* >${cmtSimplePIDSvctempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=SimplePIDSvc -version=SimplePIDSvc-00-00-16 -path=/workfs/bes/yuanchy8/7.0.5  $* >${cmtSimplePIDSvctempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtSimplePIDSvctempfile}
  unset cmtSimplePIDSvctempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtSimplePIDSvctempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtSimplePIDSvctempfile}
unset cmtSimplePIDSvctempfile
return $cmtcleanupstatus

