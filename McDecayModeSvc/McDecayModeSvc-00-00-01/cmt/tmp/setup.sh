# echo "setup McDecayModeSvc McDecayModeSvc-00-00-01 in /workfs/bes/yuanchy8/7.0.5"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtMcDecayModeSvctempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtMcDecayModeSvctempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=McDecayModeSvc -version=McDecayModeSvc-00-00-01 -path=/workfs/bes/yuanchy8/7.0.5  -no_cleanup $* >${cmtMcDecayModeSvctempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=McDecayModeSvc -version=McDecayModeSvc-00-00-01 -path=/workfs/bes/yuanchy8/7.0.5  -no_cleanup $* >${cmtMcDecayModeSvctempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtMcDecayModeSvctempfile}
  unset cmtMcDecayModeSvctempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtMcDecayModeSvctempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtMcDecayModeSvctempfile}
unset cmtMcDecayModeSvctempfile
return $cmtsetupstatus

