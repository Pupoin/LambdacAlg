# echo "setup TestRelease TestRelease-00-00-92 in /workfs/bes/yuanchy8/7.0.5"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtTestReleasetempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtTestReleasetempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=TestRelease -version=TestRelease-00-00-92 -path=/workfs/bes/yuanchy8/7.0.5  -no_cleanup $* >${cmtTestReleasetempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=TestRelease -version=TestRelease-00-00-92 -path=/workfs/bes/yuanchy8/7.0.5  -no_cleanup $* >${cmtTestReleasetempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtTestReleasetempfile}
  unset cmtTestReleasetempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtTestReleasetempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtTestReleasetempfile}
unset cmtTestReleasetempfile
return $cmtsetupstatus

