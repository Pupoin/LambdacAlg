# echo "setup SimplePIDSvc SimplePIDSvc-00-00-16 in /workfs/bes/yuanchy8/7.0.5"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtSimplePIDSvctempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtSimplePIDSvctempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=SimplePIDSvc -version=SimplePIDSvc-00-00-16 -path=/workfs/bes/yuanchy8/7.0.5  -no_cleanup $* >${cmtSimplePIDSvctempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=SimplePIDSvc -version=SimplePIDSvc-00-00-16 -path=/workfs/bes/yuanchy8/7.0.5  -no_cleanup $* >${cmtSimplePIDSvctempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtSimplePIDSvctempfile}
  unset cmtSimplePIDSvctempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtSimplePIDSvctempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtSimplePIDSvctempfile}
unset cmtSimplePIDSvctempfile
exit $cmtsetupstatus

