# echo "cleanup LambdacAlg LambdacAlg-00-00-00 in /workfs/bes/yuanchy8/7.0.5"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtLambdacAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtLambdacAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=LambdacAlg -version=LambdacAlg-00-00-00 -path=/workfs/bes/yuanchy8/7.0.5  $* >${cmtLambdacAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=LambdacAlg -version=LambdacAlg-00-00-00 -path=/workfs/bes/yuanchy8/7.0.5  $* >${cmtLambdacAlgtempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtLambdacAlgtempfile}
  unset cmtLambdacAlgtempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtLambdacAlgtempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtLambdacAlgtempfile}
unset cmtLambdacAlgtempfile
exit $cmtcleanupstatus

