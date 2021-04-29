# echo "setup LambdacAlg LambdacAlg-00-00-00 in /publicfs/ucas/user/yuanchy8/workbefs/7.0.6/LambdacAlg-etaprime"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtLambdacAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtLambdacAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=LambdacAlg -version=LambdacAlg-00-00-00 -path=/publicfs/ucas/user/yuanchy8/workbefs/7.0.6/LambdacAlg-etaprime  -no_cleanup $* >${cmtLambdacAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=LambdacAlg -version=LambdacAlg-00-00-00 -path=/publicfs/ucas/user/yuanchy8/workbefs/7.0.6/LambdacAlg-etaprime  -no_cleanup $* >${cmtLambdacAlgtempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtLambdacAlgtempfile}
  unset cmtLambdacAlgtempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtLambdacAlgtempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtLambdacAlgtempfile}
unset cmtLambdacAlgtempfile
exit $cmtsetupstatus

