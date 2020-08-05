if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/ihep.ac.cn/bes3/offline/ExternalLib/contrib/CMT/v1r20p20081118
endif
source ${CMTROOT}/mgr/setup.csh
set tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set tempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=LambdacSigmapEta -version=LambdacSigmapEta-00-00-01 -path=/workfs/bes/yuanchy8/6.6.4.p01 $* >${tempfile}; source ${tempfile}
/bin/rm -f ${tempfile}

