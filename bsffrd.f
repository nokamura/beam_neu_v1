      subroutine bsffrd
      implicit none
      include 'inc/hmffdt.inc'
      include 'inc/bsffcm.inc'
      include 'inc/bsfile.inc'
      include 'inc/crndm.inc'
      include 'inc/hmunit.inc'
      include 'inc/hmflag.inc'
      include 'inc/params.inc'
      integer ikey
      integer i,lnblnk

      do ikey=1,nkeys
         if(ckey(ikey).eq.'basic_binsize') then
            read(cval(ikey),*) basic_binsize
         else if(ckey(ikey).eq.'binsize_factor') then
            read(cval(ikey),*) binsize_factor
         else if(ckey(ikey).eq.'l') then
            read(cval(ikey),*) L
         else if(ckey(ikey).eq.'ndiv') then
            read(cval(ikey),*) ndiv
         else if(ckey(ikey).eq.'beam') then
            read(cval(ikey),*) beam
         else if(ckey(ikey).eq.'oab') then
            read(cval(ikey),*) OAB
         else if(ckey(ikey).eq.'nu_mode') then
            read(cval(ikey),*) nu_mode
         else if(ckey(ikey).eq.'detect') then
            read(cval(ikey),*) detect
         else if(ckey(ikey).eq.'v') then
            read(cval(ikey),*) V
         else if(ckey(ikey).eq.'y') then
            read(cval(ikey),*) Y
         else if(ckey(ikey).eq.'r_nu') then
            read(cval(ikey),*) r_nu
         else if(ckey(ikey).eq.'r_anu') then
            read(cval(ikey),*) r_anu
         else if(ckey(ikey).eq.'mode') then
            read(cval(ikey),*) mode
         else if(ckey(ikey).eq.'ifixl') then
            read(cval(ikey),*) ifixL
         else if(ckey(ikey).eq.'ifluc') then
            read(cval(ikey),*) ifluc
         else if(ckey(ikey).eq.'ihisto') then
            read(cval(ikey),*) ihisto
         else if(ckey(ikey).eq.'ismear') then
            read(cval(ikey),*) ismear
         else if(ckey(ikey).eq.'iccqe') then
            read(cval(ikey),*) iCCQE
         else if(ckey(ikey).eq.'iccres') then
            read(cval(ikey),*) iCCRes
         else if(ckey(ikey).eq.'inc') then
            read(cval(ikey),*) iNC
         else if(ckey(ikey).eq.'ipi0dist') then
            read(cval(ikey),*) ipi0dist
         else if(ckey(ikey).eq.'ipi0xsec') then
            read(cval(ikey),*) ipi0xsec
         else if(ckey(ikey).eq.'iminuit') then
            read(cval(ikey),*) iminuit
         else if(ckey(ikey).eq.'ipi0unc') then
            read(cval(ikey),*) ipi0unc
         else if(ckey(ikey).eq.'ihypo') then
            read(cval(ikey),*) ihypo
         else if(ckey(ikey).eq.'iminevents') then
            read(cval(ikey),*) iminevents
         else if(ckey(ikey).eq.'innm') then
            read(cval(ikey),*) iNnm
         else if(ckey(ikey).eq.'inne') then
            read(cval(ikey),*) iNne
         else if(ckey(ikey).eq.'inam') then
            read(cval(ikey),*) iNam
         else if(ckey(ikey).eq.'inae') then
            read(cval(ikey),*) iNae
         else if(ckey(ikey).eq.'isk') then
            read(cval(ikey),*) iSK
         else if(ckey(ikey).eq.'srho') then
            read(cval(ikey),*) Srho
         else if(ckey(ikey).eq.'sv') then
            read(cval(ikey),*) SV
         else if(ckey(ikey).eq.'sl') then
            read(cval(ikey),*) SL
         else if(ckey(ikey).eq.'soab') then
            read(cval(ikey),*) SOAB
         else if(ckey(ikey).eq.'semin') then
            read(cval(ikey),*) Semin
         else if(ckey(ikey).eq.'semax') then
            read(cval(ikey),*) Semax
         else if(ckey(ikey).eq.'smmin') then
            read(cval(ikey),*) Smmin
         else if(ckey(ikey).eq.'smmax') then
            read(cval(ikey),*) Smmax
         else if(ckey(ikey).eq.'sbfactor') then
            read(cval(ikey),*) Sbfactor
         else if(ckey(ikey).eq.'ioki') then
            read(cval(ikey),*) iOki
         else if(ckey(ikey).eq.'orho') then
            read(cval(ikey),*) Orho
         else if(ckey(ikey).eq.'ov') then
            read(cval(ikey),*) OV
         else if(ckey(ikey).eq.'ol') then
            read(cval(ikey),*) OL
         else if(ckey(ikey).eq.'ooab') then
            read(cval(ikey),*) OOAB
         else if(ckey(ikey).eq.'oemin') then
            read(cval(ikey),*) Oemin
         else if(ckey(ikey).eq.'oemax') then
            read(cval(ikey),*) Oemax
         else if(ckey(ikey).eq.'ommin') then
            read(cval(ikey),*) Ommin
         else if(ckey(ikey).eq.'ommax') then
            read(cval(ikey),*) Ommax
         else if(ckey(ikey).eq.'obfactor') then
            read(cval(ikey),*) Obfactor
         else if(ckey(ikey).eq.'ikr') then
            read(cval(ikey),*) iKr
         else if(ckey(ikey).eq.'krho') then
            read(cval(ikey),*) Krho
         else if(ckey(ikey).eq.'kv') then
            read(cval(ikey),*) KV
         else if(ckey(ikey).eq.'kl') then
            read(cval(ikey),*) KL
         else if(ckey(ikey).eq.'koab') then
            read(cval(ikey),*) KOAB
         else if(ckey(ikey).eq.'kemin') then
            read(cval(ikey),*) Kemin
         else if(ckey(ikey).eq.'kemax') then
            read(cval(ikey),*) Kemax
         else if(ckey(ikey).eq.'kmmin') then
            read(cval(ikey),*) Kmmin
         else if(ckey(ikey).eq.'kmmax') then
            read(cval(ikey),*) Kmmax
         else if(ckey(ikey).eq.'kbfactor') then
            read(cval(ikey),*) Kbfactor
         else if(ckey(ikey).eq.'serror') then
            read(cval(ikey),*) serror
         else if(ckey(ikey).eq.'snmax') then
            read(cval(ikey),*) snmax
         else if(ckey(ikey).eq.'ihfunc') then
            read(cval(ikey),*) ihfunc
         else if(ckey(ikey).eq.'ieff') then
            read(cval(ikey),*) ieff
         else if(ckey(ikey).eq.'imissid') then
            read(cval(ikey),*) imissid
         else if(ckey(ikey).eq.'iinput') then
            read(cval(ikey),*) iinput
         else if(ckey(ikey).eq.'ithatm') then
            read(cval(ikey),*) ithatm
         else if(ckey(ikey).eq.'ioct') then
            read(cval(ikey),*) ioct
         else if(ckey(ikey).eq.'idmatm2') then
            read(cval(ikey),*) idmatm2
         else if(ckey(ikey).eq.'emin') then
            read(cval(ikey),*) Emin
         else if(ckey(ikey).eq.'emax') then
            read(cval(ikey),*) Emax
         else if(ckey(ikey).eq.'s2sol_2') then
            read(cval(ikey),*) s2sol_2
         else if(ckey(ikey).eq.'thatm') then
            read(cval(ikey),*) thatm
         else if(ckey(ikey).eq.'s2rct_2') then
            read(cval(ikey),*) s2rct_2
         else if(ckey(ikey).eq.'dm21_2') then
            read(cval(ikey),*) dm21_2
         else if(ckey(ikey).eq.'dmatm_2') then
            read(cval(ikey),*) dmatm_2
         else if(ckey(ikey).eq.'dcp') then
            read(cval(ikey),*) dCP
         else if(ckey(ikey).eq.'fs2sol_2') then
            read(cval(ikey),*) fs2sol_2(1)
         else if(ckey(ikey).eq.'err_s2sol_2') then
            read(cval(ikey),*) fs2sol_2(2)
         else if(ckey(ikey).eq.'ifit_s2sol_2') then
            read(cval(ikey),*) fs2sol_2(3)
         else if(ckey(ikey).eq.'ichi2_s2sol_2') then
            read(cval(ikey),*) fs2sol_2(4)
         else if(ckey(ikey).eq.'fthatm') then
            read(cval(ikey),*) fthatm(1)
         else if(ckey(ikey).eq.'err_thatm') then
            read(cval(ikey),*) fthatm(2)
         else if(ckey(ikey).eq.'ifit_thatm') then
            read(cval(ikey),*) fthatm(3)
         else if(ckey(ikey).eq.'ichi2_thatm') then
            read(cval(ikey),*) fthatm(4)
         else if(ckey(ikey).eq.'fs2rct_2') then
            read(cval(ikey),*) fs2rct_2(1)
         else if(ckey(ikey).eq.'err_s2rct_2') then
            read(cval(ikey),*) fs2rct_2(2)
         else if(ckey(ikey).eq.'ifit_s2rct_2') then
            read(cval(ikey),*) fs2rct_2(3)
         else if(ckey(ikey).eq.'ichi2_s2rct_2') then
            read(cval(ikey),*) fs2rct_2(4)
         else if(ckey(ikey).eq.'fdm21_2') then
            read(cval(ikey),*) fdm21_2(1)
         else if(ckey(ikey).eq.'err_dm21_2') then
            read(cval(ikey),*) fdm21_2(2)
         else if(ckey(ikey).eq.'ifit_dm21_2') then
            read(cval(ikey),*) fdm21_2(3)
         else if(ckey(ikey).eq.'ichi2_dm21_2') then
            read(cval(ikey),*) fdm21_2(4)
         else if(ckey(ikey).eq.'fdmatm_2') then
            read(cval(ikey),*) fdmatm_2(1)
         else if(ckey(ikey).eq.'err_dmatm_2') then
            read(cval(ikey),*) fdmatm_2(2)
         else if(ckey(ikey).eq.'ifit_dmatm_2') then
            read(cval(ikey),*) fdmatm_2(3)
         else if(ckey(ikey).eq.'ichi2_dmatm_2') then
            read(cval(ikey),*) fdmatm_2(4)
         else if(ckey(ikey).eq.'fdcp') then
            read(cval(ikey),*) fdCP(1)
         else if(ckey(ikey).eq.'err_dcp') then
            read(cval(ikey),*) fdCP(2)
         else if(ckey(ikey).eq.'ifit_dcp') then
            read(cval(ikey),*) fdCP(3)
         else if(ckey(ikey).eq.'ichi2_dcp') then
            read(cval(ikey),*) fdCP(4)
         else if(ckey(ikey).eq.'mh') then
            read(cval(ikey),*) MH
         else if(ckey(ikey).eq.'rho') then
            read(cval(ikey),*) rho
         else if(ckey(ikey).eq.'fccqen') then
            read(cval(ikey),*) fCCQEn(1)
         else if(ckey(ikey).eq.'err_fccqen') then
            read(cval(ikey),*) fCCQEn(2)
         else if(ckey(ikey).eq.'fccqea') then
            read(cval(ikey),*) fCCQEa(1)
         else if(ckey(ikey).eq.'err_fccqea') then
            read(cval(ikey),*) fCCQEa(2)
         else if(ckey(ikey).eq.'fccresn') then
            read(cval(ikey),*) fCCResn(1)
         else if(ckey(ikey).eq.'err_fccresn') then
            read(cval(ikey),*) fCCResn(2)
         else if(ckey(ikey).eq.'fccresa') then
            read(cval(ikey),*) fCCResa(1)
         else if(ckey(ikey).eq.'err_fccresa') then
            read(cval(ikey),*) fCCResa(2)
         else if(ckey(ikey).eq.'fpi0') then
            read(cval(ikey),*) fpi0(1)
         else if(ckey(ikey).eq.'err_fpi0') then
            read(cval(ikey),*) fpi0(2)
         else if(ckey(ikey).eq.'fpirs') then
            read(cval(ikey),*) fpirs(1)
         else if(ckey(ikey).eq.'err_fpirs') then
            read(cval(ikey),*) fpirs(2)
         else if(ckey(ikey).eq.'fpico') then
            read(cval(ikey),*) fpico(1)
         else if(ckey(ikey).eq.'err_fpico') then
            read(cval(ikey),*) fpico(2)
         else if(ckey(ikey).eq.'frpi0') then
            read(cval(ikey),*) frpi0(1)
         else if(ckey(ikey).eq.'err_frpi0') then
            read(cval(ikey),*) frpi0(2)
         else if(ckey(ikey).eq.'fs') then
            read(cval(ikey),*) fS(1)
         else if(ckey(ikey).eq.'err_fs') then
            read(cval(ikey),*) fS(2)
         else if(ckey(ikey).eq.'fsrho') then
            read(cval(ikey),*) fSrho(1)
         else if(ckey(ikey).eq.'err_fsrho') then
            read(cval(ikey),*) fSrho(2)
         else if(ckey(ikey).eq.'fsnmn') then
            read(cval(ikey),*) fSnmn(1)
         else if(ckey(ikey).eq.'err_fsnmn') then
            read(cval(ikey),*) fSnmn(2)
         else if(ckey(ikey).eq.'fsnen') then
            read(cval(ikey),*) fSnen(1)
         else if(ckey(ikey).eq.'err_fsnen') then
            read(cval(ikey),*) fSnen(2)
         else if(ckey(ikey).eq.'fsnma') then
            read(cval(ikey),*) fSnma(1)
         else if(ckey(ikey).eq.'err_fsnma') then
            read(cval(ikey),*) fSnma(2)
         else if(ckey(ikey).eq.'fsnea') then
            read(cval(ikey),*) fSnea(1)
         else if(ckey(ikey).eq.'err_fsnea') then
            read(cval(ikey),*) fSnea(2)
         else if(ckey(ikey).eq.'fsamn') then
            read(cval(ikey),*) fSamn(1)
         else if(ckey(ikey).eq.'err_fsamn') then
            read(cval(ikey),*) fSamn(2)
         else if(ckey(ikey).eq.'fsaen') then
            read(cval(ikey),*) fSaen(1)
         else if(ckey(ikey).eq.'err_fsaen') then
            read(cval(ikey),*) fSaen(2)
         else if(ckey(ikey).eq.'fsama') then
            read(cval(ikey),*) fSama(1)
         else if(ckey(ikey).eq.'err_fsama') then
            read(cval(ikey),*) fSama(2)
         else if(ckey(ikey).eq.'fsaea') then
            read(cval(ikey),*) fSaea(1)
         else if(ckey(ikey).eq.'err_fsaea') then
            read(cval(ikey),*) fSaea(2)
         else if(ckey(ikey).eq.'seffe') then
            read(cval(ikey),*) Seffe(1)
         else if(ckey(ikey).eq.'err_seffe') then
            read(cval(ikey),*) Seffe(2)
         else if(ckey(ikey).eq.'seffmu') then
            read(cval(ikey),*) Seffmu(1)
         else if(ckey(ikey).eq.'err_seffmu') then
            read(cval(ikey),*) Seffmu(2)
         else if(ckey(ikey).eq.'spe2m') then
            read(cval(ikey),*) SPe2m(1)
         else if(ckey(ikey).eq.'err_spe2m') then
            read(cval(ikey),*) SPe2m(2)
         else if(ckey(ikey).eq.'spm2e') then
            read(cval(ikey),*) SPm2e(1)
         else if(ckey(ikey).eq.'err_spm2e') then
            read(cval(ikey),*) SPm2e(2)
         else if(ckey(ikey).eq.'fo') then
            read(cval(ikey),*) fO(1)
         else if(ckey(ikey).eq.'err_fo') then
            read(cval(ikey),*) fO(2)
         else if(ckey(ikey).eq.'forho') then
            read(cval(ikey),*) fOrho(1)
         else if(ckey(ikey).eq.'err_forho') then
            read(cval(ikey),*) fOrho(2)
         else if(ckey(ikey).eq.'fonmn') then
            read(cval(ikey),*) fOnmn(1)
         else if(ckey(ikey).eq.'err_fonmn') then
            read(cval(ikey),*) fOnmn(2)
         else if(ckey(ikey).eq.'fonen') then
            read(cval(ikey),*) fOnen(1)
         else if(ckey(ikey).eq.'err_fonen') then
            read(cval(ikey),*) fOnen(2)
         else if(ckey(ikey).eq.'fonma') then
            read(cval(ikey),*) fOnma(1)
         else if(ckey(ikey).eq.'err_fonma') then
            read(cval(ikey),*) fOnma(2)
         else if(ckey(ikey).eq.'fonea') then
            read(cval(ikey),*) fOnea(1)
         else if(ckey(ikey).eq.'err_fonea') then
            read(cval(ikey),*) fOnea(2)
         else if(ckey(ikey).eq.'foamn') then
            read(cval(ikey),*) fOamn(1)
         else if(ckey(ikey).eq.'err_foamn') then
            read(cval(ikey),*) fOamn(2)
         else if(ckey(ikey).eq.'foaen') then
            read(cval(ikey),*) fOaen(1)
         else if(ckey(ikey).eq.'err_foaen') then
            read(cval(ikey),*) fOaen(2)
         else if(ckey(ikey).eq.'foama') then
            read(cval(ikey),*) fOama(1)
         else if(ckey(ikey).eq.'err_foama') then
            read(cval(ikey),*) fOama(2)
         else if(ckey(ikey).eq.'foaea') then
            read(cval(ikey),*) fOaea(1)
         else if(ckey(ikey).eq.'err_foaea') then
            read(cval(ikey),*) fOaea(2)
         else if(ckey(ikey).eq.'oeffe') then
            read(cval(ikey),*) Oeffe(1)
         else if(ckey(ikey).eq.'err_oeffe') then
            read(cval(ikey),*) Oeffe(2)
         else if(ckey(ikey).eq.'oeffmu') then
            read(cval(ikey),*) Oeffmu(1)
         else if(ckey(ikey).eq.'err_oeffmu') then
            read(cval(ikey),*) Oeffmu(2)
         else if(ckey(ikey).eq.'ope2m') then
            read(cval(ikey),*) OPe2m(1)
         else if(ckey(ikey).eq.'err_ope2m') then
            read(cval(ikey),*) OPe2m(2)
         else if(ckey(ikey).eq.'opm2e') then
            read(cval(ikey),*) OPm2e(1)
         else if(ckey(ikey).eq.'err_opm2e') then
            read(cval(ikey),*) OPm2e(2)
         else if(ckey(ikey).eq.'fk') then
            read(cval(ikey),*) fK(1)
         else if(ckey(ikey).eq.'err_fk') then
            read(cval(ikey),*) fK(2)
         else if(ckey(ikey).eq.'fkrho') then
            read(cval(ikey),*) fKrho(1)
         else if(ckey(ikey).eq.'err_fkrho') then
            read(cval(ikey),*) fKrho(2)
         else if(ckey(ikey).eq.'fknmn') then
            read(cval(ikey),*) fKnmn(1)
         else if(ckey(ikey).eq.'err_fknmn') then
            read(cval(ikey),*) fKnmn(2)
         else if(ckey(ikey).eq.'fknen') then
            read(cval(ikey),*) fKnen(1)
         else if(ckey(ikey).eq.'err_fknen') then
            read(cval(ikey),*) fKnen(2)
         else if(ckey(ikey).eq.'fknma') then
            read(cval(ikey),*) fKnma(1)
         else if(ckey(ikey).eq.'err_fknma') then
            read(cval(ikey),*) fKnma(2)
         else if(ckey(ikey).eq.'fknea') then
            read(cval(ikey),*) fKnea(1)
         else if(ckey(ikey).eq.'err_fknea') then
            read(cval(ikey),*) fKnea(2)
         else if(ckey(ikey).eq.'fkamn') then
            read(cval(ikey),*) fKamn(1)
         else if(ckey(ikey).eq.'err_fkamn') then
            read(cval(ikey),*) fKamn(2)
         else if(ckey(ikey).eq.'fkaen') then
            read(cval(ikey),*) fKaen(1)
         else if(ckey(ikey).eq.'err_fkaen') then
            read(cval(ikey),*) fKaen(2)
         else if(ckey(ikey).eq.'fkama') then
            read(cval(ikey),*) fKama(1)
         else if(ckey(ikey).eq.'err_fkama') then
            read(cval(ikey),*) fKama(2)
         else if(ckey(ikey).eq.'fkaea') then
            read(cval(ikey),*) fKaea(1)
         else if(ckey(ikey).eq.'err_fkaea') then
            read(cval(ikey),*) fKaea(2)
         else if(ckey(ikey).eq.'keffe') then
            read(cval(ikey),*) Keffe(1)
         else if(ckey(ikey).eq.'err_keffe') then
            read(cval(ikey),*) Keffe(2)
         else if(ckey(ikey).eq.'keffmu') then
            read(cval(ikey),*) Keffmu(1)
         else if(ckey(ikey).eq.'err_keffmu') then
            read(cval(ikey),*) Keffmu(2)
         else if(ckey(ikey).eq.'kpe2m') then
            read(cval(ikey),*) KPe2m(1)
         else if(ckey(ikey).eq.'err_kpe2m') then
            read(cval(ikey),*) KPe2m(2)
         else if(ckey(ikey).eq.'kpm2e') then
            read(cval(ikey),*) KPm2e(1)
         else if(ckey(ikey).eq.'err_kpm2e') then
            read(cval(ikey),*) KPm2e(2)
         endif
      enddo

C$$$      write(lustdo,'()')
C$$$      write(lustdo,'(5x,a)') '==============================='
C$$$      write(lustdo,'(5x,a)') '  Parameters for Bases-Spring'
C$$$      write(lustdo,'(5x,a)') '==============================='
C$$$      write(lustdo,'()')

C$$$      if (ifprog.eq.ifbase) then
C$$$         if(lnblnk(bsoutf).gt.0) then
C$$$            write(lustdo,'(10x,a,a)')
C$$$     &           'BSOUTF = ',bsoutf(:lnblnk(bsoutf))
C$$$            open(lubsf, file=bsoutf,status='unknown',form='unformatted')
C$$$         endif
C$$$      elseif (ifprog.eq.ifsprg) then
C$$$         if(lnblnk(bsoutf).gt.0) then
C$$$            write(lustdo,'(10x,a,a)') 
C$$$     &           'BSOUTF = ',bsoutf(:lnblnk(bsoutf))
C$$$            open(lubsf, file=bsoutf,status='old',form='unformatted')
C$$$         else
C$$$            call msgout('bsffrd','run spring without bases data?',1)
C$$$            stop
C$$$         endif
C$$$         if(lnblnk(spevtf).gt.0) then
C$$$            write(lustdo,'(10x,a,a)') 
C$$$     &           'SPEVTF = ' ,spevtf(:lnblnk(spevtf))
C$$$            open(luevtf, file=spevtf,status='unknown',
C$$$     &           form='formatted')
C$$$         endif
C$$$      endif

C$$$      if(lnblnk(bsisrf).gt.0) then
C$$$         write(lustdo,'(10x,a,a)') 
C$$$     &        'BSISRF = ',bsisrf(:lnblnk(bsisrf))
C$$$      endif

C$$$      if(lnblnk(bshstf).gt.0) then
C$$$         write(lustdo,'(10x,a,a)') 
C$$$     &        'BSHSTF = ',bshstf(:lnblnk(bshstf))
C$$$         open(luhist, file=bshstf,status='unknown',
C$$$     &        form='formatted')
C$$$      endif

C$$$      write(lustdo,'(10x,a,i10)')   'Ncall  = ',ncall
C$$$      write(lustdo,'(10x,a,i10)')   'ITMX1  = ',itmx1
C$$$      write(lustdo,'(10x,a,i10)')   'ITMX2  = ',itmx2
C$$$      write(lustdo,'(10x,a,f10.4)') 'Acc1   = ',acc1
C$$$      write(lustdo,'(10x,a,f10.4)') 'Acc2   = ',acc2

      return
      end
