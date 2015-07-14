      Subroutine Hafc3(ca,xbin,a,x,scr,itop,ct,Nbf,Nmo,Nocc,Ncore,Nx)
      Character*16 S_nrec,S_ntot,S_nints
c        *****  Copyright (c)  1994  Polyatomics Research Inst. *****
c        *****  A. Komornicki: Polyatomics Research Institute   ***** 
c        *****  Hessian integral transformation.                *****
c        *****  Routine drives second half of transformation    *****
c        *****  This new version will generate all integrals    *****
c        *****  except:  (vv|vv)                                *****
c        *****  Calls:  Gather                                  ***** 
c        *****          Dgemm()                                 ***** 
c        *****          Xwrite                                  ***** 
c        *****          Intcoma                                 ***** 
c        *****          Xerror                                  *****
c        *****     Ca    - Molecular orbitals                   *****
c        *****     Xbin  - Input half transformed integrals     *****
c        *****     A     - Output (mo) integrals.               *****
c        *****     X     - Input integrals in square form       *****
c        *****     Scr   - Scratch array (Nbf,Nmo)              *****
c        *****     Itop  - Scatter pointer                      *****
c        *****     Ct    - Transpose of the Ca matrix.          *****
c        *****     Nbf   - Nbasis                               *****
c        *****     Nmo   - Number of vectors                    *****
c        *****     Nocc  - Number of occ vectors in trans.      *****
c        *****     Ncore - Number of frozen orbitals.           *****
      Parameter (maxb=512*96, maxa=maxb/2-1)
      Parameter (Zero=0.0E0, One=1.0e0,Tol=1.0e-9)
      Common/Bufmo/ Ncnt,jxx(maxa),xxx(maxa),xxdum
      Common/TranIO/Iuao,iusrt,Iumo,Iw
      Dimension ca(Nbf,Nmo),a(Nmo,Nmo),xbin(Nx,*),x(*),itop(*),scr(*)
      Dimension ct(Nmo,Nbf)
 1101 Format(10x,'End of Hessian integral transformation ',/,
     &       10x,'Number of records written        ', A,/,
     &       10x,'Number of non-zero integrals   ', A,/,
     &       10x,'Maximum number of integrals    ', A,/ )
 1102 Format(12x,'Actual number and maximum for each type: ',/,
     &       10x,'(oo//oo)  (vo//oo)  (vv//oo)  (vo//vo)  (vv//vo)',
     &        7x 'Total',/,
     &        9x,5(i8,2x),4x,A,/,9X,5(i8,2x),4x,A,/)
 1103 Format(10x,3(2H*-),'Error in Hessian transformation ',3(2H-*),/
     &       10x,'There appear to be no virtual orbitals - Exit ',/)
c
      Ia(kk)= kk*(kk-1)/2
c        *****  Set up pointer index.                           *****
      noooo = 0
      nvooo = 0
      nvvoo = 0
      nvovo = 0
      nvvvo = 0
      ncnt  = 0
      nrec  = 0
      nvirt = Nmo-Nocc
      Nocc1 = Nocc+1
c        *****  Set up scatter pointers:                        *****
      Nbasis = Nbf
      nddx = Nbasis*Nbasis
      ij = 0
      Do 120 i=1,Nbasis
      Do 110 j=1,Nbasis
      ij = ij+1
  110 itop(ij) = Ia(max(i,j))+min(i,j)
  120 continue
      Do 140 i=1,Nbasis
      Do 140 j=1,Nmo
  140 ct(j,i) = ca(i,j)
c        *****  Start the transformation.                       *****
      ij = 0
      Do 500 i=1,Nmo
      jend = min(Nocc,i)
      Do 500 j=1,jend
      ij = ij+1
c        *****  The following loop replaces the sort step:      *****
c        *****  Put integrals in square form:transform          *****
      Call Gather(nddx,x,xbin(1,ij),itop)
      idx  = i+Ncore
      jdx  = j+Ncore
      labij= Ishft(idx,16) +jdx
c        *****  Check for two possible cases.                   *****
      If(i.le.Nocc) then
c        *****  A(Nmo,Nmo) = ca(n,mo)*x(n,n)*ca(n,mo)           *****
         Call Dgemm('n','n',Nbasis,Nmo,Nbasis,One,x,Nbasis,ca,Nbasis,
     *                      Zero,scr,Nbasis)
         Call Dgemm('n','n',Nmo,Nmo,Nbasis,One,ct,Nmo,scr,Nbasis,
     *                      Zero,a,Nmo)
c        *****  Now fill the output buffers.                    *****
c        *****  (oc,oc//oc,oc) block first.                     *****
         Do 300 k=1,i
         lmax = k
         If(k.eq.i) lmax=j
         Do 300 l=1,lmax
         val= a(k,l)
         If(abs(val).lt.Tol) go to 300
         ncnt = ncnt+1
         noooo= noooo+1
         kdx  = k+Ncore
         ldx  = l+Ncore
         label= Ishft(labij,32)+Ishft(kdx,16)+ldx
         jxx(ncnt) = label
         xxx(ncnt) = val
         If(ncnt.ge.maxa) then
           Call Xwrite(Iumo,ncnt,maxb,ierr)
           ncnt = 0
           nrec = nrec+1
         endif
  300    continue
c        *****  Now the   (oc,oc//virt,all)  block              *****
         Do 320 k=Nocc1,Nmo
         Do 320 l=1,k
         val= a(k,l)
         If(abs(val).lt.Tol) go to 320
         ncnt = ncnt+1
         If(l.le.Nocc) then
           nvooo = nvooo+1
         else
           nvvoo = nvvoo+1
         endif
         kdx  = k+Ncore
         ldx  = l+Ncore
         label= Ishft(kdx,48)+Ishft(ldx,32)+labij
         jxx(ncnt) = label
         xxx(ncnt) = val
         If(ncnt.ge.maxa) then
           Call Xwrite(Iumo,ncnt,maxb,ierr)
           ncnt = 0
           nrec = nrec+1
         endif
  320    continue
      else
c        *****  A(nvirt,j) = ca(n,nvirt)*x(n,n)*ca(n,j)         *****
         Call Dgemm('t','n',nvirt,Nbasis,Nbasis,One,
     *                      ca(1,nocc1),Nbasis,x,Nbasis,
     *                      Zero,scr,nvirt)
         Call Dgemm('n','n',nvirt,j,Nbasis,One,scr,nvirt,
     *                      ca,Nbasis,Zero,a,Nmo)
c        *****  Now fill the output buffers.                    *****
c        *****  (virt,oc//virt,oc) block                        *****
         Do 420 l=1,j
         kend = Nmo
         If(l.eq.j) kend=i
         Do 420 k=Nocc1,kend
         val= a(k-Nocc,l)
         If(abs(val).lt.Tol) go to 420
         ncnt = ncnt+1
         nvovo= nvovo+1
         kdx  = k+Ncore
         ldx  = l+Ncore
         labkl= Ishft(kdx,16) +ldx
         If(labkl.gt.labij) then
            lab = Ishft(labkl,32) +labij
         else
            lab = Ishft(labij,32) +labkl
         endif
         jxx(ncnt) = lab
         xxx(ncnt) = val
         If(ncnt.ge.maxa) then
           Call Xwrite(Iumo,ncnt,maxb,ierr)
           ncnt = 0
           nrec = nrec+1
         endif
  420    continue
c        *****  Final part:  Get (vv/vo) + (vo/vv)              *****
c        ***** A(nvirt,nvirt) = Scr(nv,Nbasis)ca(Nbasis,nvirt)  *****
         Call Dgemm('n','n',nvirt,nvirt,Nbasis,One,scr,nvirt,
     *                      ca(1,nocc1),Nbasis,Zero,a,Nmo)
c        *****  Now fill the output buffers.                    *****
c        *****  (virt,oc//virt,virt) block                      *****
         Do 430 k=nocc1,Nmo
         Do 430 l=nocc1,k
         val= a(k-nocc,l-nocc)
         If(abs(val).lt.Tol) go to 430
         ncnt = ncnt+1
         nvvvo= nvvvo+1
         kdx  = k+Ncore
         ldx  = l+Ncore
         labkl= Ishft(kdx,16) +ldx
         If(labkl.gt.labij) then
            lab = Ishft(labkl,32) +labij
         else
            lab = Ishft(labij,32) +labkl
         endif
         jxx(ncnt) = lab
         xxx(ncnt) = val
         If(ncnt.ge.maxa) then
           Call Xwrite(Iumo,ncnt,maxb,ierr)
           ncnt = 0
           nrec = nrec+1
         endif
  430    continue
      endif
  500 continue
c        *****  Final flush of output buffers.                  *****
      Call Xwrite(Iumo,ncnt,maxb,ierr)
      nrec = nrec+1
      ntot = (nrec-1)*maxa+ncnt
      ncp  = Ia(Nocc+1)
      nvp  = Ia(nvirt+1)
      moooo= ncp*(ncp+1)/2
      mvooo= ncp*Nocc*nvirt
      mvovo= Nocc*nvirt*(Nocc*nvirt+1)/2
      mvvoo= ncp*nvp
      mvvvo= (nvirt*nocc)*nvirt*(nvirt+1)/2
      nints= moooo+ mvooo+ mvvoo+ mvovo+ mvvvo
      Call Intcoma(nrec,S_nrec,ierr)
      Call Intcoma(ntot,S_ntot,ierr)
      Call Intcoma(nints,S_nints,ierr)
      Write(Iw,1101) S_nrec,S_ntot,S_nints
      ntot = noooo+ nvooo+ nvvoo+ nvovo+ nvvvo
      Call Intcoma(ntot,S_ntot,ierr)
      Write(Iw,1102) noooo, nvooo, nvvoo, nvovo, nvvvo, S_ntot,
     &               moooo, mvooo, mvvoo, mvovo, mvvvo, S_nints
c        *****  Check for error in transformation:              *****
      If(ntot.eq.noooo) then
        Write(Iw,1103)
        Call Xerror('Hafc3')
      endif
      Return
      End
