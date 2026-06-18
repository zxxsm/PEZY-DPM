c*******************************************************************
c* DPM                                                             *
c* Short description                                               *
c*   Transport of MeV electrons in a 3D voxel geometry             *
c*                                                                 *
c* Copyright (c) 2000,2001                                         *
c* The University of Barcelona                                     *
c*                                                                 *
c* Permission to use, copy, modify, distribute and sell this       *
c* software and its documentation for any purpose is hereby granted*
c* without fee, provided that the above copyright notice appear in *
c* all copies and that both that copyright notice and this         *
c* permission notice appear in supporting documentation.  The      *
c* University of Barcelona makes no representations about the      *
c* suitability of this software for any purpose. It is provided    *
c* "as is" without express or implied warranty.                    *
c*                                                                 *
c* Copyright (c) 2000,2001                                         *
c* Polytechnical University of Catalonia                           *
c*                                                                 *
c* Permission to use, copy, modify, distribute and sell this       *
c* software and its documentation for any purpose is hereby granted*
c* without fee, provided that the above copyright notice appear in *
c* all copies and that both that copyright notice and this         *
c* permission notice appear in supporting documentation.  The      *
c* Polytechnical University of Catalonia makes no representations  *
c* about the suitability of this software for any purpose. It is   *
c* provided "as is" without express or implied warranty.           *
c*                                                                 *
c* Copyright (c) 2000,2001                                         *
c* The University of Michigan                                      *
c*                                                                 *
c* Permission to use, copy, modify, distribute and sell this       *
c* software and its documentation for any purpose is hereby granted*
c* without fee, provided that the above copyright notice appear in *
c* all copies and that both that copyright notice and this         *
c* permission notice appear in supporting documentation.  The      *
c* University of Michigan makes no representations about the       *
c* suitability of this software for any purpose. It is provided    *
c* "as is" without express or implied warranty.                    *
c*                                                                 *
c* Dependencies:                                                   *
c*   -> exports common /dpmpart/  to geometry routines             *
c*   -> imports /dpmvox/ from libgeom.f                            *
c*                                                                 *
c* Last updates                                                    *
c*   2001-02-09   JS  Fully Uncoupled Event Logic (FUEL) &         *
c*                    multi-material q(u;E) surface implemented.   *
c*   2000-11-29   SW  Modified for release.                        *
c*   2000-03-08   SW                                               *
c*   1999-02-26   JS  Created.                                     *
c*******************************************************************
c ******************************************************************
c * File: dpm.F
c * Author: Zhuxiaoxong
c * Email: zhuxiaoxiong_cn@163.com
c * Description: 与原始的dpm.f程序功能保持一致；
c *              增加了MPI的并行功能，将数据读取后和初始化后进行分发
c ******************************************************************
      program main
      implicit none
#ifdef MPI_P
      include 'mpif.h'
#endif

c     一些基本参数值
      integer*4 nxvox,nyvox,nzvox,nxyz
      parameter (nxvox=2**6,nyvox=2**6,nzvox=150)
      parameter (nxyz=nxvox*nyvox*nzvox)
      
      integer scndry,loopbk
      real*8 time0,cputim,mc2,twomc2,lapstime
      parameter (mc2=510.9991d3,twomc2=2.0d0*mc2)

c     s1-/dpmsim/
c     此结构体都是从dpm.in中读取
c     maxhis：最大历史次数
c     atime：只在控制函数command函数中调用
      integer*4 maxhis
      real*8 atime
      common /dpmsim/ atime,maxhis

c     s2-/dpmsrc/
c     此结构体都是从dpm.in中读取
c     esrc:源能量（eV）-------->Source energy (eV)
c     eabs:电子吸收能量-------->Electron absorption energy
c     eabsph:光子吸收能量-------->Photon absorption energy
c     param: 方形源束尺寸；用0表示笔形射束（厘米）-------->Square source beam size; use 0 for pencil beam (cm):
c     intype:光子/电子类型-------->Particle Type (-1 for e-, 0 for photons)
      integer intype
      real*8 esrc,eabs,eabsph,param
      common /dpmsrc/ esrc,eabs,eabsph,param,intype

c     s3-/dpmmat/
c 	在rmater中完成读取和计算
      integer maxmat,nmat
      parameter (maxmat=5)
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat

c     s4-/dpmscsr/
c     读取散射强度数据并设置插值矩阵，rstep函数读取计算
      integer nscsr
      parameter (nscsr=2**9)
      real*8 escsr,scssp,scsspa,scsspb,scsspc,scsspd,idless
      common /dpmscsr/ escsr(nscsr),scssp(nscsr),
     &                 scsspa(nscsr),scsspb(nscsr),
     &                 scsspc(nscsr),scsspd(nscsr),idless

c 	s5-/dpmscpw/ 
c     从文件中读取第一个 TMFP 并设置插值矩阵，rscpw函数读取计算
      integer nscp
      parameter (nscp=2**11)
      real*8 escp,scpsp,scpspa,scpspb,idlesc
      common /dpmscpw/ escp(nscp),scpsp(nscp),
     &       scpspa(maxmat,nscp),scpspb(maxmat,nscp),idlesc

c 	s6-/xdpmbw/
c     从文件中读取筛选参数数据并设置插值矩阵，rbw函数读取计算
      integer nbw
      parameter (nbw=2**9)
      real*8 ebw,bwsp,bwspa,bwspb,bwspc,bwspd,idlebw
      common /xdpmbw/ ebw(nbw),bwsp(nbw),
     &       bwspa(nbw),bwspb(nbw),
     &       bwspc(nbw),bwspd(nbw),idlebw

c     s7-/xdpmq2D/
c     从文件中读取 q 表面数据并设置插值矩阵，rqsurf函数读取计算
      integer nuq,neq
      parameter (nuq=2**7,neq=2**8)
      real q,le0q,idleq,iduq
      common /xdpmq2D/ q(nuq,neq),le0q,idleq,iduq

c     s8-/dpmlam/
c     从文件中读取莫勒相互作用的 lambda 并设置插值矩阵，rlammo函数读取计算
      integer nlam
      parameter (nlam=2**8)
      real*8 elam,lamsp,lamspa,lamspb,lamspc,lamspd,idlela
      common /dpmlam/ elam(nlam),lamsp(nlam),
     &  lamspa(nlam),lamspb(nlam),lamspc(nlam),lamspd(nlam),idlela

c	s9-/dpmlbr/
c     从文件中读取 lambda_radiative 数据并设置插值矩阵，rlabre函数读取计算   
      integer nlra
      parameter (nlra=2**9)
      real*8 elra,lrasp,lraspa,lraspb,idlelr
      common /dpmlbr/ elra(nlra),lrasp(nlra),
     &                lraspa(maxmat,nlra),lraspb(maxmat,nlra),
     &                idlelr

c     s10-/dpmstpw/
c     从文件中读取 StopPower_restr数据，并设置插值矩阵，rrstpw函数读取计算      
      integer nst
      parameter (nst=2**10)
      real*8 est,stsp,stspa,stspb,idlest
      common /dpmstpw/ est(nst),stsp(nst),
     &       stspa(maxmat,nst),stspb(maxmat,nst),idlest

c     s11-/xbre/
c     使用 PENELOPE 例程读取 brems 常数以进行采样，rbcon函数读取计算
      real*8 f0,bcb
      COMMON/XBRE/F0(MAXMAT,3),BCB(MAXMAT)

c     s12-/dpmlph/
c     从文件中读取光子总反向平均自由路径数据并设置插值矩阵，rlamph函数读取计算
      integer nlaph
      parameter (nlaph=2**11)
      real*8 elaph,lamph,lampha,lamphb,lamphc,lamphd,idleph
      common /dpmlph/ elaph(nlaph),lamph(nlaph),
     &       lampha(maxmat,nlaph),lamphb(maxmat,nlaph),
     &       lamphc(maxmat,nlaph),lamphd(maxmat,nlaph),idleph

c     s13-/dpmcmp/
c     从文件中读取康普顿反平均自由路径数据并设置插值矩阵，rcompt函数读取计算
      integer ncmpt
      parameter (ncmpt=2**8)
      real*8 ecmpt,compt,compta,comptb,comptc,comptd,idlecp
      common /dpmcmp/ ecmpt(ncmpt),compt(ncmpt),
     &       compta(maxmat,ncmpt),comptb(maxmat,ncmpt),
     &       comptc(maxmat,ncmpt),comptd(maxmat,ncmpt),idlecp

c     s14-/dpmpte/
c     从文件中读取光电逆平均自由路径数据并设置插值矩阵，rphote函数读取计算
c     dpm计算中没有调用后续相关函数ipheip
c     可以考虑去掉
      integer nphte
      parameter (nphte=2**12)
      real*8 ephte,phote,photea,photeb,photec,photed,idlepe
      common /dpmpte/ ephte(nphte),phote(nphte),
     &       photea(maxmat,nphte),photeb(maxmat,nphte),
     &       photec(maxmat,nphte),photed(maxmat,nphte),idlepe
      

c     s15-/dpmpap/
c     从输入文件读取体素几何图形，rpairp函数读取计算
      integer npair
      parameter (npair=2**10)
      real*8 epair,pairp,pairpa,pairpb,pairpc,pairpd,idlepp
      common /dpmpap/ epair(npair),pairp(npair),
     &       pairpa(maxmat,npair),pairpb(maxmat,npair),
     &       pairpc(maxmat,npair),pairpd(maxmat,npair),idlepp

c     s16-/dpmvox/
c     从文件中读取逆配对生产平均自由路径数据并设置插值矩阵，rvoxg函数读取计算
c     注：/dpmvox/ 与 /cgeom3/ 同时出现时，需要：
c     /dpmvox/ dx,dy,dz  ----------> /dpmvox/ ddx,ddy,ddz
      real*8 ddx,ddy,ddz,xmid,ymid,zmid,dens
      integer Unxvox,Unyvox,Unzvox,mat
      common /dpmvox/ ddx,ddy,ddz,xmid,ymid,zmid
     +          ,dens(nxyz),mat(nxyz),Unxvox,Unyvox,Unzvox

c     s17-/dpmrpt/
c     读取以 cm 为单位的矩形 RoI 并设置体素索引值，in文件并计算得到   
      integer nxini,nxfin,nyini,nyfin,nzini,nzfin
      common /dpmrpt/ nxini,nxfin,nyini,nyfin,nzini,nzfin

c     s18-/cutoff/
c     在iniion 和 inibre中计算得到
      real*8 wcion,wcbre
      common /cutoff/ wcion,wcbre

c     s19-/dpmwck/
c     在iniwck中计算得到
      integer nwck
      parameter (nwck=2**12)
      real*8 a0wck,a1wck,idlewk,wcke0
      common /dpmwck/ a0wck(nwck),a1wck(nwck),idlewk,wcke0

c     s20-/dpmsub/
c     在inisub 中计算得到
      real*8 subden,subfac,substp
      common /dpmsub/ subden,subfac,substp

c     t1-/dpmstck/
c 	记录1个原生粒子所派生的所有粒子信息
c     在以下函数中被调用
c 	putmol()：通过莫勒相互作用产生一个新的二次电子，并将其状态存储在二次堆栈中
c     putbre()：创建一个新的二次轫致辐射光子，并将其状态存储在二次堆栈中
c     photon()：传输光子，直到它逃离宇宙或其能量降到 EabsPhoton 以下为止，电子对效应中创建二级电子
c     putcom()：通过康普顿相互作用产生一个次级电子，并将其状态存储在二级堆栈中
c     putann()：为湮灭光子创建二级粒子，将其状态存储在二级堆栈中
c     scndry()：从二级堆栈中检索一个粒子，并将其状态填充到当前粒子公共堆栈中
c 	PROPERTY：private 
c  
c     sptype 派生粒子类型 -1电子/0光子
c     senerg 派生粒子能量
c     在init中初始化了nsec=0;
      integer sptype
      integer*4 maxns,nsec,sxvox,syvox,szvox,sabsvx
      real*8 swght,senerg,svx,svy,svz,sx,sy,sz
      parameter (maxns=2**7)
      common /dpmstck/ swght(maxns),senerg(maxns),
     & svx(maxns),svy(maxns),svz(maxns),sx(maxns),sy(maxns),
     & sz(maxns),
     & sptype(maxns),sxvox(maxns),syvox(maxns),szvox(maxns),
     & sabsvx(maxns),nsec

c     t2-/dpmesc/
c     记录1个原生粒子所派生的所有粒子的计数累加和；
c     在init中初始化
      integer*4 lasthi
      real*8 escore,escor2,etmp
      common /dpmesc/ escore(nxyz),escor2(nxyz),etmp(nxyz),
     &                lasthi(nxyz)

c     t3-/RSEED/
c     程序原始初始的随机种子，在init；
c     在并行版本中被 分段随机种子替代；
c     可以考虑去掉
      integer*4 iseed1,iseed2,icnt
      COMMON/RSEED/ISEED1,ISEED2,ICNT

c     t4-/sd/
c     新的分段随机种子数组，nsd表示种子分段数量
c     段数要求>=  总线程数；
c     例如在4个SC3s上计算，
      integer*4  nsd, Iseed
      parameter(nsd=100000)
      common /sd/ Iseed(2, nsd)	  

c     p1-/dpmpart/
c     记录当前计算粒子的速度位置信息
c     在source中初始化
      integer ptype
      integer*4 xvox,yvox,zvox,absvox,nhist
      real*8 energy,vx,vy,vz,x,y,z
      common /dpmpart/ energy,vx,vy,vz,x,y,z,
     &                 xvox,yvox,zvox,absvox,nhist,ptype

c     p2-/cgeom3/
c     在inigeo中生成dx/dy/dz/idx/idy/idz
c     其中dx/dy/dz与/dpmvox/中的ddx/ddy/ddz相等，
c     idx/idy/idz 记录的是dx/dy/dz的倒数
c     后续计算中，只有ivx/ivy/ivz发生变化
c     ivx/ivy/ivz 记录的是vx/vy/vz的倒数
      real*8 dx,dy,dz,idx,idy,idz,ivx,ivy,ivz
      common /cgeom3/ dx,dy,dz,idx,idy,idz,ivx,ivy,ivz   

c     p3-/dpmjmp/
c     记录电子发生作用的信息
      integer matid,event,modeel
      real*8 fuelel,fuelmo,fuelbr,burnel,burnmo,burnbr
      real*8 smax,voxden,ebefor,fuelxt
      common /dpmjmp/ fuelel,fuelmo,fuelbr,fuelxt,
     &                burnel,burnmo,burnbr,smax,voxden,
     &                ebefor,matid,event,modeel

c     p4-/cgo/
c     记录电子发生作用的时移动穿越哪个截面
      integer*4 index,dvox
      common /cgo/ index,dvox

c     DEBUG参数
      integer*4 cnt1,cnt2,cnt3;
      common /callcnt/ cnt1,cnt2,cnt3
	  
c     s21-/mpi_para/
c     MPI相关参数
      integer myrank,pro_num,ierr
      common /mpi_para/ myrank,pro_num,ierr

#ifdef MPI_P
c     MPI初始化
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,pro_num,ierr)
#else
c     没有MPI，也用到了这部分参数
      myrank=0
      pro_num=1
      ierr=0
#endif
#ifdef DEBUG
      write(*,*) "myrank=",myrank,"pro_num=",pro_num
#endif

      if (myrank.eq.0) then
            call initim
      endif
c     读取各类文件并初始化
      call init

c     读取分段随机种子文件
      call initsd

c      call sim(atime, energy, ddx,
c     &  esrc, dx, swght(1),
c     &  escore(1), a0wck(1), ecmpt(1),
c     &  nxini, elaph(1), epair(1), matden(1),
c     &  fuelel, subden, escsr(1), elam(1),
c     &  index, est(1), escp(1), elra(1),
c     &  ebw(1), q(1,1), wcion, F0(1,1),
c     &  ISEED1,Iseed(1,1),cnt1, lapstime, myrank)
     
      call sim(atime,esrc,matden(1),escsr(1),
     & escp(1),ebw(1),q(1,1),elam(1),
     & elra(1),est(1),F0(1,1),elaph(1),
     & ecmpt(1),epair(1),ddx,
     & nxini,wcion,a0wck(1),subden,
     & ISEED1,Iseed(1,1),myrank,
     & cnt1,escore(1),lapstime)
c    & energy,swght(1),dx,fuelel,index)

C      PRINT *, "q(nuq,neq) =", q(1,1), "in END"
      nhist = maxhis

      if (myrank.eq.0) then
            call report(nhist,lapstime)
      endif

      if (myrank.eq.0) then
            call endtim
      endif
#ifdef MPI_P
      call MPI_FINALIZE(ierr)
#endif
      stop
      end

      subroutine initsd
c*******************************************************************
c*    Get all seeds for all parallel threads                       *
c*                                                                 *
c*    Output:                                                      *
c*      Initialization seed for all threads                        *
c*    Author: Qinglin Wang(wangqinglin.thu@gmail.com)              *
c*******************************************************************
      implicit none
#ifdef MPI_P
      include 'mpif.h'
#endif

c     t4-/sd/
c     新的分段随机种子数组，nsd表示种子分段数量
c     段数要求>=  总线程数；
c     例如在4个SC3s上计算，
      integer*4  nsd, Iseed, i, strt
      parameter(nsd=100000, strt=30)
      common /sd/ Iseed(2, nsd)

c     s21-/mpi_para/
c     MPI相关参数
      integer myrank,pro_num,ierr
      common /mpi_para/ myrank,pro_num,ierr
      
c     *****************myrank ==0 start********************
      if (myrank.eq.0) then
      open(15, file='data/seeds/seed1.d')
      open(16, file='data/seeds/seed2.d')

      do i=1, strt
        read(15,*)
        read(16,*)
      enddo

      do i=1, nsd
        read(15,*) Iseed(1,i)
        read(16,*) Iseed(2,i)
c        write(*,'(7x,i11,4x,i11)') Iseed(1,i), Iseed(2,i)
      enddo

      close(15)
      close(16)
      endif

#ifdef MPI_P      
      call MPI_BCAST(Iseed(1,1),2*nsd,MPI_INT,0,
     & MPI_COMM_WORLD,ierr)
#endif

#ifdef DEBUG
      write(*,*) "myrank=",myrank,"Iseed(1)=",Iseed(1,1),
     & "Iseed(e)=",Iseed(2,nsd)
#endif
      end

      subroutine source
c*******************************************************************
c*    Creates a new primary particle state                         *
c*                                                                 *
c*    Output:                                                      *
c*      ptype -> -1 if it is an electron, 0 when it is a photon    *
c*      energy -> kinetic energy                                   *
c*      {vx,vy,vz} -> direction of flight                          *
c*      {x,y,z} -> position                                        *
c*      vox -> voxel#                                              *
c*    Comments:                                                    *
c*      -> It is this routine's responsability to make sure that   *
c*         all dynamic variables are assigned valid values; in     *
c*         particular, kinetic energy must be in the interval      *
c*         (Emin,Emax) defined by the material data generated with *
c*         predpm.                                                 *
c*******************************************************************
      implicit none
      integer ptype
      integer*4 xvox,yvox,zvox,absvox,nhist
      real*8 energy,vx,vy,vz,x,y,z
      common /dpmpart/ energy,vx,vy,vz,x,y,z,
     &                 xvox,yvox,zvox,absvox,nhist,ptype
      integer*4 nxvox,nyvox,nzvox,nxyz
      parameter (nxvox=2**6,nyvox=2**6,nzvox=150)
      parameter (nxyz=nxvox*nyvox*nzvox)
      integer mat
      real*8 dens,dx,dy,dz,xmid,ymid,zmid
      integer Unxvox,Unyvox,Unzvox
      common /dpmvox/ dx,dy,dz,xmid,ymid,zmid
     +          ,dens(nxyz),mat(nxyz),Unxvox,Unyvox,Unzvox
      integer intype
      real*8 esrc,eabs,eabsph,param
      common /dpmsrc/ esrc,eabs,eabsph,param,intype

      real*8 rng,zero
      parameter (zero=1.d-30)

c     *** SOURCE 1:
c       * Monoenergetic, normal, square parallel beam, side='param':
      ptype = intype
      energy = esrc
      vx = 0.d0
      vy = 0.d0
      vz = 1.d0
      x = xmid+param*(rng()-0.5d0)
      y = ymid+param*(rng()-0.5d0)
      z = +zero
      call where

cc     *** SOURCE 2: ICCR2000 benchmark
cc       * Centered square beam, side='param', point source @ z=-100 cm:
c      integer dim,ndata
c      parameter (dim=2**10)
c      real*8 e,pcum,hfside,cosmin
c      common /nrcsrc/ e(dim),pcum(dim),hfside,cosmin,ndata
c
c      character*80 buffer
c      integer i,seeki
c      real*8 costhe,sinthe,phi,zsrc,s,pi,twopi,prob,eps
c      parameter (eps=1.0d-8,zsrc=-100.0d0)
c      parameter (pi=3.1415926535897932d0,twopi=2.0d0*pi)
c
c      ptype = intype
c
c      if (nhist.lt.2) then
cc       *** Init:
c        hfside = 0.5d0*param
c        cosmin = -zsrc/dsqrt(zsrc**2+0.5d0*param**2)
c      endif
cc     *** Sample direction in a cone and reject if not in the square field:
c 10   continue
c        costhe = cosmin+(1.0d0-cosmin)*rng()
c        sinthe = dsqrt(1.0d0-costhe**2)
c        phi = twopi*rng()
c        vx = sinthe*dcos(phi)
c        vy = sinthe*dsin(phi)
c        vz = costhe
c        x = xmid
c        y = ymid
c        z = zsrc
cc       *** Find intersection with z=0 side of the universe:
c        s = -z/vz
c        x = x+s*vx
c        y = y+s*vy
c        z = +zero
c      if (dabs(x-xmid).gt.hfside.or.dabs(y-ymid).gt.hfside) goto 10
c      call where
c      if (absvox.eq.0) then
c        write(*,*)
c     &  'source:error: Particle not in universe! v{x,y,z}&vox{x,y,z}:'
c        write(*,'(3(1x,1pe10.3),3(1x,i6))') vx,vy,vz,xvox,yvox,zvox
c        write(*,*)
c     &  '              check source() for errors.'
c        stop
c      endif
c
c      if (nhist.lt.2.and.intype.eq.0) then
cc       *** Init; for photons, read spectrum from file:
c        write(*,*) ' '
c        write(*,*) ' '
c        write(*,*) 'source: Reading spectrum from spectrum.spc'
c        open(1,file='spectrum.spc')
c        read(1,*) buffer
c        read(1,*) ndata
c        if (ndata.gt.dim) then
c          write(*,*) 'source:error: too many channels in spectrum file;'
c          write(*,*) '              enlarge dim parameter in source().'
c          stop
c        endif
c        read(1,*) e(1),pcum(1)
c        if (e(1).ne.0.0d0.or.pcum(1).ne.0.0d0) then
c          write(*,*) 'source:error: First spectrum data must be 0.'
c          stop
c        endif
c        do 20 i=2,ndata
c          read(1,*) e(i),prob
cc         *** NRCC spectrum given in MeV; transform to eV:
c          e(i) = e(i)*1.0d6
c          if (prob.lt.0.0d0) then
c            write(*,*) 'source:error: Negative probability found @eV:'
c            write(*,*) e(i)
c            stop
c          endif
c          pcum(i) = pcum(i-1)+prob
c 20     continue
c        close(1)
c        write(*,*) 'source: Data successfully read: ',ndata
c        if (e(ndata).gt.esrc) then
c          write(*,*) 'source:error: Spectrum energy is too high:'
c          write(*,'(1x,2(1pe10.3,1x))') e(ndata),esrc
c          stop
c        endif
c        if (dabs(pcum(ndata)-1.0d0).gt.eps) then
c          write(*,*) 'source: Spectrum PDF not normalized:'
c          write(*,'(1pe10.3)') pcum(ndata)
c          write(*,*) '        renormalizing to 1.'
c          do 30 i=1,ndata
c            pcum(i) = pcum(i)/pcum(ndata)
c 30       continue
c        endif
c      endif
c
cc     *** Sample photon E from spectrum; E is fixed for electrons:
c      if (ptype.eq.0) then
c 40     continue
c          i = seeki(pcum,rng(),ndata)
c          energy = e(i)+(e(i+1)-e(i))*rng()
c        if (energy.lt.eabsph) goto 40
c      else
c        energy = esrc
c      endif
      end


      subroutine score(edep)
c*******************************************************************
c*    Deposites energy in the corresponding counters               *
c*                                                                 *
c*    Input:                                                       *
c*      edep -> energy being deposited (eV)                        *
c*    Comments:                                                    *
c*******************************************************************
      implicit none
      real*8 edep
      integer ptype
      integer*4 xvox,yvox,zvox,absvox,nhist
      real*8 energy,vx,vy,vz,x,y,z
      common /dpmpart/ energy,vx,vy,vz,x,y,z,
     &                 xvox,yvox,zvox,absvox,nhist,ptype
      integer*4 nxvox,nyvox,nzvox,nxyz
      parameter (nxvox=2**6,nyvox=2**6,nzvox=150)
      parameter (nxyz=nxvox*nyvox*nzvox)
      integer*4 lasthi
      real*8 escore,escor2,etmp
      common /dpmesc/ escore(nxyz),escor2(nxyz),etmp(nxyz),
     &                lasthi(nxyz)
      integer nxini,nxfin,nyini,nyfin,nzini,nzfin
      common /dpmrpt/ nxini,nxfin,nyini,nyfin,nzini,nzfin
      integer*4 cnt1,cnt2,cnt3;
      common /callcnt/ cnt1,cnt2,cnt3

      cnt1=cnt1+1
c     *** Do not score if outside the RoI:
      if (xvox.lt.nxini.or.xvox.gt.nxfin.or.yvox.lt.nyini.or.
     &    yvox.gt.nyfin.or.zvox.lt.nzini.or.zvox.gt.nzfin) return

      if (lasthi(absvox).ne.nhist) then
c       *** Transfer energy to final counters:
        escore(absvox) = escore(absvox)+etmp(absvox)
        escor2(absvox) = escor2(absvox)+etmp(absvox)**2
        etmp(absvox) = edep
        lasthi(absvox) = nhist
        cnt2=cnt2+1
      else
        etmp(absvox) = etmp(absvox)+edep
        cnt3=cnt3+1
      endif
      end


      subroutine dumpe
c*******************************************************************
c*    Dumps all tmp counters into the corresponding final counters *
c*                                                                 *
c*    Input:                                                       *
c*    Comments:                                                    *
c*      -> must be called prior to printing any report when        *
c*         exact-variance is in use                                *
c*******************************************************************
      implicit none
      integer*4 nxvox,nyvox,nzvox,nxyz
      parameter (nxvox=2**6,nyvox=2**6,nzvox=150)
      parameter (nxyz=nxvox*nyvox*nzvox)
      integer*4 lasthi
      real*8 escore,escor2,etmp
      common /dpmesc/ escore(nxyz),escor2(nxyz),etmp(nxyz),
     &                lasthi(nxyz)
      integer*4 i

      do 10 i=1,nxyz
        escore(i) = escore(i)+etmp(i)
        escor2(i) = escor2(i)+etmp(i)**2
        etmp(i) = 0.d0
 10   continue
      end


      subroutine report(n,cputim)
c*******************************************************************
c*    Writes a report of the required quantities                   *
c*                                                                 *
c*    Input:                                                       *
c*      n -> No of histories                                       *
c*      cputim -> elapsed CPU time to simulate 'n' histories       *
c*    Comments:                                                    *
c*      -> init() must be called before first call                 *
c*      -> 'cputim' does not include init and report times;        *
c*         enter 0 or neg. if not available.                       *
c*******************************************************************
      implicit none
      integer*4 n
      real*8 cputim

      integer*4 nxvox,nyvox,nzvox,nxyz
      parameter (nxvox=2**6,nyvox=2**6,nzvox=150)
      parameter (nxyz=nxvox*nyvox*nzvox)
      integer mat
      real*8 dens,dx,dy,dz,xmid,ymid,zmid
      integer Unxvox,Unyvox,Unzvox
      common /dpmvox/ dx,dy,dz,xmid,ymid,zmid
     &          ,dens(nxyz),mat(nxyz),Unxvox,Unyvox,Unzvox
      integer*4 lasthi
      real*8 escore,escor2,etmp
      common /dpmesc/ escore(nxyz),escor2(nxyz),etmp(nxyz),
     &                lasthi(nxyz)
      integer intype
      real*8 esrc,eabs,eabsph,param
      common /dpmsrc/ esrc,eabs,eabsph,param,intype
      integer nxini,nxfin,nyini,nyfin,nzini,nzfin
      common /dpmrpt/ nxini,nxfin,nyini,nyfin,nzini,nzfin
      integer*4 cnt1,cnt2,cnt3;
      common /callcnt/ cnt1,cnt2,cnt3
	  
      integer*4 i,j,k,absvox,getabs,nsig
      real*8 q,in,eff,sigma,voxvol,maxq,avesig

      write(*,*) ' '
      write(*,*) ' '
      write(*,*) 'report: Dose in target voxels'
      write(*,*)
     & '        Voxel interval {nx,ny,nz} of the Region of Interest:'
      write(*,'(4x,3(2(i5,1x),3x))')
     &  nxini,nxfin,nyini,nyfin,nzini,nzfin
      in = 1.d0/n
      voxvol = dx*dy*dz

c     *** Get the voxel at which the dose is maximum:
      maxq = 0.0d0
      do 75 k=nzini,nzfin
        do 85 j=nyini,nyfin
          do 95 i=nxini,nxfin
            absvox = getabs(i,j,k)
            if (escore(absvox).gt.maxq) maxq = escore(absvox)
 95       continue
 85     continue
 75   continue
c     *** And take half of the max as the reference value:
      maxq = maxq*0.5d0
      nsig = 0
      avesig = 0.0d0

      write(*,*) ' '
      write(*,*) ' Uncertainty (+-) at 2 sigma'
      write(*,*)
     &' x(cm)      y(cm)      z(cm)      dose(MeV/g)  +-'
      write(*,*)
     &'----------------------------------------------------'
      do 70 k=nzini,nzfin
        do 80 j=nyini,nyfin
          do 90 i=nxini,nxfin
            absvox = getabs(i,j,k)
            q = escore(absvox)*in
            sigma = escor2(absvox)*in-q**2
            sigma = dsqrt(dmax1(sigma*in,0.0d0))
c           write(*,*)escore(absvox),escor2(absvox),in,q,voxvol
            q  = q/(voxvol*dens(absvox)*1.0d6)
            sigma = sigma/(voxvol*dens(absvox)*1.0d6)
            write(*,'(1x,4(1pe10.3,1x),1pe8.1)')
     &       (i-0.5d0)*dx,(j-0.5d0)*dy,(k-0.5d0)*dz,q,2.d0*sigma
c           *** Calculate mean rel. variance for voxels with dose > Dmax/2:
            if (escore(absvox).gt.maxq) then
              avesig = avesig+(sigma/q)**2
              nsig = nsig+1
            endif
 90       continue
 80     continue
 70   continue
      avesig = dsqrt(avesig/nsig)

c     *** Call to force RNG to report:
      write(*,*) ' '
      call inirng(-666)
      write(*,*) ' '
      write(*,*) 'report: Performance statistics:'
      write(*,*) '  No of histories simulated:'
      write(*,'(4x,i12)') n
      if (cputim.gt.0.0d0) then
        write(*,*) '  CPU non-init and non-report time [t] (s):'
        write(*,'(4x,1pe10.3)') cputim
        write(*,*) '  Performance (ms/history):'
        write(*,'(4x,1pe10.3)') cputim/n*1.0d3
      else
        write(*,*) '  CPU non-init and non-report time [t] (s):'
        write(*,*) '     Sorry, time functions not available'
        write(*,*) '  Performance (ms/history):'
        write(*,*) '     Sorry, time functions not available'
      endif
      write(*,*) '  Average relative sigma from report above (%):'
      write(*,'(4x,1pe10.3)') avesig*100.0d0
      write(*,*) '  Intrinsic efficiency [N*sigma^2]^-1:'
      eff = n*avesig**2
      if (eff.gt.0.0d0) then
        write(*,'(4x,1pe10.3)') 1.0d0/eff
      else
        write(*,'(4x,1pe10.3)') 0.0d0
      endif
      write(*,*) '  Absolute efficiency [t*sigma^2]^-1 (s^-1):'
      eff = max(cputim*avesig**2,1.0d-30)
      if (cputim.gt.0.0d0) then
        write(*,'(4x,1pe10.3)') 1.0d0/eff
      else
        write(*,*) '     Sorry, time functions not available'
      endif
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) 'Have a nice day. '
      write(*,*) ' '
	  
c*******************************Debug******************************
      write(*,*) 'cnt1=', cnt1	  
      write(*,*) 'cnt2=', cnt2
      write(*,*) 'cnt3=', cnt3
      write(*,*) ' '
      end


      integer function loopbk(n)
c*******************************************************************
c*    Determines whether the simulation is done or not             *
c*                                                                 *
c*    Input:                                                       *
c*      n -> current no. of histories completed                    *
c*    Output:                                                      *
c*      -> 1 if done, 0 else                                       *
c*    Comments:                                                    *
c*******************************************************************
      implicit none
      integer*4 n

      integer*4 maxhis
      real*8 atime
      common /dpmsim/ atime,maxhis

      real*8 elaps,cputim

      loopbk = 0
      if (atime.gt.0.d0) then
        if (elaps().gt.atime) loopbk = 1
      else
        if (cputim().gt.-atime) loopbk = 1
      endif
      if (n.ge.maxhis) loopbk = 1
      end


      subroutine comand(n,nscan)
c*******************************************************************
c*    Reads and executes commands from an external file, allowing  *
c*    in-flight steering of the simulation.                        *
c*                                                                 *
c*    Input:                                                       *
c*      n -> current no. of histories already simulated.           *
c*      nscan -> command file is read every 'nscan' histories.     *
c*      -> From file command.in.  Possible codes are:              *
c*          0     -> Stop reading and get back to work             *
c*          1     -> Write report to file progress.out             *
c*          2 <N> -> Reset No of histories                         *
c*          3 <T> -> Reset simulation time, in s (if defined)      *
c*          4 <string> -> write <string> to stdout                 *
c*          5 <k> -> Following commands apply only to processor #k *
c*    Output:                                                      *
c*      -> Command is executed, command.in  reset to 0             *
c*    Comments:                                                    *
c*      -> Commands are read sequentially until command 0 is found;*
c*         all command arguments (such as <N>) must be preceded    *
c*         by a line feed (as in '2 <return> 200000')              *
c*      -> <string> must not exceed 80 chars.                      *
c*      -> Command 5 is intended for parallel processing; if not   *
c*         given, actions are assumed to refer to all processors;  *
c*         when monoprocessing, all messages are passed to the     *
c*         CPU in use.
c*      -> Note that the simulation may be stopped immediately by  *
c*         resetting <N> to 0.                                     *
c*      -> Be extremely careful before saving file command.in ;    *
c*         a syntax error is likely to cause the program to abort. *
c*******************************************************************
      implicit none
      integer*4 n,nscan

      integer*4 maxhis
      real*8 atime
      common /dpmsim/ atime,maxhis

      integer*4 nlast
      common /comlst/ nlast
      character*80 buffer
      integer*4 iscom(5),i,com,comsum,cpuid,xnp
      real*8 xatime

c     *** Scan for commands once every nscan histories:
      if (n.lt.10) nlast = 0
      if (n-nlast.lt.nscan) return
      nlast = n

c     *** Clear iscom() array:
      do 5 i=1,5
        iscom(i) = 0
 5    continue
      comsum = 0

c     *** Parse command file:
      open(1,file='data/input/command.in',err=20)
 10   continue
        read(1,*,err=20) com
        if (com.le.0) then
          goto 20
        else if (com.eq.1) then
          call progre(n)
        else if (com.eq.2) then
          read(1,*,err=20) xnp
          maxhis = xnp
        else if (com.eq.3) then
          read(1,*,err=20) xatime
          atime = xatime
        else if (com.eq.4) then
          read(1,'(a80)',err=20) buffer
        else if (com.eq.5) then
          read(1,*,err=20) cpuid
        else
          write(*,*) 'comand:WARNING: invalid command code; ignored:'
          write(*,*) com
        endif
        if (com.le.5) iscom(com)=1
        comsum = comsum+com
      goto 10
 20   continue
      close(1)

      if (comsum.eq.0) return
c     *** Reset command file:
      open(1,file='data/input/command.in')
      write(1,*) ' 0'
      write(1,*) '  '
      write(1,*) '*** Codes:  '
      write(1,*)
     & '0               -> Stop reading and get back to work'
      write(1,*)
     & '1               -> Write report to progress.out'
      write(1,*)
     & '2 <CR> <N>      -> Reset No of histories to N'
      write(1,*)
     & '3 <CR> <T>      -> Reset simulation time to T'
      write(1,*)
     & '4 <CR> <string> -> Write <string> (80 chars max) to stdout'
      write(1,*)
     & '5 <CR> <n>      -> Next commands only for processor #n'
      write(1,*) ' '
      write(1,*) '<CR> stands for carriage return'
      write(1,*) ' '
      close(1)

      write(*,*) 'comand: Command received when Nhist was:'
      write(*,'(1x,i12)') n
      write(*,*) '  Command description:'
      if (iscom(1).eq.1) then
        write(*,*) '  Report written to progress.out'
      endif
      if (iscom(2).eq.1) then
        write(*,*) '  Max Nhist reset to:'
        write(*,'(1x,i12)') maxhis
      endif
      if (iscom(3).eq.1) then
        write(*,*) '  Max time reset to:'
        write(*,'(1x,1pe10.3)') atime
      endif
      if (iscom(4).eq.1) then
        write(*,*) '*** Message from command.in  follows:'
        write(*,'(1x,a80)') buffer
      endif
      if (iscom(5).eq.1) then
        write(*,*)
     &      'comand:WARNING: Parallel processing inquiry; ignored.'
      endif
      write(*,*) '  '
      end


      subroutine progre(n)
c*******************************************************************
c*    This progress report will be printed out when requested by   *
c*    command().                                                   *
c*******************************************************************
      implicit none
      integer*4 n

      integer*4 nxvox,nyvox,nzvox,nxyz
      parameter (nxvox=2**6,nyvox=2**6,nzvox=150)
      parameter (nxyz=nxvox*nyvox*nzvox)
      integer mat
      real*8 dens,dx,dy,dz,xmid,ymid,zmid
      integer Unxvox,Unyvox,Unzvox
      common /dpmvox/ dx,dy,dz,xmid,ymid,zmid
     &          ,dens(nxyz),mat(nxyz),Unxvox,Unyvox,Unzvox
      integer*4 lasthi
      real*8 escore,escor2,etmp
      common /dpmesc/ escore(nxyz),escor2(nxyz),etmp(nxyz),
     &                lasthi(nxyz)

      integer*4 seed1,seed2,getabs,absvox,nxmid,nymid,nzmid
      real*8 q,sigma,in

      in = 1.0d0/n
      nxmid = Unxvox/2+1
      nymid = Unyvox/2+1
      nzmid = Unzvox/2+1
      absvox = getabs(nxmid,nymid,nzmid)
      q = escore(absvox)*in
      sigma = escor2(absvox)*in-q**2
      sigma = dsqrt(max(sigma*in,0.0d0))
      if (q.gt.0.0d0) then
        sigma = 100.0d0*sigma/q
      else
        sigma = 0.0d0
      endif

      open(2,file='data/output/progress.out')
      write(2,*) 'Relative uncertainty in central voxel so far (%):'
      write(2,'(1x,f7.2)') sigma
      write(2,*) 'No of histories simulated so far:'
      write(2,'(1x,i12)') n
      seed1 = -1.0d0
      call inirng(seed1)
      seed2 = -1.0d0
      call seed2n(seed2)
      write(2,*) 'Random seeds at this point:'
      write(2,'(1x,2(i12,1x))') seed1,seed2
      write(2,*) ' '
      close(2)
      end


      subroutine electr
c*******************************************************************
c*    Transports an electron until it either escapes from the      *
c*    universe or its energy drops below Eabs                      *
c*                                                                 *
c*    Input:                                                       *
c*      electron initial state                                     *
c*    Output:                                                      *
c*      deposits energy in counters and updates secondary stack    *
c*******************************************************************
      implicit none
      integer ptype
      integer*4 nxvox,nyvox,nzvox,nxyz
      parameter (nxvox=2**6,nyvox=2**6,nzvox=150)
      parameter (nxyz=nxvox*nyvox*nzvox)
      integer mat
      real*8 dens,dx,dy,dz,xmid,ymid,zmid
      integer Unxvox,Unyvox,Unzvox
      common /dpmvox/ dx,dy,dz,xmid,ymid,zmid
     +          ,dens(nxyz),mat(nxyz),Unxvox,Unyvox,Unzvox
      integer*4 xvox,yvox,zvox,absvox,nhist
      real*8 energy,vx,vy,vz,x,y,z
      common /dpmpart/ energy,vx,vy,vz,x,y,z,
     &                 xvox,yvox,zvox,absvox,nhist,ptype
      integer maxmat,nmat
      parameter (maxmat=5)
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat
      integer matid,event,modeel
      real*8 fuelel,fuelmo,fuelbr,burnel,burnmo,burnbr
      real*8 smax,voxden,ebefor,fuelxt
      common /dpmjmp/ fuelel,fuelmo,fuelbr,fuelxt,
     &                burnel,burnmo,burnbr,smax,voxden,
     &                ebefor,matid,event,modeel
      integer intype
      real*8 esrc,eabs,eabsph,param
      common /dpmsrc/ esrc,eabs,eabsph,param,intype
      real*8 subden,subfac,substp
      common /dpmsub/ subden,subfac,substp
      integer*4 maxhis
      real*8 atime
      common /dpmsim/ atime,maxhis

      real*8 costhe,de,inve
      real*8 pi,twopi,rng,lamoip,sammo,stepip,sbrems
      parameter (pi=3.1415926535897932d0,twopi=2.0d0*pi)

c     *** Loading fuel before take off:
      ebefor = energy
      fuelel = stepip(energy)
      fuelxt = fuelel*rng()
      fuelel = fuelel-fuelxt
      inve = -1.d0/energy
      fuelmo = -lamoip(inve)*dlog(rng())
      fuelbr = -dlog(rng())

c     *** Set switches to steer the flight:
      modeel = 1
      event = -1

c     *** Repeat for every flight stop:
 20   continue
        call flight

        if (event.eq.2) then
c         *** Elastic stop for refueling:
          if (modeel.eq.0) then
c           *** End of second scattering substep, no real interaction:
            ebefor = energy
            fuelel = stepip(energy)
            fuelxt = fuelel*rng()
            fuelel = fuelel-fuelxt
            modeel = 1
            event = 20
          else
c           *** Elastic scattering event:
c           *** opt-qSingleMat ON: activate next line and deact next+1:
            call samsca(ebefor,costhe)
c           * call xsamsca(mat(absvox),ebefor,costhe)
            call rotate(vx,vy,vz,costhe,twopi*rng())
            fuelel = fuelxt
            modeel = 0
          endif

        else if (event.eq.3) then
c         *** Moller stop for refueling; note that:
c           *  no Moller is simulated below energy=2*cutoff;
c           *  energy/2 > de > cutoff(moller) >= Eabs  holds, then
c           *  EnergyFinal = energy-de can never be below Eabs;
c           *  sammo() can also return de=0.
          de = energy*sammo(energy)
          if (de.gt.0.d0) then
            call putmol(de)
            energy = energy-de
          endif
          fuelmo = -lamoip(-1.d0/energy)*dlog(rng())

        else if (event.eq.4) then
c         *** Bremsstrahlung stop for refueling; note that
c           * de > cutoff(brem) >= Eabsph  OR  de=0 --see sambre():
          de = sbrems(energy,mat(absvox))
          if (de.gt.0.d0) then
            call putbre(de)
            energy = energy-de
          endif
          if (energy.lt.eabs) then
c           *** Determine if subEabs transport is needed:
            if (voxden.gt.subden) then
              call score(energy)
            else
              call subabs
            endif
            return
          endif
          fuelbr = -dlog(rng())

        else
c         *** Other event values indicate that history ended:
          return
        endif

      goto 20
      end


      subroutine flight
c*******************************************************************
c*    Transports the particle following a rectiliniar trajectory,  *
c*    taking care of interface crossings and keeping track of      *
c*    energy losses in the corresponding counters (using CSDA)     *
c*                                                                 *
c*    Input:                                                       *
c*      e- initial state                                           *
c*      Fuel variables affecting the flight                        *
c*      event -> event that stopped flight last time; in addition  *
c*                to the output codes,                             *
c*                -1 new particle                                  *
c*                20 end of elastic step                           *
c*    Output:                                                      *
c*      e- final state                                             *
c*      Remaining fuel vars                                        *
c*      event -> kind of event that causes flight to stop:         *
c*                1 run out of energy; absorbed                    *
c*                2 run out of elastic fuel                        *
c*                3 run out of Moller fuel                         *
c*                4 run out of bremsstrahlung fuel                 *
c*               99 escaped from universe                          *
c*    Comments:                                                    *
c*      -> this routine does NOT transport correctly e- below Eabs *
c*      -> internally, 'event' can also have the value 0 to        *
c*         indicate that the particle hit a voxel boundary; the    *
c*         flight is resumed without a stop.                       *
c*      -> the 'fuel vars' are contained in /dpmjmp/.              *
c*******************************************************************
      implicit none
      integer ptype
      integer*4 xvox,yvox,zvox,absvox,nhist
      real*8 energy,vx,vy,vz,x,y,z
      common /dpmpart/ energy,vx,vy,vz,x,y,z,
     &                 xvox,yvox,zvox,absvox,nhist,ptype
      integer*4 nxvox,nyvox,nzvox,nxyz
      parameter (nxvox=2**6,nyvox=2**6,nzvox=150)
      parameter (nxyz=nxvox*nyvox*nzvox)
      integer mat
      real*8 dens,dx,dy,dz,xmid,ymid,zmid
      integer Unxvox,Unyvox,Unzvox
      common /dpmvox/ dx,dy,dz,xmid,ymid,zmid
     +          ,dens(nxyz),mat(nxyz),Unxvox,Unyvox,Unzvox
      integer maxmat,nmat
      parameter (maxmat=5)
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat
      integer matid,event,modeel
      real*8 fuelel,fuelmo,fuelbr,burnel,burnmo,burnbr
      real*8 smax,voxden,ebefor,fuelxt
      common /dpmjmp/ fuelel,fuelmo,fuelbr,fuelxt,
     &                burnel,burnmo,burnbr,smax,voxden,
     &                ebefor,matid,event,modeel
      integer intype
      real*8 esrc,eabs,eabsph,param
      common /dpmsrc/ esrc,eabs,eabsph,param,intype
      real*8 subden,subfac,substp
      common /dpmsub/ subden,subfac,substp

      real*8 inters,scpwip,sback,s,rstpip,ilabip,dedx,newe,news
      real*8 de,infuel

c     *** Init according to the dynamic vars changed:
      if (event.eq.2.or.event.eq.-1) then
        call setv
        smax = inters()
      endif
      event = 0

c     *** Loop until it runs out of fuel:
 10   continue
        s = smax

c       *** Calculate fuel burn rate in the current voxel:
        matid = mat(absvox)
        voxden = dens(absvox)
        dedx = rstpip(matid,energy)*voxden
        newe = max(energy-0.5d0*dedx*s,eabs)
        burnel = scpwip(matid,newe)*voxden
        burnmo = zmass(matid)*voxden
        burnbr = voxden*ilabip(matid,-1.0d0/newe)

c       *** Burn Moller fuel:
        fuelmo = fuelmo-s*burnmo
        if (fuelmo.lt.0.d0) then
          sback = -fuelmo/burnmo
          s = s-sback
          fuelmo = 0.d0
          event = 3
        endif

c       *** Burn bremss fuel:
        fuelbr = fuelbr-s*burnbr
        if (fuelbr.lt.0.d0) then
          sback = -fuelbr/burnbr
          s = s-sback
          fuelmo = fuelmo+sback*burnmo
          fuelbr = 0.d0
          event = 4
        endif

c       *** Burn elastic fuel:
        infuel = fuelel
        fuelel = fuelel-s*burnel
        if (fuelel.lt.0.d0) then
c         *** Refine calculation of scattering 1st MFP:
          news = infuel/(scpwip(matid,energy)*voxden)
          newe = max(energy-0.5d0*dedx*news,eabs)
          news = infuel/(scpwip(matid,newe)*voxden)
          if (news.gt.s) news = s
          sback = s-news
          s = news
          fuelmo = fuelmo+sback*burnmo
          fuelbr = fuelbr+sback*burnbr
          fuelel = 0.d0
          event = 2
        endif

c       *** Accounting for continous energy loss:
        newe = max(energy-0.5d0*dedx*s,eabs)
        de = s*rstpip(matid,newe)*voxden
        energy = energy-de
        call score(de)
        if (energy.lt.eabs) then
c         *** Determine if subEabs transport is needed:
          if (voxden.gt.subden) then
            call score(energy)
          else
            call subabs()
          endif
          event = 1
          return
        endif

c       *** Move the electron:
        x = x+s*vx
        y = y+s*vy
        z = z+s*vz
        smax = smax-s

c       *** Check whether an interaction has not ocurred:
        if (event.eq.0) then
          call chvox()
          if (absvox.eq.0) then
            event = 99
            return
          endif
          smax = inters()
          goto 10
        endif

c       *** Otherwise, run out of fuel so return:
      end


      subroutine subabs
c*******************************************************************
c*    Transport of e- below the nominal Eabs until absorption      *
c*                                                                 *
c*    Input:                                                       *
c*      {x,y,z} -> take off location                               *
c*      {vx,vy,vz} -> direction of flight                          *
c*      energy -> initial kinetic energy                           *
c*    Output:                                                      *
c*      escore -> energy deposited by CSDA                         *
c*    Comments:                                                    *
c*      -> the particle flies following a straight path.           *
c*      -> uses a constant StopPow extracted from reference mat.   *
c*******************************************************************
      implicit none
      integer ptype
      integer*4 xvox,yvox,zvox,absvox,nhist
      real*8 energy,vx,vy,vz,x,y,z
      common /dpmpart/ energy,vx,vy,vz,x,y,z,
     &                 xvox,yvox,zvox,absvox,nhist,ptype
      integer*4 nxvox,nyvox,nzvox,nxyz
      parameter (nxvox=2**6,nyvox=2**6,nzvox=150)
      parameter (nxyz=nxvox*nyvox*nzvox)
      integer mat
      real*8 dens,dx,dy,dz,xmid,ymid,zmid
      integer Unxvox,Unyvox,Unzvox
      common /dpmvox/ dx,dy,dz,xmid,ymid,zmid
     +          ,dens(nxyz),mat(nxyz),Unxvox,Unyvox,Unzvox
      real*8 subden,subfac,substp
      common /dpmsub/ subden,subfac,substp

      real*8 de,s,inters,voxden

 10   continue
        voxden = dens(absvox)
c       *** Determine if further transport is needed:
        if (voxden.gt.energy*subfac) then
          call score(energy)
          return
        endif
        s = inters()
        de = substp*voxden*s
        energy = energy-de
        call score(de)
        call chvox()
        if (absvox.eq.0) then
          return
        endif
        x = x+s*vx
        y = y+s*vy
        z = z+s*vz
      goto 10
      end


      subroutine samsca(e,mu)
c*******************************************************************
c*    Samples cos(theta) according to the G&S distribution.        *
c*    Uses interpolated data for bw and the q surface and this     *
c*    latter quantity to perform a rejection procedure.            *
c*                                                                 *
c*    Input:                                                       *
c*      e -> kinetic energy in eV                                  *
c*    Output:                                                      *
c*      mu -> polar angle -cos(theta)-                             *
c*******************************************************************
      implicit none
      real*8 e,mu

      real*8 bw,onebw,u,rng,xbwip,xq2Dip,ie

      ie = -1.d0/e
      bw = xbwip(ie)
      onebw = 1.d0+bw
 10   continue
      u = rng()
      mu = (onebw-u*(onebw+bw))/(onebw-u)
      if (rng().gt.xq2Dip(u,ie)) goto 10
      end


      subroutine xsamsca(matid,e,mu)
c*******************************************************************
c*    Samples cos(theta) according to the G&S distribution.        *
c*    Uses interpolated data for bw and the q surface and this     *
c*    latter quantity to perform a rejection procedure.            *
c*    Multi-material version.                                      *
c*                                                                 *
c*    Input:                                                       *
c*      matid -> material id#                                      *
c*      e -> kinetic energy in eV                                  *
c*    Output:                                                      *
c*      mu -> polar angle -cos(theta)-                             *
c*    Comments:                                                    *
c*      -> Identical to samsca() but uses a q surface for each mat.*
c*******************************************************************
      implicit none
      integer matid
      real*8 e,mu

      real*8 bw,onebw,u,rng,bwip,q2Dip,ie

      ie = -1.d0/e
      bw = bwip(matid,ie)
      onebw = 1.d0+bw
 10   continue
      u = rng()
      mu = (onebw-u*(onebw+bw))/(onebw-u)
      if (rng().gt.q2Dip(matid,u,ie)) goto 10
      end


      subroutine putmol(elost)
c*******************************************************************
c*    Creates a new secondary electron from a Moller interaction   *
c*    and stores its state in the secondary stack                  *
c*                                                                 *
c*    Input:                                                       *
c*      elost -> kinetic energy of the secondary electron (eV)     *
c*******************************************************************
      implicit none
      real*8 elost
      integer ptype
      integer*4 xvox,yvox,zvox,absvox,nhist
      real*8 energy,vx,vy,vz,x,y,z
      common /dpmpart/ energy,vx,vy,vz,x,y,z,
     &                 xvox,yvox,zvox,absvox,nhist,ptype
      integer sptype
      integer*4 maxns,nsec,sxvox,syvox,szvox,sabsvx
      real*8 swght,senerg,svx,svy,svz,sx,sy,sz
      parameter (maxns=2**7)
      common /dpmstck/ swght(maxns),senerg(maxns),
     & svx(maxns),svy(maxns),svz(maxns),sx(maxns),sy(maxns),sz(maxns),
     & sptype(maxns),sxvox(maxns),syvox(maxns),szvox(maxns),
     & sabsvx(maxns),nsec
      real*8 rng,twopi,pi,secmo
      parameter (pi=3.1415926535897932d0,twopi=2.d0*pi)

      if (nsec.ge.maxns) then
        write(*,*)
     &   'putmol:error: Stack is full, enlarge cutoffs or maxns'
        stop
      endif
      nsec = nsec+1
      sptype(nsec) = -1
      senerg(nsec) = elost
      svx(nsec) = vx
      svy(nsec) = vy
      svz(nsec) = vz
      sx(nsec) = x
      sy(nsec) = y
      sz(nsec) = z
      sxvox(nsec) = xvox
      syvox(nsec) = yvox
      szvox(nsec) = zvox
      sabsvx(nsec) = absvox
      call rotate(svx(nsec),svy(nsec),svz(nsec),
     &            secmo(elost,energy),rng()*twopi)
      end


      subroutine putbre(elost)
c*******************************************************************
c*    Creates a new secondary bremsstrahlung photon and stores its *
c*    state in the secondary stack                                 *
c*                                                                 *
c*    Input:                                                       *
c*      elost -> energy of the secondary photon (eV)               *
c*******************************************************************
      implicit none
      real*8 elost
      integer ptype
      integer*4 xvox,yvox,zvox,absvox,nhist
      real*8 energy,vx,vy,vz,x,y,z
      common /dpmpart/ energy,vx,vy,vz,x,y,z,
     &                 xvox,yvox,zvox,absvox,nhist,ptype
      integer sptype
      integer*4 maxns,nsec,sxvox,syvox,szvox,sabsvx
      real*8 swght,senerg,svx,svy,svz,sx,sy,sz
      parameter (maxns=2**7)
      common /dpmstck/ swght(maxns),senerg(maxns),
     & svx(maxns),svy(maxns),svz(maxns),sx(maxns),sy(maxns),sz(maxns),
     & sptype(maxns),sxvox(maxns),syvox(maxns),szvox(maxns),
     & sabsvx(maxns),nsec
      real*8 rng,twopi,pi,mc2,mc2sq2,angsq2
      parameter (pi=3.1415926535897932d0,twopi=2.d0*pi)
      parameter (mc2=510.9991d3,mc2sq2=7.07106781d-1*mc2)

      nsec = nsec+1
      if (nsec.gt.maxns) then
        write(*,*)
     &   'putbre:error: Stack is full, enlarge cutoffs or maxns'
        stop
      endif
      sptype(nsec) = 0
      senerg(nsec) = elost
      svx(nsec) = vx
      svy(nsec) = vy
      svz(nsec) = vz
      sx(nsec) = x
      sy(nsec) = y
      sz(nsec) = z
      sxvox(nsec) = xvox
      syvox(nsec) = yvox
      szvox(nsec) = zvox
      sabsvx(nsec) = absvox

c     *** Polar angle set to mean value, i.e. no angular distribution used,
c       * besides, small angle approx is used for cosine:
      angsq2 = mc2sq2/(energy+mc2)
      call rotate(svx(nsec),svy(nsec),svz(nsec),
     &            1.d0-angsq2*angsq2,rng()*twopi)
      end


      subroutine photon()
c*******************************************************************
c*    Transports a photon until it either escapes from the         *
c*    universe or its energy drops below EabsPhoton                *
c*                                                                 *
c*    Input:                                                       *
c*      photon initial state                                       *
c*    Output:                                                      *
c*      deposits energy in counters and updates secondary stack    *
c*******************************************************************
      implicit none
      integer ptype
      integer*4 xvox,yvox,zvox,absvox,nhist
      real*8 energy,vx,vy,vz,x,y,z
      common /dpmpart/ energy,vx,vy,vz,x,y,z,
     &                 xvox,yvox,zvox,absvox,nhist,ptype
      integer sptype
      integer*4 maxns,nsec,sxvox,syvox,szvox,sabsvx
      real*8 swght,senerg,svx,svy,svz,sx,sy,sz
      parameter (maxns=2**7)
      common /dpmstck/ swght(maxns),senerg(maxns),
     & svx(maxns),svy(maxns),svz(maxns),sx(maxns),sy(maxns),sz(maxns),
     & sptype(maxns),sxvox(maxns),syvox(maxns),szvox(maxns),
     & sabsvx(maxns),nsec
      integer*4 nxvox,nyvox,nzvox,nxyz
      parameter (nxvox=2**6,nyvox=2**6,nzvox=150)
      parameter (nxyz=nxvox*nyvox*nzvox)
      integer mat
      real*8 dens,dx,dy,dz,xmid,ymid,zmid
      integer Unxvox,Unyvox,Unzvox
      common /dpmvox/ dx,dy,dz,xmid,ymid,zmid
     +          ,dens(nxyz),mat(nxyz),Unxvox,Unyvox,Unzvox
      integer intype
      real*8 esrc,eabs,eabsph,param
      common /dpmsrc/ esrc,eabs,eabsph,param,intype
      integer*4 lasthi
      real*8 escore,escor2,etmp
      common /dpmesc/ escore(nxyz),escor2(nxyz),etmp(nxyz),
     &                lasthi(nxyz)

      integer j
      real*8 itphip,icptip,ipapip,lamwck,lamden,prob,randno
      real*8 s,rng,efrac,costhe,de,twopi,pi,phi,lammin,mc2,twomc2
      real*8 epair(2)
      parameter (pi=3.1415926535897932d0,twopi=2.d0*pi)
      parameter (mc2=510.9991d3,twomc2=2.d0*mc2)


c     *** Loop until it either escapes or is absorbed:
 10   continue
c       *** Get lambda from the minimum lambda at the current energy:
        lammin = lamwck(energy)
        s = -lammin*dlog(rng())
        x = x+s*vx
        y = y+s*vy
        z = z+s*vz
        call where
        if (absvox.eq.0) then
          return
        endif
c       *** Apply Woodcock trick:
        lamden = lammin*dens(absvox)
        prob = 1.d0-lamden*itphip(mat(absvox),energy)
        randno = rng()

c       *** No real event; continue jumping:
        if (randno.lt.prob) goto 10

c       *** Compton:
        prob = prob+lamden*icptip(mat(absvox),energy)
        if (randno.lt.prob) then
c         *** opt-IncohScat-ON ->  Activate sincoh() and deact comsam():
c         * call sincoh(mat(absvox),energy,efrac,costhe)
          call comsam(energy,efrac,costhe)
          de = energy*(1.d0-efrac)
          phi = twopi*rng()
          if (de.lt.eabs) then
            call score(de)
          else
c           *** Create a secondary electron in the stack:
            call putcom(de,efrac,phi+pi,costhe)
          endif
          energy = energy-de
          if (energy.lt.eabsph) then
            call score(energy)
            return
          endif
          call rotate(vx,vy,vz,costhe,phi)
          goto 10
        endif

c       *** Pair production:
        prob = prob+lamden*ipapip(mat(absvox),energy)
        if (randno.lt.prob) then
          epair(1) = rng()*(energy-twomc2)
          epair(2) = energy-twomc2-epair(1)
          do 20 j=1,2
            if(epair(j).gt.eabs) then
              nsec=nsec+1
              if (nsec.gt.maxns) then
                write(*,*)
     &            'photon:error: Stack is full, enlarge maxns.'
                stop
              endif
              senerg(nsec) = epair(j)
              svx(nsec) = vx
              svy(nsec) = vy
              svz(nsec) = vz
              sx(nsec) = x
              sy(nsec) = y
              sz(nsec) = z
              sxvox(nsec) = xvox
              syvox(nsec) = yvox
              szvox(nsec) = zvox
              sabsvx(nsec) = absvox
              if(j.eq.1) then
                sptype(nsec) = -1
              else
                sptype(nsec) = +1
              endif
            else
              call score(epair(j))
              if(j.eq.2) call putann
            endif
 20       continue
          return
        endif

c     *** In any other case, photoelectric absorption occurs:
      call score(energy)
      end


      subroutine putcom(elost,efrac,phi,costhe)
c*******************************************************************
c*    Creates a secondary electron from a Compton interaction and  *
c*    stores its state in the secondary stack                      *
c*                                                                 *
c*    Input:                                                       *
c*      elost -> energy of the secondary electron being created    *
c*      efrac -> fraction of initial energy kept by 2nd photon     *
c*      phi -> athimutal angle                                     *
c*      costhe -> photon scattering angle
c*******************************************************************
      implicit none
      real*8 elost,efrac,phi,costhe
      integer ptype
      integer*4 xvox,yvox,zvox,absvox,nhist
      real*8 energy,vx,vy,vz,x,y,z
      common /dpmpart/ energy,vx,vy,vz,x,y,z,
     &                 xvox,yvox,zvox,absvox,nhist,ptype
      integer sptype
      integer*4 maxns,nsec,sxvox,syvox,szvox,sabsvx
      real*8 swght,senerg,svx,svy,svz,sx,sy,sz
      parameter (maxns=2**7)
      common /dpmstck/ swght(maxns),senerg(maxns),
     & svx(maxns),svy(maxns),svz(maxns),sx(maxns),sy(maxns),sz(maxns),
     & sptype(maxns),sxvox(maxns),syvox(maxns),szvox(maxns),
     & sabsvx(maxns),nsec
      real*8 comele

      nsec = nsec+1
      if (nsec.gt.maxns) then
        write(*,*)
     &   'putcom:error: Stack is full, enlarge cutoffs or maxns.'
        stop
      endif
      sptype(nsec) = -1
      senerg(nsec) = elost
      svx(nsec) = vx
      svy(nsec) = vy
      svz(nsec) = vz
      sx(nsec) = x
      sy(nsec) = y
      sz(nsec) = z
      sxvox(nsec) = xvox
      syvox(nsec) = yvox
      szvox(nsec) = zvox
      sabsvx(nsec) = absvox
      call rotate(svx(nsec),svy(nsec),svz(nsec),
     &               comele(energy,efrac,costhe),phi)
      end


      integer function scndry()
c*******************************************************************
c*    Retrieves a particle from the secondary stack and fills up   *
c*    the current particle common with its state                   *
c*                                                                 *
c*    Output:                                                      *
c*      -> 0 if there was no particle in the stack, 1 else         *
c*******************************************************************
      implicit none
      integer ptype
      integer*4 xvox,yvox,zvox,absvox,nhist
      real*8 energy,vx,vy,vz,x,y,z
      common /dpmpart/ energy,vx,vy,vz,x,y,z,
     &                 xvox,yvox,zvox,absvox,nhist,ptype
      integer sptype
      integer*4 maxns,nsec,sxvox,syvox,szvox,sabsvx
      real*8 swght,senerg,svx,svy,svz,sx,sy,sz
      parameter (maxns=2**7)
      common /dpmstck/ swght(maxns),senerg(maxns),
     & svx(maxns),svy(maxns),svz(maxns),sx(maxns),sy(maxns),sz(maxns),
     & sptype(maxns),sxvox(maxns),syvox(maxns),szvox(maxns),
     & sabsvx(maxns),nsec

      if (nsec.eq.0) then
        scndry = 0
        return
      endif
      scndry = 1
      energy = senerg(nsec)
      vx = svx(nsec)
      vy = svy(nsec)
      vz = svz(nsec)
      x = sx(nsec)
      y = sy(nsec)
      z = sz(nsec)
      xvox = sxvox(nsec)
      yvox = syvox(nsec)
      zvox = szvox(nsec)
      absvox = sabsvx(nsec)
      ptype = sptype(nsec)
      nsec = nsec-1
      end


      real*8 function stepip(e)
c*******************************************************************
c*    3spline interpolation for scattering strength as a function  *
c*    of kinetic energy; this quantity is related to the step      *
c*    length                                                       *
c*                                                                 *
c*    Input:                                                       *
c*      e -> kinetic energy in eV                                  *
c*    Output:                                                      *
c*      K = scattering strength = integ{ds/lambda1(s)}             *
c*    Comments:                                                    *
c*      -> init() must be called before first call                 *
c*******************************************************************
      implicit none
      real*8 e
      integer nscsr
      parameter (nscsr=2**9)
      real*8 escsr,scssp,scsspa,scsspb,scsspc,scsspd,idless
      common /dpmscsr/ escsr(nscsr),scssp(nscsr),
     &                 scsspa(nscsr),scsspb(nscsr),
     &                 scsspc(nscsr),scsspd(nscsr),idless
      integer i

      i = idless*(e-escsr(1))+1
      stepip = scsspa(i)+e*(scsspb(i)+e*(scsspc(i)+e*scsspd(i)))
      end


      real*8 function scpwip(matid,e)
c*******************************************************************
c*    Inverse 1st transport MFP --linear interpolation             *
c*                                                                 *
c*    Input:                                                       *
c*      matid -> material id#                                      *
c*      e -> kinetic energy in eV                                  *
c*    Output:                                                      *
c*      -> lambda_1^{-1} in cm^2/g                                 *
c*    Comments:                                                    *
c*      -> init() must be called before first call                 *
c*******************************************************************
      implicit none
      integer matid
      real*8 e
      integer maxmat,nscp
      parameter (maxmat=5,nscp=2**11)
      real*8 escp,scpsp,scpspa,scpspb,idlesc
      common /dpmscpw/ escp(nscp),scpsp(nscp),
     &       scpspa(maxmat,nscp),scpspb(maxmat,nscp),idlesc
      integer i

      i = idlesc*(e-escp(1))+1
      scpwip = scpspa(matid,i)+e*scpspb(matid,i)
      end


      real*8 function bwip(matid,ie)
c*******************************************************************
c*    3spline interpolation for bw as a function of mat & energy   *
c*                                                                 *
c*    Input:                                                       *
c*      matid -> material id#                                      *
c*      ie -> -1/energy in eV^-1  --kinetic energy--               *
c*    Output:                                                      *
c*      bw, broad screening parameter that gets flattest q surf    *
c*    Comments:                                                    *
c*      -> init() must be called before first call                 *
c*******************************************************************
      implicit none
      integer matid
      real*8 ie

      integer maxmat
      parameter (maxmat=5)
      integer nbw
      parameter (nbw=2**9)
      real*8 ebw,bwsp,bwspa,bwspb,bwspc,bwspd,idlebw
      common /dpmbw/ ebw(nbw),bwsp(nbw),
     &       bwspa(maxmat,nbw),bwspb(maxmat,nbw),
     &       bwspc(maxmat,nbw),bwspd(maxmat,nbw),idlebw

      integer i

      i = idlebw*(ie-ebw(1))+1
      bwip = bwspa(matid,i)+ie*(bwspb(matid,i)+ie*(bwspc(matid,i)+
     &       ie*bwspd(matid,i)))
      end


      real*8 function xbwip(ie)
c*******************************************************************
c*    3spline interpolation for bw as a function energy            *
c*                                                                 *
c*    Input:                                                       *
c*      ie -> -1/energy in eV^-1  --kinetic energy--               *
c*    Output:                                                      *
c*      bw, broad screening parameter that gets flattest q surf    *
c*    Comments:                                                    *
c*      -> init() must be called before first call                 *
c*    Comments:                                                    *
c*      -> Identical to bwip() but uses only the reference mat.    *
c*******************************************************************
      implicit none
      real*8 ie

      integer nbw
      parameter (nbw=2**9)
      real*8 ebw,bwsp,bwspa,bwspb,bwspc,bwspd,idlebw
      common /xdpmbw/ ebw(nbw),bwsp(nbw),
     &       bwspa(nbw),bwspb(nbw),
     &       bwspc(nbw),bwspd(nbw),idlebw

      integer i

      i = idlebw*(ie-ebw(1))+1
      xbwip = bwspa(i)+ie*(bwspb(i)+ie*(bwspc(i)+ie*bwspd(i)))
      end


      real*8 function q2Dip(matid,u,ie)
c*******************************************************************
c*    Linearly interpolated q(u;energy) surface                    *
c*                                                                 *
c*    Input:                                                       *
c*      matid -> material id#                                      *
c*      u -> angular variable                                      *
c*      ie -> -1/energy in eV^-1  --kinetic energy--               *
c*    Comments:                                                    *
c*      -> init() must be called before first call                 *
c*******************************************************************
      implicit none
      integer matid
      real*8 u,ie

      integer maxmat
      parameter (maxmat=5)
      integer nuq,neq
      parameter (nuq=2**7,neq=2**8)
      real q,le0q,idleq,iduq
      common /dpmq2D/ q(maxmat,nuq,neq),le0q,idleq,iduq

      integer i,j
      real*8 ru,ou,rle,ole,ouole

      ru = u*iduq
      i = ru
      ou = ru-i
      i = i+1
      rle = idleq*(ie-le0q)
      j = rle
      ole = rle-j
      j = j+1
      ouole = ou*ole
      q2Dip = q(matid,i,j)*(1.d0-ou-ole+ouole)+
     &        q(matid,i+1,j)*(ou-ouole)+q(matid,i,j+1)*(ole-ouole)+
     &        q(matid,i+1,j+1)*ouole
      end


      real*8 function xq2Dip(u,ie)
c*******************************************************************
c*    Linearly interpolated q(u;energy) surface                    *
c*                                                                 *
c*    Input:                                                       *
c*      u -> angular variable                                      *
c*      ie -> -1/energy in eV^-1  --kinetic energy--               *
c*    Comments:                                                    *
c*      -> init() must be called before first call                 *
c*    Comments:                                                    *
c*      -> Identical to q2Dip() but uses only the reference mat.   *
c*******************************************************************
      implicit none
      real*8 u,ie

      integer nuq,neq
      parameter (nuq=2**7,neq=2**8)
      real q,le0q,idleq,iduq
      common /xdpmq2D/ q(nuq,neq),le0q,idleq,iduq

      integer i,j
      real*8 ru,ou,rle,ole,ouole

      ru = u*iduq
      i = ru
      ou = ru-i
      i = i+1
      rle = idleq*(ie-le0q)
      j = rle
      ole = rle-j
      j = j+1
      ouole = ou*ole
      xq2Dip = q(i,j)*(1.d0-ou-ole+ouole)+
     &         q(i+1,j)*(ou-ouole)+q(i,j+1)*(ole-ouole)+
     &         q(i+1,j+1)*ouole
      end


      real*8 function lamoip(ie)
c*******************************************************************
c*    Moller mean free path, 3-spline interpolation                *
c*                                                                 *
c*    Input:                                                       *
c*      ie -> -1/energy in eV^-1  --kinetic energy--               *
c*    Output:                                                      *
c*      mean free path in g/cm^2                                   *
c*    Comments:                                                    *
c*      -> init() must be called before first call                 *
c*******************************************************************
      implicit none
      real*8 ie
      integer nlam
      parameter (nlam=2**8)
      real*8 elam,lamsp,lamspa,lamspb,lamspc,lamspd,idlela
      common /dpmlam/ elam(nlam),lamsp(nlam),
     &  lamspa(nlam),lamspb(nlam),lamspc(nlam),lamspd(nlam),idlela
      integer i
      real*8 inf,de
      parameter (inf=1.d30)

      de = ie-elam(1)
      if (de.gt.0.d0) then
        i = idlela*de+1
        lamoip = 1.d0/
     &    (lamspa(i)+ie*(lamspb(i)+ie*(lamspc(i)+ie*lamspd(i))))
      else
        lamoip = +inf
      endif
      end


      real*8 function ilabip(matid,ie)
c*******************************************************************
c*    Bremsstrahlung inverse MFP --linear interpolation            *
c*                                                                 *
c*    Input:                                                       *
c*      matid -> material id number                                *
c*      ie -> -1/energy in eV^-1  --kinetic energy--               *
c*    Output:                                                      *
c*      inverse mean free path in cm^2/g                           *
c*    Comments:                                                    *
c*      -> init() must be called before first call                 *
c*******************************************************************
      implicit none
      integer matid
      real*8 ie

      integer nlra,maxmat
      parameter (nlra=2**9,maxmat=5)
      real*8 elra,lrasp,lraspa,lraspb,idlelr
      common /dpmlbr/ elra(nlra),lrasp(nlra),
     &                lraspa(maxmat,nlra),lraspb(maxmat,nlra),
     &                idlelr
      integer i
      real*8 de,zero
      parameter (zero=1.0d-30)

      de = ie-elra(1)
      if (de.gt.0.d0) then
        i = idlelr*de+1
        ilabip = lraspa(matid,i)+ie*lraspb(matid,i)
      else
        ilabip = zero
      endif
      end


      real*8 function rstpip(matid,e)
c*******************************************************************
c*    Restricted stopping power --linear interpolation             *
c*                                                                 *
c*    Input:                                                       *
c*      matid -> material id#                                      *
c*      e -> energy in eV  --kinetic energy--                      *
c*    Output:                                                      *
c*      StopPow in eV*cm^2/g                                       *
c*    Comments:                                                    *
c*      -> init() must be called before first call                 *
c*******************************************************************
      implicit none
      integer matid
      real*8 e
      integer maxmat,nst
      parameter (maxmat=5,nst=2**10)
      real*8 est,stsp,stspa,stspb,idlest
      common /dpmstpw/ est(nst),stsp(nst),
     &       stspa(maxmat,nst),stspb(maxmat,nst),idlest
      integer i

      i = idlest*(e-est(1))+1
      rstpip = stspa(matid,i)+e*stspb(matid,i)
      end


      real*8 function itphip(matid,e)
c*******************************************************************
c*    Photon total inverse mean free path --3spline interpolation  *
c*                                                                 *
c*    Input:                                                       *
c*      matid -> material id#                                      *
c*      e -> kinetic energy in eV                                  *
c*    Output:                                                      *
c*      Total inverse mean free path in cm^2/g                     *
c*    Comments:                                                    *
c*      -> init() must be called before first call                 *
c*******************************************************************
      implicit none
      integer matid
      real*8 e
      integer maxmat,nlaph
      parameter (maxmat=5,nlaph=2**11)
      real*8 elaph,lamph,lampha,lamphb,lamphc,lamphd,idleph
      common /dpmlph/ elaph(nlaph),lamph(nlaph),
     &       lampha(maxmat,nlaph),lamphb(maxmat,nlaph),
     &       lamphc(maxmat,nlaph),lamphd(maxmat,nlaph),idleph
      integer i

      i = idleph*(e-elaph(1))+1
      itphip = lampha(matid,i)+e*(lamphb(matid,i)+e*(lamphc(matid,i)+
     &         e*lamphd(matid,i)))
      end


      real*8 function icptip(matid,e)
c*******************************************************************
c*    Inverse Compton mean free path --3spline interpolation       *
c*                                                                 *
c*    Input:                                                       *
c*      matid -> material id#                                      *
c*      e -> kinetic energy in eV                                  *
c*    Output:                                                      *
c*      Inverse total mean free path in cm^2/g                     *
c*    Comments:                                                    *
c*      -> init() must be called before first call                 *
c*******************************************************************
      implicit none
      integer matid
      real*8 e
      integer maxmat,ncmpt
      parameter (maxmat=5,ncmpt=2**8)
      real*8 ecmpt,compt,compta,comptb,comptc,comptd,idlecp
      common /dpmcmp/ ecmpt(ncmpt),compt(ncmpt),
     &       compta(maxmat,ncmpt),comptb(maxmat,ncmpt),
     &       comptc(maxmat,ncmpt),comptd(maxmat,ncmpt),idlecp
      integer i

      i = idlecp*(e-ecmpt(1))+1
      icptip = compta(matid,i)+e*(comptb(matid,i)+e*(comptc(matid,i)+
     &         e*comptd(matid,i)))
      end


      real*8 function ipheip(matid,e)
c*******************************************************************
c*    Inverse photoelectric mean free path --3spline interpol      *
c*                                                                 *
c*    Input:                                                       *
c*      matid -> material id#                                      *
c*      e -> kinetic energy in eV                                  *
c*    Output:                                                      *
c*      Inverse total mean free path in cm^2/g                     *
c*    Comments:                                                    *
c*      -> init() must be called before first call                 *
c*******************************************************************
      implicit none
      integer matid
      real*8 e
      integer maxmat,nphte
      parameter (maxmat=5,nphte=2**12)
      real*8 ephte,phote,photea,photeb,photec,photed,idlepe
      common /dpmpte/ ephte(nphte),phote(nphte),
     &       photea(maxmat,nphte),photeb(maxmat,nphte),
     &       photec(maxmat,nphte),photed(maxmat,nphte),idlepe
      integer i

      i = idlepe*(e-ephte(1))+1
      ipheip = photea(matid,i)+e*(photeb(matid,i)+e*(photec(matid,i)+
     &         e*photed(matid,i)))
      end


      real*8 function ipapip(matid,e)
c*******************************************************************
c*    Inverse pair production mean free path --3spline interpol    *
c*                                                                 *
c*    Input:                                                       *
c*      matid -> material id#                                      *
c*      e -> kinetic energy in eV                                  *
c*    Output:                                                      *
c*      Inverse total mean free path in cm^2/g                     *
c*    Comments:                                                    *
c*      -> init() must be called before first call                 *
c*******************************************************************
      implicit none
      integer matid
      real*8 e
      integer maxmat,npair
      parameter (maxmat=5,npair=2**10)
      real*8 epair,pairp,pairpa,pairpb,pairpc,pairpd,idlepp
      common /dpmpap/ epair(npair),pairp(npair),
     &       pairpa(maxmat,npair),pairpb(maxmat,npair),
     &       pairpc(maxmat,npair),pairpd(maxmat,npair),idlepp
      integer i

      if (e.gt.epair(1)) then
        i = idlepp*(e-epair(1))+1
        ipapip = pairpa(matid,i)+e*(pairpb(matid,i)+
     &           e*(pairpc(matid,i)+e*pairpd(matid,i)))
      else
        ipapip = 0.d0
        return
      endif
      end


      real*8 function lamwck(e)
c*******************************************************************
c*    Mean free path prepared to play the Woodcock trick           *
c*                                                                 *
c*    Input:                                                       *
c*      e -> kinetic energy in eV                                  *
c*    Output:                                                      *
c*      Minimum mean free path in cm                               *
c*    Comments:                                                    *
c*      -> iniwck() must be called before first call               *
c*******************************************************************
      implicit none
      real*8 e
      integer nwck
      parameter (nwck=2**12)
      real*8 a0wck,a1wck,idlewk,wcke0
      common /dpmwck/ a0wck(nwck),a1wck(nwck),idlewk,wcke0
      integer i

      i = (e-wcke0)*idlewk+1
      lamwck = a0wck(i)+a1wck(i)*e
      end


      subroutine init
c*******************************************************************
c*    Initializes the dpm system                                   *
c*******************************************************************
      implicit none
#ifdef MPI_P
      include 'mpif.h'
#endif

      integer*4 nxvox,nyvox,nzvox,nxyz
      parameter (nxvox=2**6,nyvox=2**6,nzvox=150)
      parameter (nxyz=nxvox*nyvox*nzvox)
      
c     s1-/dpmsim/
c     从dpm.in中读取，需要bcast
c     起始地址:atime
c     长度:1*DOUBLE + 1*INT =3*INT
      integer*4 maxhis
      real*8 atime
      common /dpmsim/ atime,maxhis

c     s2-/dpmsrc/
c     从dpm.in中读取，需要bcast
c     起始地址:esrc
c     长度:4*DOUBLE + 1*INT =9*INT
      integer intype
      real*8 esrc,eabs,eabsph,param
      common /dpmsrc/ esrc,eabs,eabsph,param,intype

c     s3-/dpmmat/
c     从rmater函数读取，需要bcast
c     起始地址:matden(1)
c     长度:3*maxmat*DOUBLE+2*DOUBLE+1*INT=(6*maxmat+5)*INT
      integer maxmat,nmat
      parameter (maxmat=5)
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat

c     s4-/dpmscsr/
c     读取散射强度数据并设置插值矩阵，rstep函数读取，需要bcast
c     起始地址:escsr(1)
c     长度:6*nscsr*DOUBLE+1*DOUBLE=(12*nscsr+2)*INT
      integer nscsr
      parameter (nscsr=2**9)
      real*8 escsr,scssp,scsspa,scsspb,scsspc,scsspd,idless
      common /dpmscsr/ escsr(nscsr),scssp(nscsr),
     &                 scsspa(nscsr),scsspb(nscsr),
     &                 scsspc(nscsr),scsspd(nscsr),idless

c     s5-/dpmscpw/
c     从文件中读取第一个 TMFP 并设置插值矩阵，rscpw函数读取，需要bcast
c     起始地址:escp(1)
c     长度:2*nscp*DOUBLE+2*nscp*maxmat*DOUBLE+1*DOUBLE=((4+4*maxmat)*nscp+2)*INT
      integer nscp
      parameter (nscp=2**11)
      real*8 escp,scpsp,scpspa,scpspb,idlesc
      common /dpmscpw/ escp(nscp),scpsp(nscp),
     &       scpspa(maxmat,nscp),scpspb(maxmat,nscp),idlesc

c     s6-/xdpmbw/
c     从文件中读取筛选参数数据并设置插值矩阵，rbw函数读取，需要bcast
c     如果调用xrbw可能需要修改
c     起始地址:ebw(1)
c     长度:6*nbw*DOUBLE+1*DOUBLE=(12*nbw+2)*INT
      integer nbw
      parameter (nbw=2**9)
      real*8 ebw,bwsp,bwspa,bwspb,bwspc,bwspd,idlebw
      common /xdpmbw/ ebw(nbw),bwsp(nbw),
     &       bwspa(nbw),bwspb(nbw),
     &       bwspc(nbw),bwspd(nbw),idlebw

c     s7-/xdpmq2D/
c     从文件中读取 q 表面数据并设置插值矩阵，rqsurf函数读取，需要bcast
c     起始地址:q(1)
c     长度:neq*nuq*FLOAT+3*FLOAT=(neq*nuq+3)*INT
      integer nuq,neq
      parameter (nuq=2**7,neq=2**8)
      real q,le0q,idleq,iduq
      common /xdpmq2D/ q(nuq,neq),le0q,idleq,iduq

c     s8-/dpmlam/
c     从文件中读取莫勒相互作用的 lambda 并设置插值矩阵，rlammo函数读取，需要bcast
c     起始地址:elam(1)
c     长度:6*nlam*DOUBLE+1*DOUBLE=(12*nlam+2)*INT
      integer nlam
      parameter (nlam=2**8)
      real*8 elam,lamsp,lamspa,lamspb,lamspc,lamspd,idlela
      common /dpmlam/ elam(nlam),lamsp(nlam),
     &  lamspa(nlam),lamspb(nlam),lamspc(nlam),lamspd(nlam),idlela

c     s9-/dpmlbr/
c     从文件中读取 lambda_radiative 数据并设置插值矩阵，rlabre函数读取，需要bcast
c     起始地址:elra(1)
c     长度:(2+2*maxmat)*nlra*DOUBLE+1*DOUBLE=((4+4*maxmat)*nlra+2)*INT
      integer nlra
      parameter (nlra=2**9)
      real*8 elra,lrasp,lraspa,lraspb,idlelr
      common /dpmlbr/ elra(nlra),lrasp(nlra),
     &                lraspa(maxmat,nlra),lraspb(maxmat,nlra),
     &                idlelr

c     s10-/dpmstpw/
c     从文件中读取 StopPower_restr数据，并设置插值矩阵，rrstpw函数读取，需要bcast
c     起始地址:est(1)
c     长度:(2+2*maxmat)*nst*DOUBLE+1*DOUBLE=((4+4*maxmat)*nst+2)*INT
      integer nst
      parameter (nst=2**10)
      real*8 est,stsp,stspa,stspb,idlest
      common /dpmstpw/ est(nst),stsp(nst),
     &       stspa(maxmat,nst),stspb(maxmat,nst),idlest

c     s11-/xbre/
c     使用 PENELOPE 例程读取 brems 常数以进行采样，rbcon函数读取，需要bcast
c     起始地址:f0(1)
c     长度:(4*maxmat)*DOUBLE=(8*maxmat)*INT
      real*8 f0,bcb
      COMMON/xbre/f0(maxmat,3),bcb(maxmat)

c     s12-/dpmlph/
c     从文件中读取光子总反向平均自由路径数据并设置插值矩阵，rlamph函数读取，需要bcast
c     起始地址:elaph(1)
c     长度:(2+4*maxmat)*nlaph*DOUBLE+1*DOUBLE=((4+8*maxmat)*nlaph+2)*INT
      integer nlaph
      parameter (nlaph=2**11)
      real*8 elaph,lamph,lampha,lamphb,lamphc,lamphd,idleph
      common /dpmlph/ elaph(nlaph),lamph(nlaph),
     &       lampha(maxmat,nlaph),lamphb(maxmat,nlaph),
     &       lamphc(maxmat,nlaph),lamphd(maxmat,nlaph),idleph

c     s13-/dpmcmp/
c     从文件中读取康普顿反平均自由路径数据并设置插值矩阵，rcompt函数读取，需要bcast
c     起始地址:ecmpt(1)
c     长度:(2+4*maxmat)*ncmpt*DOUBLE+1*DOUBLE=((4+8*maxmat)*ncmpt+2)*INT
      integer ncmpt
      parameter (ncmpt=2**8)
      real*8 ecmpt,compt,compta,comptb,comptc,comptd,idlecp
      common /dpmcmp/ ecmpt(ncmpt),compt(ncmpt),
     &       compta(maxmat,ncmpt),comptb(maxmat,ncmpt),
     &       comptc(maxmat,ncmpt),comptd(maxmat,ncmpt),idlecp

c     s14-/dpmpte/
c     从文件中读取光电逆平均自由路径数据并设置插值矩阵，rphote函数读取，需要bcast
c     起始地址:ephte(1)
c     长度:(2+4*maxmat)*nphte*DOUBLE+1*DOUBLE=((4+8*maxmat)*nphte+2)*INT
      integer nphte
      parameter (nphte=2**12)
      real*8 ephte,phote,photea,photeb,photec,photed,idlepe
      common /dpmpte/ ephte(nphte),phote(nphte),
     &       photea(maxmat,nphte),photeb(maxmat,nphte),
     &       photec(maxmat,nphte),photed(maxmat,nphte),idlepe

c     s15-/dpmpap/
c     从输入文件读取体素几何图形，rpairp函数读取，需要bcast
c     起始地址:epair(1)
c     长度:(2+4*maxmat)*npair*DOUBLE+1*DOUBLE=((4+8*maxmat)*npair+2)*INT
      integer npair
      parameter (npair=2**10)
      real*8 epair,pairp,pairpa,pairpb,pairpc,pairpd,idlepp
      common /dpmpap/ epair(npair),pairp(npair),
     &       pairpa(maxmat,npair),pairpb(maxmat,npair),
     &       pairpc(maxmat,npair),pairpd(maxmat,npair),idlepp

c     s16-/dpmvox/
c     从文件中读取逆配对生产平均自由路径数据并设置插值矩阵，rvoxg函数读取，需要bcast
c     起始地址:ddx
c     长度:(6+nxyz)*DOUBLE+(3+nxyz)*INT=((6+nxyz)*2+(3+nxyz))*INT=(15+3*nxyz)*INT
      real*8 ddx,ddy,ddz,xmid,ymid,zmid,dens
      integer Unxvox,Unyvox,Unzvox,mat
      common /dpmvox/ ddx,ddy,ddz,xmid,ymid,zmid
     +          ,dens(nxyz),mat(nxyz),Unxvox,Unyvox,Unzvox

c     s17-/dpmrpt/
c     读取以 cm 为单位的矩形 RoI 并设置体素索引值，读取文件并计算得到，需要bcast
c     起始地址:nxini
c     长度:6*INT     
      integer nxini,nxfin,nyini,nyfin,nzini,nzfin
      common /dpmrpt/ nxini,nxfin,nyini,nyfin,nzini,nzfin

c     s18-/cutoff/
c     在iniion 和 inibre中生成;不用广播
c     real*8 wcion,wcbre
c     common /cutoff/ wcion,wcbre

c     s19-/dpmwck/
c     在iniwck中生成;不用广播
c     integer nwck
c     parameter (nwck=2**12)
c     real*8 a0wck,a1wck,idlewk,wcke0
c     common /dpmwck/ a0wck(nwck),a1wck(nwck),idlewk,wcke0

c     s20-/dpmsub/
c     在inisub中生成;不用广播
c     real*8 subden,subfac,substp
c     common /dpmsub/ subden,subfac,substp

c     t1-/dpmstck/
c 	记录1个原生粒子所派生的所有粒子信息
c     在init中初始化了nsec=0;可以自己归零生成;不用广播
      integer sptype
      integer*4 maxns,nsec,sxvox,syvox,szvox,sabsvx
      real*8 swght,senerg,svx,svy,svz,sx,sy,sz
      parameter (maxns=2**7)
      common /dpmstck/ swght(maxns),senerg(maxns),
     & svx(maxns),svy(maxns),svz(maxns),sx(maxns),sy(maxns),sz(maxns),
     & sptype(maxns),sxvox(maxns),syvox(maxns),szvox(maxns),
     & sabsvx(maxns),nsec

c     t2-/dpmesc/
c     计数器，每个进程可以自己归零初始化;不用广播
      integer*4 lasthi
      real*8 escore,escor2,etmp
      common /dpmesc/ escore(nxyz),escor2(nxyz),etmp(nxyz),
     &                lasthi(nxyz)

c     t3-/RSEED/
c     程序原始初始的随机种子，需要bcast，用于验证；
c     在并行版本中被 分段随机种子替代；
c     起始地址:iseed1
c     长度:3*INT 
      integer*4 iseed1,iseed2,icnt
      COMMON/RSEED/ISEED1,ISEED2,ICNT

c     p2-/cgeom3/
c     在inigeo中生成;不用广播
c     其中dx/dy/dz与/dpmvox/中的ddx/ddy/ddz相等，
c     idx/idy/idz 记录的是dx/dy/dz的倒数
c     后续计算中，只有ivx/ivy/ivz发生变化
c     ivx/ivy/ivz 记录的是vx/vy/vz的倒数
c     real*8 dx,dy,dz,idx,idy,idz,ivx,ivy,ivz
c     common /cgeom3/ dx,dy,dz,idx,idy,idz,ivx,ivy,ivz

c     一些临时变量
      character*100 buffer
      character*40 fname
      character*20 prefix
      integer len
      integer*4 i,bytes,seed1,seed2
      real*8 kbytes,emin,emax,eminph,refz,refz2,refmas,refden
      real*8 stdspw,xini,xfin,yini,yfin,zini,zfin
      parameter (stdspw=2.0d6)

c     s21-/mpi_para/
c     MPI相关参数
      integer myrank,pro_num,ierr
      common /mpi_para/ myrank,pro_num,ierr
      
      kbytes = 0.0d0
c     *****************myrank ==0 start********************
      if (myrank.eq.0) then

      write(*,*) ' '
      write(*,*) 'init: Reading simulation setup from STDIN;'
      write(*,*) '      information from this stream follows:'
c     s1-生成结构体 /dpmsim/
      read(*,'(a100)') buffer
      write(*,*) buffer
      read(*,*) maxhis
      write(*,'(1x,i10)') maxhis

      read(*,'(a100)') buffer
      write(*,*) buffer
      read(*,*) atime
      write(*,'(1x,1pe10.3)') atime
c     s2-生成结构体 /dpmsrc/
      read(*,'(a100)') buffer
      write(*,*) buffer
      read(*,*) intype
      write(*,'(1x,i2)') intype

      read(*,'(a100)') buffer
      write(*,*) buffer
      read(*,*) esrc
      write(*,'(1x,1pe10.3)') esrc

      read(*,'(a100)') buffer
      write(*,*) buffer
      read(*,*) param
      write(*,'(1x,1pe10.3)') param

      read(*,'(a100)') buffer
      write(*,*) buffer
      read(*,*) eabs
      write(*,'(1x,1pe10.3)') eabs

      read(*,'(a100)') buffer
      write(*,*) buffer
      read(*,*) eabsph
      write(*,'(1x,1pe10.3)') eabsph

c     *** Read the prefix with the name of the data files, truncate it:
      read(*,'(a100)') buffer
      write(*,*) buffer
      call getnam(5,prefix,len)
      write(*,*) prefix(1:len)

c     *** Read material data
      fname= prefix(1:len)//'.matter'
c     s3-生成结构体 /dpmmat/
      call rmater(fname,emin,eminph,emax,
     &            refz,refz2,refmas,refden)
      write(*,*) ' '
      write(*,*) ' '
      if (esrc.ge.emax.or.eabs.lt.emin.or.eabsph.lt.eminph) then
        write(*,*) 'init:error: Esrc or Eabs out of range.'
        stop
      endif
c     *** These restrictions save some 'if's in the e- transport routine:
      if (wcc.lt.eabs.or.wcb.lt.eabsph) then
        write(*,*) 'init:error: Cutoffs cannot be less than Eabs'
        stop
      endif
      fname= prefix(1:len)//'.step'
c     s4-生成结构体 /dpmscsr/
      call rstep(fname,bytes)  
      kbytes = kbytes+bytes*1.d-3
      fname= prefix(1:len)//'.scpw'
c     s5-生成结构体 /dpmscpw/
      call rscpw(fname,bytes)
      kbytes = kbytes+bytes*1.d-3
      fname= prefix(1:len)//'.bw'
c     *** opt-qSingleMat ON: activate next line and deact next+1:
c     s6-生成结构体 /xdpmbw/
      call rbw(fname,bytes)
c     * call xrbw(fname,bytes)
      kbytes = kbytes+bytes*1.d-3
      fname= prefix(1:len)//'.q'
c     *** opt-qSingleMat ON: activate next line and deact next+1:
c     s7-生成结构体 /xdpmq2D/
      call rqsurf(fname,bytes)
c     * call xrqsurf(fname,bytes)
      kbytes = kbytes+bytes*1.d-3
      fname= prefix(1:len)//'.lammo'
c     s8-生成结构体 /dpmlam/
      call rlammo(fname,bytes)
      kbytes = kbytes+bytes*1.d-3
      fname= prefix(1:len)//'.lambre'
c     s9-生成结构体 /dpmlbr/
      call rlabre(fname,bytes)
      kbytes = kbytes+bytes*1.d-3
      fname= prefix(1:len)//'.rstpw'
c     s10-生成结构体 /dpmstpw/
      call rrstpw(fname,bytes)
      kbytes = kbytes+bytes*1.d-3
      fname= prefix(1:len)//'.brecon'
c     s11-生成结构体 /xbre/
      call rbcon(fname,bytes)
      kbytes = kbytes+bytes*1.d-3
      fname= prefix(1:len)//'.lamph'
c     s12-生成结构体 /dpmlph/
      call rlamph(fname,bytes)
      kbytes = kbytes+bytes*1.d-3
      fname= prefix(1:len)//'.compt'
c     s13-生成结构体 /dpmcmp/
      call rcompt(fname,bytes)
      kbytes = kbytes+bytes*1.d-3
      fname= prefix(1:len)//'.phote'
c     s14-生成结构体 /dpmpte/
      call rphote(fname,bytes)
      kbytes = kbytes+bytes*1.d-3
      fname= prefix(1:len)//'.pairp'
c     s15-生成结构体 /dpmpap/
      call rpairp(fname,bytes)
      kbytes = kbytes+bytes*1.d-3

c     *** Read the geometry file name, and then the geometry
      read(*,'(a100)') buffer
      read(*,'(a100)') buffer
c     16-生成结构体 /dpmvox/
      call rvoxg(buffer)

c     *** Read the rectangular RoI in cm and set voxel index values:
c     s17-生成结构体 /dpmrpt/
      read(*,'(a100)') buffer
      read(*,*) xini,xfin,yini,yfin,zini,zfin
      nxini = xini/ddx+1
      if (nxini.lt.1) nxini = 1
      if (nxini.gt.Unxvox) nxini = Unxvox
      nxfin = xfin/ddx+1
      if (nxfin.lt.1) nxfin = 1
      if (nxfin.gt.Unxvox) nxfin = Unxvox
      nyini = yini/ddy+1
      if (nyini.lt.1) nyini = 1
      if (nyini.gt.Unyvox) nyini = Unyvox
      nyfin = yfin/ddy+1
      if (nyfin.lt.1) nyfin = 1
      if (nyfin.gt.Unyvox) nyfin = Unyvox
      nzini = zini/ddz+1
      if (nzini.lt.1) nzini = 1
      if (nzini.gt.Unzvox) nzini = Unzvox
      nzfin = zfin/ddz+1
      if (nzfin.lt.1) nzfin = 1
      if (nzfin.gt.Unzvox) nzfin = Unzvox

c     *** Read rng() seed and load:
c     *** 原始的seed 设置为默认值0 0***************
      write(*,*) ' '
      write(*,*) ' '
      read(*,'(a100)') buffer
      write(*,*) buffer
      read(*,*) seed1,seed2
      write(*,'(1x,2(i10,1x))') seed1,seed2
      call inirng(seed1)
      call seed2n(seed2)
c     *** One more call to force RNG to report:
      call inirng(-666)
c     *****************myrank ==0 end*********************
      endif

#ifdef MPI_P
c     *****************mpi_bcast_start********************
c     s1-common /dpmsim/
      call MPI_BCAST(atime,3,MPI_INT,0,MPI_COMM_WORLD,ierr)

c     s2-common /dpmsrc/
      call MPI_BCAST(esrc,9,MPI_INT,0,MPI_COMM_WORLD,ierr)

c     s3-common /dpmmat/
      call MPI_BCAST(matden(1),maxmat*6+5,MPI_INT,0,
     & MPI_COMM_WORLD,ierr)
      
c     s4-common /dpmscsr/
      call MPI_BCAST(escsr(1),12*nscsr+2,MPI_INT,0,
     & MPI_COMM_WORLD,ierr)
      
c     s5-common /dpmscpw/
      call MPI_BCAST(escp(1),(4+4*maxmat)*nscp+2,MPI_INT,0,
     & MPI_COMM_WORLD,ierr)
      
c     s6-common /xdpmbw/
      call MPI_BCAST(ebw(1),12*nbw+2,MPI_INT,0,
     & MPI_COMM_WORLD,ierr)
      
c     s7-common /xdpmq2D/
      call MPI_BCAST(q(1,1),neq*nuq+3,MPI_INT,0,
     & MPI_COMM_WORLD,ierr)
      
c     s8-common /dpmlam/
      call MPI_BCAST(elam(1),12*nlam+2,MPI_INT,0,
     & MPI_COMM_WORLD,ierr)

c     s9-common /dpmlbr/
      call MPI_BCAST(elra(1),(4+4*maxmat)*nlra+2,MPI_INT,0,
     & MPI_COMM_WORLD,ierr)
      
c     s10-common /dpmstpw/
      call MPI_BCAST(est(1),(4+4*maxmat)*nst+2,MPI_INT,0,
     & MPI_COMM_WORLD,ierr)
      
c     s11-common /xbre/
      call MPI_BCAST(f0(1,1),8*maxmat,MPI_INT,0,
     & MPI_COMM_WORLD,ierr)

c     s12-common /dpmlph/
      call MPI_BCAST(elaph(1),(4+8*maxmat)*nlaph+2,MPI_INT,0,
     & MPI_COMM_WORLD,ierr)
      
c     s13-common /dpmcmp/
      call MPI_BCAST(ecmpt(1),(4+8*maxmat)*ncmpt+2,MPI_INT,0,
     & MPI_COMM_WORLD,ierr)
      
c     s14-common /dpmpte/
      call MPI_BCAST(ephte(1),(4+8*maxmat)*nphte+2,MPI_INT,0,
     & MPI_COMM_WORLD,ierr)

c     s15-common /dpmpap/
      call MPI_BCAST(epair(1),(4+8*maxmat)*npair+2,MPI_INT,0,
     & MPI_COMM_WORLD,ierr)

c     s16-common /dpmvox/
      call MPI_BCAST(ddx,15+3*nxyz,MPI_INT,0,MPI_COMM_WORLD,ierr)

c     s17-common /dpmrpt/
      call MPI_BCAST(nxini,6,MPI_INT,0,MPI_COMM_WORLD,ierr)

c     t3-common /RSEED/
      call MPI_BCAST(iseed1,3,MPI_INT,0,MPI_COMM_WORLD,ierr)
#endif
#ifdef DEBUG
      write(*,*) "myrank=",myrank,"atime=",atime,"maxhis=",maxhis
      write(*,*) "myrank=",myrank,"esrc=",esrc,"intype=",intype
      write(*,*) "myrank=",myrank,"matden=",matden(1),"nmat=",nmat
      write(*,*) "myrank=",myrank,"escsr=",escsr(1),"idless=",idless
      write(*,*) "myrank=",myrank,"escp=",escp(1),"idlesc=",idlesc
      write(*,*) "myrank=",myrank,"ebw=",ebw(1),"idlebw=",idlebw
      write(*,*) "myrank=",myrank,"q=",q(1,1),"iduq=",iduq
      write(*,*) "myrank=",myrank,"elam=",elam(1),"idlela=",idlela
      write(*,*) "myrank=",myrank,"elra=",elra(1),"idlelr=",idlelr
      write(*,*) "myrank=",myrank,"est=",est(1),"idlest=",idlest
      write(*,*) "myrank=",myrank,"f0=",f0(1,1),"bcb=",bcb(maxmat)
      write(*,*) "myrank=",myrank,"elaph=",elaph(1),"idleph=",idleph
      write(*,*) "myrank=",myrank,"ecmpt=",ecmpt(1),"idlecp=",idlecp
      write(*,*) "myrank=",myrank,"ephte=",ephte(1),"idlepe=",idlepe
      write(*,*) "myrank=",myrank,"epair=",epair(1),"idlepp=",idlepp
      write(*,*) "myrank=",myrank,"ddx=",ddx,"Unzvox=",Unzvox
      write(*,*) "myrank=",myrank,"nxini=",nxini,"nzfin=",nzfin
      write(*,*) "myrank=",myrank,"iseed1=",iseed1,"iseed2=",iseed2
#endif
c     *****************mpi_bcast_end*********************
#ifdef MPI_P
c     *** refz/refz2/refmas/refden ***
      call MPI_BCAST(refz,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(refz2,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(refmas,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(refden,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
c     *****  eminph/emax/bytes  ******
      call MPI_BCAST(eminph,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(emax,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(bytes,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
#endif
c     s18-生成结构体/cutoff/ /mater/ /factor/ /cbre/
c     只有/cutoff/参与了dpm的后续计算
      call iniion(refz,refmas,wcc)
      call inibre(refz2,refmas,wcb)

c     p2-生成结构体/cgeom3/
      call inigeo(ddx,ddy,ddz,nxvox,nyvox,nzvox)

c     *** iniwck must be called after reading esrc & eabsph:
c     s19-生成结构体/dpmwck/
      call iniwck(eminph,emax,bytes)

c     *** inisub called after eabs has been read:
c     s20-生成结构体/dpmsub/
      call inisub(refden)

c     *****************myrank ==0 start********************
      if (myrank.eq.0) then
      kbytes = kbytes+bytes*1.d-3
      write(*,*) ' '
      write(*,*) ' '
      write(*,*)
     &  'init: kB by interpolation arrays and secondary stack:'
      kbytes = kbytes+(8*8+1*4+5*2)*maxns*1.d-3
      write(*,'(1x,1pe8.1)') kbytes
      write(*,*) 'init: kB by 3D energy counters:'
      kbytes = (3*8+1*4)*nxyz*1.d-3
      write(*,'(1x,1pe8.1)') kbytes
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) 'init: Done.'
c     *****************myrank ==0 end*********************
      endif

c     *** Cleaning counters:
c     初始化t1-/dpmstck/中的nsec
c     初始化t2-/dpmesc/中escore(nxyz),escor2(nxyz),etmp(nxyz),lasthi(nxyz)
      nsec = 0
      do 10 i=1,nxyz
        escore(i) = 0.d0
        escor2(i) = 0.d0
        etmp(i) = 0.d0
        lasthi(i) = 0
 10   continue
      end


      subroutine rmater(fname,emin,eminph,emax,
     &                  refz,refz2,refmas,refden)
c*******************************************************************
c*    Reads material data from file                                *
c*                                                                 *
c*    Output:                                                      *
c*      fname -> input file name                                   *
c*      [Emin,Eminph,Emax] -> interval where data will be gen (eV) *
c*      refz -> total atomic no of the reference material          *
c*      refz2 -> atomic no^2 of the reference material             *
c*      refmas -> atomic weight of the reference material          *
c*      refden -> density of the reference material (g/cm^3)       *
c*******************************************************************
      implicit none
      character*40 fname
      real*8 emin,eminph,emax,refz,refz2,refmas,refden
      integer maxmat
      parameter (maxmat=5)
      integer nmat
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat

      character*100 buffer
      integer i,j,nelem,len
      real*8 atno,atno2,mass,shigh,slow,ecross
      character*40 mname

      write(*,*) ' '
      write(*,*) ' '
      write(*,*) 'rmater: Reading ', fname
      write(*,*) '        information from this file follows:'
      open(1,file=fname)
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      write(*,*) buffer
      read(1,'(a100)') buffer
      write(*,*) buffer
      read(1,*) eminph,emin,emax
      write(*,'(1x,3(1pe10.3,1x))') eminph,emin,emax
      read(1,'(a100)') buffer
      write(*,*) buffer
      read(1,*) wcc,wcb
      write(*,'(1x,2(1pe10.3,1x))') wcc,wcb
      read(1,'(a100)') buffer
      write(*,*) buffer
      read(1,*) shigh,slow,ecross
      write(*,'(1x,3(1pe10.3,1x))') shigh,slow,ecross

      read(1,'(a100)') buffer
      write(*,*) buffer
      read(1,*) nmat
      write(*,'(1x,1(i3,1x))') nmat
      if (nmat.gt.maxmat) then
        write(*,*) 'rmater:error: Too many materials.'
        stop
      endif

      do 10 i=1,nmat
        read(1,'(a100)') buffer
        write(*,*) buffer
c       *** Read name of material, remove trailing blanks:
        call getna2(1,mname,len)
        read(1,'(a100)') buffer
        write(*,*) 'MATERIAL: ',mname(1:len)
        read(1,*) matden(i)
        write(*,'(1x,1pe14.7)') matden(i)
        read(1,'(a100)') buffer
        write(*,*) buffer
        read(1,*) nelem
        write(*,'(1x,i3)') nelem
        do 20 j=1,nelem
          read(1,'(a100)') buffer
          write(*,*) buffer
 20     continue
        read(1,'(a100)') buffer
        write(*,*) buffer
        read(1,*) atno,atno2
        write(*,'(1x,2(1pe14.7,1x))') atno,atno2
        read(1,'(a100)') buffer
        write(*,*) buffer
        read(1,*) mass
        write(*,'(1x,1(1pe14.7,1x))') mass
c       *** Take mat# 1 as the reference material:
        if (i.eq.1) then
          refz = atno
          refz2 = atno2
          refmas = mass
        endif
        read(1,'(a100)') buffer
        write(*,*) buffer
        read(1,*) zmass(i),z2mass(i)
        write(*,'(1x,2(1pe14.7,1x))') zmass(i),z2mass(i)
c       *** opt-IncohScat-ON -> activate inscati():
c       * call inscati(mname,len,i)
 10   continue
      close(1)

c     *** Get reference density:
      refden = matden(1)

      write(*,*) ' '
      write(*,*) 'rmater: Done.'
      end


      subroutine rstep(fname,bytes)
c*******************************************************************
c*    Reads scattering strength data and sets up interpolation     *
c*    matrices                                                     *
c*                                                                 *
c*    Input:                                                       *
c*      fname -> input file name                                   *
c*    Output:                                                      *
c*      bytes -> memory filled up by interpolation arrays          *
c*******************************************************************
      implicit none
      character*40 fname
      integer*4 bytes

      integer nscsr
      parameter (nscsr=2**9)
      real*8 escsr,scssp,scsspa,scsspb,scsspc,scsspd,idless
      common /dpmscsr/ escsr(nscsr),scssp(nscsr),
     &                 scsspa(nscsr),scsspb(nscsr),
     &                 scsspc(nscsr),scsspd(nscsr),idless

      character*100 buffer
      integer i,ndata
      real*8 shigh,slow,ecross

      write(*,*) 'rstep: Reading ',fname
      open(1,file=fname)
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      write(*,*) 'rstep: ',buffer
      read(1,*) shigh,slow,ecross
      write(*,'(1x,3(1pe10.3,1x))') shigh,slow,ecross
      read(1,'(a100)') buffer
      read(1,*) ndata
      if (ndata.ne.nscsr) then
        write(*,*) 'rstep:error: Array dim do not match:'
        write(*,'(1x,2(i5,1x))') ndata,nscsr
        stop
      endif
      read(1,'(a100)') buffer
c     *** Prepare interpolation:
      do 20 i=1,nscsr
        read(1,*) escsr(i),scssp(i)
 20   continue
      close(1)
      idless = (nscsr-1)/(escsr(nscsr)-escsr(1))
      call spline(escsr,scssp,scsspa,scsspb,scsspc,scsspd,
     &  0.d0,0.d0,nscsr)
      bytes = 8*(6*nscsr)
      end


      subroutine rscpw(fname,bytes)
c*******************************************************************
c*    Reads 1st TMFP from file and sets up interpolation matrices  *
c*                                                                 *
c*    Input:                                                       *
c*      fname -> input file name                                   *
c*    Output:                                                      *
c*      bytes -> memory filled up by interpolation arrays          *
c*******************************************************************
      implicit none
      character*40 fname
      integer*4 bytes

      integer maxmat,nscp
      parameter (maxmat=5,nscp=2**11)
      real*8 escp,scpsp,scpspa,scpspb,idlesc
      common /dpmscpw/ escp(nscp),scpsp(nscp),
     &       scpspa(maxmat,nscp),scpspb(maxmat,nscp),idlesc
      integer nmat
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat

      character*100 buffer
      integer i,j,ndata
      real*8 dummya(nscp),dummyb(nscp)

      write(*,*) 'rscpw: Reading ', fname
      open(1,file=fname)
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer

      do 20 j=1,nmat
        read(1,'(a100)') buffer
        read(1,*) ndata
        if (ndata.ne.nscp) then
          write(*,*) 'rscpw:error: Array dim do not match:'
          write(*,'(1x,2(i5,1x))') ndata,nscp
          stop
        endif
        read(1,'(a100)') buffer
c       *** Prepare interpolation:
        do 10 i=1,nscp
          read(1,*) escp(i),scpsp(i)
 10     continue
        read(1,'(a100)') buffer
        do 30 i=1,nscp-1
          dummyb(i) = (scpsp(i+1)-scpsp(i))/(escp(i+1)-escp(i))
          dummya(i) = scpsp(i)-escp(i)*dummyb(i)
  30    continue
c       *** Loading dummy arrays into multimaterial sp matrices:
        do 15 i=1,nscp
          scpspa(j,i) = dummya(i)
          scpspb(j,i) = dummyb(i)
 15     continue
 20   continue
      close(1)
      idlesc = (nscp-1)/(escp(nscp)-escp(1))
      bytes = 8*(maxmat*nscp*2+nscp*2)
      end


      subroutine rbw(fname,bytes)
c*******************************************************************
c*    Reads screening parameter data from file and sets up         *
c*    interpolation matrices                                       *
c*                                                                 *
c*    Input:                                                       *
c*      fname -> input file name                                   *
c*    Output:                                                      *
c*      bytes -> memory filled up by interpolation arrays          *
c*******************************************************************
      implicit none
      character*40 fname
      integer*4 bytes

      integer maxmat
      parameter (maxmat=5)
      integer nmat
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat
      integer nbw
      parameter (nbw=2**9)
      real*8 ebw,bwsp,bwspa,bwspb,bwspc,bwspd,idlebw
      common /xdpmbw/ ebw(nbw),bwsp(nbw),
     &       bwspa(nbw),bwspb(nbw),
     &       bwspc(nbw),bwspd(nbw),idlebw

      character*100 buffer
      integer i,ndata
      real*8 dummya(nbw),dummyb(nbw),dummyc(nbw),dummyd(nbw)

      write(*,*) 'rbw: Reading ', fname
      open(1,file=fname)
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      read(1,*) ndata
      read(1,'(a100)') buffer
      if (ndata.ne.nbw) then
        write(*,*) 'rbw:error: Array dim do not match:'
        write(*,'(1x,2(i5,1x))') ndata,nbw
        stop
      endif
c     *** Prepare interpolation:
      do 10 i=1,nbw
        read(1,*) ebw(i),bwsp(i)
        ebw(i) = -1.d0/ebw(i)
 10   continue
      read(1,'(a100)') buffer
      call spline(ebw,bwsp,dummya,dummyb,dummyc,dummyd,
     &              0.d0,0.d0,nbw)
      do 15 i=1,nbw
        bwspa(i) = dummya(i)
        bwspb(i) = dummyb(i)
        bwspc(i) = dummyc(i)
        bwspd(i) = dummyd(i)
 15   continue

      close(1)
      idlebw = (nbw-1)/(ebw(nbw)-ebw(1))
      bytes = 8*(nbw*6)
      end


      subroutine xrbw(fname,bytes)
c*******************************************************************
c*    Reads screening parameter data from file and sets up         *
c*    interpolation matrices. Multi-material version.              *
c*                                                                 *
c*    Input:                                                       *
c*      fname -> input file name                                   *
c*    Output:                                                      *
c*      bytes -> memory filled up by interpolation arrays          *
c*    Comments:                                                    *
c*      -> This routine is identical to rbw() except that it loads *
c*         different data for different materials.                 *
c*******************************************************************
      implicit none
      character*40 fname
      integer*4 bytes

      integer maxmat
      parameter (maxmat=5)
      integer nmat
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat
      integer nbw
      parameter (nbw=2**9)
      real*8 ebw,bwsp,bwspa,bwspb,bwspc,bwspd,idlebw
      common /dpmbw/ ebw(nbw),bwsp(nbw),
     &       bwspa(maxmat,nbw),bwspb(maxmat,nbw),
     &       bwspc(maxmat,nbw),bwspd(maxmat,nbw),idlebw

      character*100 buffer
      integer i,j,ndata
      real*8 dummya(nbw),dummyb(nbw),dummyc(nbw),dummyd(nbw)

      write(*,*) 'xrbw: Reading ', fname
      open(1,file=fname)
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      do 20 j=1,nmat
        read(1,'(a100)') buffer
        read(1,*) ndata
        read(1,'(a100)') buffer
        if (ndata.ne.nbw) then
          write(*,*) 'xrbw:error: Array dim do not match:'
          write(*,'(1x,2(i5,1x))') ndata,nbw
          stop
        endif
c       *** Prepare interpolation:
        do 10 i=1,nbw
          read(1,*) ebw(i),bwsp(i)
          ebw(i) = -1.d0/ebw(i)
 10     continue
        read(1,'(a100)') buffer
        call spline(ebw,bwsp,dummya,dummyb,dummyc,dummyd,
     &              0.d0,0.d0,nbw)
c       *** Loading dummy arrays into multimaterial matrices:
        do 15 i=1,nbw
          bwspa(j,i) = dummya(i)
          bwspb(j,i) = dummyb(i)
          bwspc(j,i) = dummyc(i)
          bwspd(j,i) = dummyd(i)
 15     continue
 20   continue
      close(1)
      idlebw = (nbw-1)/(ebw(nbw)-ebw(1))
      bytes = 8*(maxmat*nbw*4+nbw*2)
      end


      subroutine rqsurf(fname,bytes)
c*******************************************************************
c*    Reads q surface data data from file and sets up interpolation*
c*    matrices                                                     *
c*                                                                 *
c*    Input:                                                       *
c*      fname -> input file name                                   *
c*    Output:                                                      *
c*      bytes -> memory filled up by interpolation arrays          *
c*******************************************************************
      implicit none
      character*40 fname
      integer*4 bytes

      integer maxmat
      parameter (maxmat=5)
      integer nuq,neq
      parameter (nuq=2**7,neq=2**8)
      real q,le0q,idleq,iduq
      common /xdpmq2D/ q(nuq,neq),le0q,idleq,iduq
      integer nmat
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat

      character*100 buffer
      integer i,j,nu,ne
      real*8 e0,e1,u,qmax,effic,qval

      write(*,*) 'rqsurf: Reading ', fname
      open(1,file=fname)
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      read(1,*) nu,ne
      if ((nu.ne.nuq).or.(ne.ne.neq)) then
        write(*,*) 'rqsurf:error: Array dim do not match:'
        write(*,'(1x,4(i5,1x))') nu,nuq,ne,neq
        stop
      endif
      read(1,'(a100)') buffer
      read(1,*) qmax
      read(1,'(a100)') buffer
      effic = 1.d0/qmax
      write(*,*) 'rqsurf: q rejection efficiency(%):'
      write(*,'(1x,1pe8.1)') 100.0*effic
      do 10 i=1,neq
        read(1,*) e1
        if (i.eq.1) e0 = e1
        do 20 j=1,nuq
          read(1,*) u,qval
c         *** Incorporating the efficiency in the q function:
          q(j,i) = qval*effic
 20     continue
 10   continue
      read(1,'(a100)') buffer

      close(1)
      le0q = -1.d0/e0
      idleq = (neq-1)/(-1.d0/e1-le0q)
      iduq = nuq-1.d0
      bytes = 4*nuq*neq
      end


      subroutine xrqsurf(fname,bytes)
c*******************************************************************
c*    Reads q surface data data from file and sets up interpolation*
c*    matrices. Multi-material version.                            *
c*                                                                 *
c*    Input:                                                       *
c*      fname -> input file name                                   *
c*    Output:                                                      *
c*      bytes -> memory filled up by interpolation arrays          *
c*    Comments:                                                    *
c*      -> This routine is identical to rqsurf() except that it    *
c*         loads different data for different materials.           *
c*******************************************************************
      implicit none
      character*40 fname
      integer*4 bytes

      integer maxmat
      parameter (maxmat=5)
      integer nuq,neq
      parameter (nuq=2**7,neq=2**8)
      real q,le0q,idleq,iduq
      common /dpmq2D/ q(maxmat,nuq,neq),le0q,idleq,iduq
      integer nmat
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat

      character*100 buffer
      integer i,j,k,nu,ne
      real*8 e0,e1,u,qmax,effic,qval

      write(*,*) 'xrqsurf: Reading ', fname
      open(1,file=fname)
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      do 30 k=1,nmat
        read(1,'(a100)') buffer
        read(1,*) nu,ne
        if ((nu.ne.nuq).or.(ne.ne.neq)) then
          write(*,*) 'xrqsurf:error: Array dim do not match:'
          write(*,'(1x,4(i5,1x))') nu,nuq,ne,neq
          stop
        endif
        read(1,'(a100)') buffer
        read(1,*) qmax
        read(1,'(a100)') buffer
        effic = 1.d0/qmax
        write(*,*) 'xrqsurf: Mat id & q rejection efficiency(%):'
        write(*,'(1x,i3,1x,1pe8.1)') k,100.0*effic
        do 10 i=1,neq
          read(1,*) e1
          if (i.eq.1) e0 = e1
          do 20 j=1,nuq
            read(1,*) u,qval
c           *** Incorporating the efficiency in the q function:
            q(k,j,i) = qval*effic
 20       continue
 10     continue
        read(1,'(a100)') buffer
 30   continue
      close(1)
      le0q = -1.d0/e0
      idleq = (neq-1)/(-1.d0/e1-le0q)
      iduq = nuq-1.d0
      bytes = 4*(maxmat*nuq*neq)
      end


      subroutine rlammo(fname,bytes)
c*******************************************************************
c*    Reads lambda for Moller interactions from file and sets up   *
c*    interpolation matrices                                       *
c*                                                                 *
c*    Input:                                                       *
c*      fname -> input file name                                   *
c*    Output:                                                      *
c*      bytes -> memory filled up by interpolation arrays          *
c*******************************************************************
      implicit none
      character*40 fname
      integer*4 bytes

      integer nlam
      parameter (nlam=2**8)
      real*8 elam,lamsp,lamspa,lamspb,lamspc,lamspd,idlela
      common /dpmlam/ elam(nlam),lamsp(nlam),
     &  lamspa(nlam),lamspb(nlam),lamspc(nlam),lamspd(nlam),idlela

      character*100 buffer
      integer i,ndata

      write(*,*) 'rlammo: Reading ', fname
      open(1,file=fname)
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      read(1,*) ndata
      if (ndata.ne.nlam) then
        write(*,*) 'rlammo:error: Array dim do not match:'
        write(*,'(1x,2(i5,1x))') ndata,nlam
        stop
      endif
      read(1,'(a100)') buffer
c     *** Prepare interpolation:
      do 10 i=1,nlam
        read(1,*) elam(i),lamsp(i)
        elam(i) = -1.d0/elam(i)
        lamsp(i) = 1.d0/lamsp(i)
 10   continue
      read(1,'(a100)') buffer
      close(1)
      idlela = (nlam-1)/(elam(nlam)-elam(1))
      call spline(elam,lamsp,lamspa,lamspb,lamspc,lamspd,
     &  0.d0,0.d0,nlam)
      bytes = 8*6*nlam
      end


      subroutine rlabre(fname,bytes)
c*******************************************************************
c*    Reads lambda_radiative data from file and sets up            *
c*    interpolation matrices                                       *
c*                                                                 *
c*    Input:                                                       *
c*      fname -> input file name                                   *
c*    Output:                                                      *
c*      bytes -> memory filled up by interpolation arrays          *
c*    Comments:                                                    *
c*      -> The array employed to interpolate the Bremss MFP is     *
c*         prepared to store multi-material data, even though this *
c*         is not needed if the mono-material Sel&Ber DCS option   *
c*         is active.                                              *
c*******************************************************************
      implicit none
      character*40 fname
      integer*4 bytes

      integer nlra,maxmat
      parameter (nlra=2**9,maxmat=5)
      real*8 elra,lrasp,lraspa,lraspb,idlelr
      common /dpmlbr/ elra(nlra),lrasp(nlra),
     &                lraspa(maxmat,nlra),lraspb(maxmat,nlra),
     &                idlelr
      integer nmat
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat

      character*100 buffer
      integer i,ndata,j
      real*8 dummya(nlra),dummyb(nlra)

      write(*,*) 'rlabre: Reading ',fname
      open(1,file=fname)
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      do 20 j=1,nmat
        read(1,'(a100)') buffer
        read(1,*) ndata
        if (ndata.ne.nlra) then
          write(*,*) 'rlabre:error: Array dim do not match:'
          write(*,'(1x,2(i5,1x))') ndata,nlra
          stop
        endif
        read(1,'(a100)') buffer
c       *** Prepare interpolation:
        do 10 i=1,nlra
          read(1,*) elra(i),lrasp(i)
          elra(i) = -1.d0/elra(i)
          lrasp(i) = 1.d0/lrasp(i)
 10     continue
        read(1,'(a100)') buffer
        do 30 i=1,nlra-1
          dummyb(i) = (lrasp(i+1)-lrasp(i))/(elra(i+1)-elra(i))
          dummya(i) = lrasp(i)-elra(i)*dummyb(i)
  30    continue
        do 15 i=1,nlra
          lraspa(j,i) = dummya(i)
          lraspb(j,i) = dummyb(i)
 15     continue
 20   continue
      close(1)
      idlelr = (nlra-1)/(elra(nlra)-elra(1))
      bytes = 8*(2*nlra*maxmat+2*nlra)
      end


      subroutine rbcon(fname,bytes)
c*******************************************************************
c*    Reads brems constants for sampling with PENELOPE routine     *
c*                                                                 *
c*    Input:                                                       *
c*      fname -> input file name                                   *
c*    Output:                                                      *
c*      bytes -> memory filled up by interpolation arrays          *
c*******************************************************************
      implicit none
      character*40 fname
      integer*4 bytes

      integer maxmat
      parameter (maxmat=5)
      integer nmat
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat
      real*8 f0,bcb
      COMMON/XBRE/F0(MAXMAT,3),BCB(MAXMAT)

      character*100 buffer
      integer i,j

      write(*,*) 'rbcon: Reading ', fname
      open(1,file=fname)
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      do 20 j=1,nmat
        read(1,'(a100)') buffer
        read(1,*) (f0(j,i),i=1,3),bcb(j)
 20   continue
      close(1)
      bytes = nmat*4
      end


      subroutine rrstpw(fname,bytes)
c*******************************************************************
c*    Reads StopPower_restr  data from file and sets up            *
c*    interpolation matrices                                       *
c*                                                                 *
c*    Input:                                                       *
c*      fname -> input file name                                   *
c*    Output:                                                      *
c*      bytes -> memory filled up by interpolation arrays          *
c*******************************************************************
      implicit none
      character*40 fname
      integer*4 bytes

      integer maxmat,nst
      parameter (maxmat=5,nst=2**10)
      real*8 est,stsp,stspa,stspb,idlest
      common /dpmstpw/ est(nst),stsp(nst),
     &       stspa(maxmat,nst),stspb(maxmat,nst),idlest
      integer nmat
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat

      character*100 buffer
      integer i,j,ndata
      real*8 dummya(nst),dummyb(nst)

      write(*,*) 'rrstpw: Reading ', fname
      open(1,file=fname)
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      do 20 j=1,nmat
        read(1,'(a100)') buffer
        read(1,*) ndata
        if (ndata.ne.nst) then
          write(*,*) 'rrstpw:error: Array dim do not match:'
          write(*,'(1x,2(i5,1x))') ndata,nst
          stop
        endif
        read(1,'(a100)') buffer
        do 10 i=1,nst
          read(1,*) est(i),stsp(i)
 10     continue
        read(1,'(a100)') buffer
        do 30 i=1,nst-1
          dummyb(i) = (stsp(i+1)-stsp(i))/(est(i+1)-est(i))
          dummya(i) = stsp(i)-est(i)*dummyb(i)
  30    continue
c       *** Loading dummy arrays into multimaterial sp matrices:
        do 15 i=1,nst
          stspa(j,i) = dummya(i)
          stspb(j,i) = dummyb(i)
 15     continue
 20   continue
      close(1)
      idlest = (nst-1)/(est(nst)-est(1))
      bytes = 8*(maxmat*nst*2+nst*2)
      end


      subroutine rlamph(fname,bytes)
c*******************************************************************
c*    Reads photon total inverse mean free path data from file and *
c*    sets up interpolation matrices                               *
c*                                                                 *
c*    Input:                                                       *
c*      fname -> input file name                                   *
c*    Output:                                                      *
c*      bytes -> memory filled up by interpolation arrays          *
c*******************************************************************
      implicit none
      character*40 fname
      integer*4 bytes

      integer maxmat,nlaph
      parameter (maxmat=5,nlaph=2**11)
      real*8 elaph,lamph,lampha,lamphb,lamphc,lamphd,idleph
      common /dpmlph/ elaph(nlaph),lamph(nlaph),
     &       lampha(maxmat,nlaph),lamphb(maxmat,nlaph),
     &       lamphc(maxmat,nlaph),lamphd(maxmat,nlaph),idleph
      integer nmat
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat

      character*100 buffer
      integer i,j,ndata
      real*8 dummya(nlaph),dummyb(nlaph),dummyc(nlaph),dummyd(nlaph)

      write(*,*) 'rlamph: Reading ', fname
      open(1,file=fname)
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      do 20 j=1,nmat
        read(1,'(a100)') buffer
        read(1,*) ndata
        if (ndata.ne.nlaph) then
          write(*,*) 'rlamph:error: Array dim do not match:'
          write(*,'(1x,2(i5,1x))') ndata,nlaph
          stop
        endif
        read(1,'(a100)') buffer
c       *** Preparing interpolation:
        do 10 i=1,nlaph
          read(1,*) elaph(i),lamph(i)
 10     continue
        read(1,'(a100)') buffer
        call spline(elaph,lamph,dummya,dummyb,dummyc,dummyd,
     &              0.d0,0.d0,nlaph)
c       *** Loading dummy arrays into multimaterial sp matrices:
        do 15 i=1,nlaph
          lampha(j,i) = dummya(i)
          lamphb(j,i) = dummyb(i)
          lamphc(j,i) = dummyc(i)
          lamphd(j,i) = dummyd(i)
 15     continue
 20   continue
      close(1)
      idleph = (nlaph-1)/(elaph(nlaph)-elaph(1))
      bytes = 8*(maxmat*nlaph*4+nlaph*2)
      end


      subroutine rcompt(fname,bytes)
c*******************************************************************
c*    Reads Compton inverse mean free path data from file and sets *
c*    up interpolation matrices                                    *
c*                                                                 *
c*    Input:                                                       *
c*      fname -> input file name                                   *
c*    Output:                                                      *
c*      bytes -> memory filled up by interpolation arrays          *
c*******************************************************************
      implicit none
      character*40 fname
      integer*4 bytes

      integer maxmat,ncmpt
      parameter (maxmat=5,ncmpt=2**8)
      real*8 ecmpt,compt,compta,comptb,comptc,comptd,idlecp
      common /dpmcmp/ ecmpt(ncmpt),compt(ncmpt),
     &       compta(maxmat,ncmpt),comptb(maxmat,ncmpt),
     &       comptc(maxmat,ncmpt),comptd(maxmat,ncmpt),idlecp
      integer nmat
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat

      character*100 buffer
      integer i,j,ndata
      real*8 dummya(ncmpt),dummyb(ncmpt),dummyc(ncmpt),dummyd(ncmpt)

      write(*,*) 'rcompt: Reading ', fname
      open(1,file=fname)
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      do 20 j=1,nmat
        read(1,'(a100)') buffer
        read(1,*) ndata
        if (ndata.ne.ncmpt) then
          write(*,*) 'rcompt:error: Array dim do not match:'
          write(*,'(1x,2(i5,1x))') ndata,ncmpt
          stop
        endif
        read(1,'(a100)') buffer
c       *** Preparing interpolation:
        do 10 i=1,ncmpt
          read(1,*) ecmpt(i),compt(i)
 10     continue
        read(1,'(a100)') buffer
        call spline(ecmpt,compt,dummya,dummyb,dummyc,dummyd,
     &              0.d0,0.d0,ncmpt)
c       *** Loading dummy arrays into multimaterial sp matrices:
        do 15 i=1,ncmpt
          compta(j,i) = dummya(i)
          comptb(j,i) = dummyb(i)
          comptc(j,i) = dummyc(i)
          comptd(j,i) = dummyd(i)
 15     continue
 20   continue
      close(1)
      idlecp = (ncmpt-1)/(ecmpt(ncmpt)-ecmpt(1))
      bytes = 8*(maxmat*ncmpt*4+ncmpt*2)
      end


      subroutine rphote(fname,bytes)
c*******************************************************************
c*    Reads photoelectric inverse mean free path data from file and*
c*    sets up interpolation matrices                               *
c*                                                                 *
c*    Input:                                                       *
c*      fname -> input file name                                   *
c*    Output:                                                      *
c*      bytes -> memory filled up by interpolation arrays          *
c*******************************************************************
      implicit none
      character*40 fname
      integer*4 bytes

      integer maxmat,nphte
      parameter (maxmat=5,nphte=2**12)
      real*8 ephte,phote,photea,photeb,photec,photed,idlepe
      common /dpmpte/ ephte(nphte),phote(nphte),
     &       photea(maxmat,nphte),photeb(maxmat,nphte),
     &       photec(maxmat,nphte),photed(maxmat,nphte),idlepe
      integer nmat
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat

      character*100 buffer
      integer i,j,ndata
      real*8 dummya(nphte),dummyb(nphte),dummyc(nphte),dummyd(nphte)

      write(*,*) 'rphote: Reading ', fname
      open(1,file=fname)
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer
      do 20 j=1,nmat
        read(1,'(a100)') buffer
        read(1,*) ndata
        if (ndata.ne.nphte) then
          write(*,*) 'rphote:error: Array dim do not match:'
          write(*,'(1x,2(i5,1x))') ndata,nphte
          stop
        endif
        read(1,'(a100)') buffer
c       *** Preparing interpolation
        do 10 i=1,nphte
          read(1,*) ephte(i),phote(i)
 10     continue
        read(1,'(a100)') buffer
        call spline(ephte,phote,dummya,dummyb,dummyc,dummyd,
     &              0.d0,0.d0,nphte)
c       *** Loading dummy arrays into multimaterial sp matrices:
        do 15 i=1,nphte
          photea(j,i) = dummya(i)
          photeb(j,i) = dummyb(i)
          photec(j,i) = dummyc(i)
          photed(j,i) = dummyd(i)
 15     continue
 20   continue
      close(1)
      idlepe = (nphte-1)/(ephte(nphte)-ephte(1))
      bytes = 8*(maxmat*nphte*4+nphte*2)
      end


      subroutine rpairp(fname,bytes)
c*******************************************************************
c*    Reads inverse pair production mean free path data from file  *
c*    and sets up interpolation matrices                           *
c*                                                                 *
c*    Input:                                                       *
c*      fname -> input file name                                   *
c*    Output:                                                      *
c*      bytes -> memory filled up by interpolation arrays          *
c*******************************************************************
      implicit none
      character*40 fname
      integer*4 bytes

      integer maxmat,npair
      parameter (maxmat=5,npair=2**10)
      real*8 epair,pairp,pairpa,pairpb,pairpc,pairpd,idlepp
      common /dpmpap/ epair(npair),pairp(npair),
     &       pairpa(maxmat,npair),pairpb(maxmat,npair),
     &       pairpc(maxmat,npair),pairpd(maxmat,npair),idlepp
      integer nmat
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat

      character*100 buffer
      integer i,j,ndata
      real*8 dummya(npair),dummyb(npair),dummyc(npair),dummyd(npair)

      write(*,*) 'rpairp: Reading ', fname
      open(1,file=fname)
      read(1,'(a100)') buffer
      read(1,'(a100)') buffer

      do 20 j=1,nmat
        read(1,'(a100)') buffer
        read(1,*) ndata
        if (ndata.ne.npair) then
          write(*,*) 'rpairp:error: Array dim do not match:'
          write(*,'(1x,2(i5,1x))') ndata,npair
          stop
        endif
        read(1,'(a100)') buffer
c       *** Preparing interpolation:
        do 10 i=1,npair
          read(1,*) epair(i),pairp(i)
 10     continue
        read(1,'(a100)') buffer
        call spline(epair,pairp,dummya,dummyb,dummyc,dummyd,
     &              0.d0,0.d0,npair)
c       *** Loading dummy arrays into multimaterial sp matrices:
        do 15 i=1,npair
          pairpa(j,i) = dummya(i)
          pairpb(j,i) = dummyb(i)
          pairpc(j,i) = dummyc(i)
          pairpd(j,i) = dummyd(i)
 15     continue
 20   continue
      close(1)
      idlepp = (npair-1)/(epair(npair)-epair(1))
      bytes = 8*(maxmat*npair*4+npair*2)
      end


      subroutine rvoxg(fname)
c*******************************************************************
c*    Reads voxel geometry from an input file                      *
c*                                                                 *
c*    Input:                                                       *
c*      fname -> input file name                                   *
c*    Comments:                                                    *
c*      -> rmater must be called first to set nmat.                *
c*******************************************************************
      implicit none
      character*100 fname

      integer*4 nxvox,nyvox,nzvox,nxyz
      parameter (nxvox=2**6,nyvox=2**6,nzvox=150)
      parameter (nxyz=nxvox*nyvox*nzvox)
      integer mat
      real*8 dens,dx,dy,dz,xmid,ymid,zmid
      integer Unxvox,Unyvox,Unzvox
      common /dpmvox/ dx,dy,dz,xmid,ymid,zmid
     +          ,dens(nxyz),mat(nxyz),Unxvox,Unyvox,Unzvox
      integer maxmat
      parameter (maxmat=5)
      integer nmat
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat

      character*100 buffer
      integer*4 i,j,k,absvox,getabs

      write(*,*) ' '
      write(*,*) ' '
      write(*,*) 'rvoxg: Reading ', fname
      write(*,*) '       information from this file follows:'
      open(1,file=fname)
      read(1,'(a100)') buffer
      write(*,*) buffer
      read(1,'(a100)') buffer
      write(*,*) buffer
      read(1,'(a100)') buffer
      write(*,*) buffer
      read(1,*) Unxvox,Unyvox,Unzvox
      write(*,'(1x,3(i10,1x))') Unxvox,Unyvox,Unzvox
      if (Unxvox.gt.nxvox.or.Unyvox.gt.nyvox.or.Unzvox.gt.nzvox) then
        write(*,*) 'rvoxg:error: Too many voxels! Increase nxyz.'
        stop
      endif
      read(1,'(a100)') buffer
      write(*,*) buffer
      read(1,*) dx,dy,dz
      write(*,'(1x,3(1pe10.3,1x))') dx,dy,dz

      xmid = dx*Unxvox/2.0d0
      ymid = dy*Unyvox/2.0d0
      zmid = dz*Unzvox/2.0d0
c     *** Load material and density for each voxel:
      read(1,'(a100)') buffer
      do 10 i=1,Unxvox
        do 20 j=1,Unyvox
          do 30 k=1,Unzvox
            absvox = getabs(i,j,k)
            read(1,*) mat(absvox),dens(absvox)
            if (mat(absvox).gt.nmat) then
              write(*,*) 'rvoxg:error: Mat# too large:'
              write(*,'(1x,i6)') mat(absvox)
              stop
            endif
 30       continue
 20     continue
 10   continue
      close(1)
      end


      subroutine iniwck(eminph,emax,bytes)
c*******************************************************************
c*    Finds information used to transport photons with the Woodcock*
c*    technique                                                    *
c*                                                                 *
c*    Input:                                                       *
c*      eminph -> minimum photon energy in data files (eV)         *
c*      emax -> maximum photon energy in data files (eV)           *
c*    Output                                                       *
c*      bytes -> space allocated for arrays                        *
c*    Comments:                                                    *
c*      -> common /dpmsrc/ must be loaded previously               *
c*      -> rlamph() must be called previously                      *
c*      -> rvoxg() must be called first                            *
c*      -> emax reduced to avoid reaching the end of interpol table*
c*******************************************************************
      implicit none
      integer*4 bytes
      real*8 eminph,emax

      integer maxmat,nlaph
      parameter (maxmat=5,nlaph=2**11)
      integer*4 nxvox,nyvox,nzvox,nxyz
      parameter (nxvox=2**6,nyvox=2**6,nzvox=150)
      parameter (nxyz=nxvox*nyvox*nzvox)
      integer mat
      real*8 dens,dx,dy,dz,xmid,ymid,zmid
      integer Unxvox,Unyvox,Unzvox
      common /dpmvox/ dx,dy,dz,xmid,ymid,zmid
     +          ,dens(nxyz),mat(nxyz),Unxvox,Unyvox,Unzvox
      integer nwck
      parameter (nwck=2**12)
      real*8 a0wck,a1wck,idlewk,wcke0
      common /dpmwck/ a0wck(nwck),a1wck(nwck),idlewk,wcke0
      integer intype
      real*8 esrc,eabs,eabsph,param
      common /dpmsrc/ esrc,eabs,eabsph,param,intype
      integer nmat
      real*8 matden,zmass,z2mass,wcc,wcb
      common /dpmmat/ matden(maxmat),zmass(maxmat),z2mass(maxmat),
     &                wcc,wcb,nmat

      integer i,j,granul,ke
      parameter (granul=3)
      integer*4 vox
      real*8 maxden(maxmat),de,e,ymax,ymin,ycanbe
      real*8 itphip,ylast,me,avedif,maxdif,gle,fip,freal,dif,eps
      parameter (eps=1.d-10)

c     MPI相关参数
      integer myrank,pro_num,ierr
      common /mpi_para/ myrank,pro_num,ierr
      
      if (myrank.eq.0) then
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) 'iniwck: Started.'
      endif
c     *** Find the largest density for each present material:
      do 10 i=1,maxmat
        maxden(i) = 0.d0
 10   continue
      do 20 vox=1,Unxvox*Unyvox*Unzvox
        if (dens(vox).gt.maxden(mat(vox))) maxden(mat(vox))=dens(vox)
 20   continue

c     *** Prepare linear interpolation:
      wcke0 = eminph
      de = (emax*(1.d0-eps)-wcke0)/nwck
      idlewk = 1.d0/de
      do 30 i=1,nwck+1
        e = wcke0+de*(i-1)
        ymax = 0.d0
        do 40 j=1,nmat
          ycanbe = itphip(j,e)*maxden(j)
          if (ycanbe.gt.ymax) ymax = ycanbe
 40     continue
        ymin = 1.d0/ymax
        if (i.eq.1) then
          ylast = ymin
          goto 30
        endif
        a1wck(i-1) = (ymin-ylast)*idlewk
        a0wck(i-1) = ylast-(e-de)*a1wck(i-1)
        ylast = ymin
 30   continue
      bytes = 8*2*nwck

c     *** Check accuracy inside mesh:
      me = wcke0
      avedif = 0.d0
      maxdif = 0.d0
      do 50 i=1,nwck
        e = wcke0+de*(i-1)
        do 60 ke=1,granul
          gle = e+ke*de/(granul+1)
          fip = a0wck(i)+a1wck(i)*gle
          ymax = 0.d0
          do 70 j=1,nmat
            ycanbe = itphip(j,gle)*maxden(j)
            if (ycanbe.gt.ymax) ymax = ycanbe
 70       continue
          freal = 1.d0/ymax
          dif = dabs(1.d0-fip/freal)
          if (maxdif.lt.dif) then
            maxdif = dif
            me = gle
          endif
          avedif = avedif+dif
 60     continue
 50   continue
      avedif = 100.d0*avedif/(nwck*granul)
      maxdif = 100.d0*maxdif

      if (myrank.eq.0) then
      write(*,*)
     &'iniwck: ndata,granularity,mean_dif(%),max_dif(%) at(eV):'
      write(*,'(1x,i5,1x,i3,1x,3(1pe10.3,1x))')
     &  nwck,granul,avedif,maxdif,me
      endif
      end


      subroutine inisub(refden)
c*******************************************************************
c*    Initializes data used when transporting electrons below      *
c*    Eabs and reads StopPower from file to set up interpolation   *
c*    matrices                                                     *
c*                                                                 *
c*    Input:                                                       *
c*      refden -> reference material density (g/cm^3)              *
c*    Comments:                                                    *
c*      -> must be called after reading common /dpmsrc/            *
c*******************************************************************
      implicit none
      real*8 refden

      integer intype
      real*8 esrc,eabs,eabsph,param
      common /dpmsrc/ esrc,eabs,eabsph,param,intype
      real*8 subden,subfac,substp
      common /dpmsub/ subden,subfac,substp
c     MPI相关参数
      integer myrank,pro_num,ierr
      common /mpi_para/ myrank,pro_num,ierr
      
c     *** Density and StopPower for SubEabs transport are set quite arbitrarily:
      subden = refden/10.0d0
      substp = 2.0d6
      subfac = subden/eabs
      if (myrank.eq.0) then
      write(*,*)
     &  'inisub: Minimum density for SubEabs motion (g/cm^3):'
      write(*,'(1x,1pe10.3)') subden
      write(*,*)
     &  'inisub: StopPower for SubEabs motion (eV*cm^2/g):'
      write(*,'(1x,1pe10.3)') substp
      endif

      end


      subroutine putann
c*******************************************************************
c*    Creates secondaries for annhilation photons,                 *
c*    stores states in the secondary stack                         *
c*******************************************************************
      implicit none
      integer ptype
      integer*4 xvox,yvox,zvox,absvox,nhist
      real*8 energy,vx,vy,vz,x,y,z
      common /dpmpart/ energy,vx,vy,vz,x,y,z,
     &                 xvox,yvox,zvox,absvox,nhist,ptype
      integer sptype
      integer*4 maxns,nsec,sxvox,syvox,szvox,sabsvx
      real*8 swght,senerg,svx,svy,svz,sx,sy,sz
      parameter (maxns=2**7)
      common /dpmstck/ swght(maxns),senerg(maxns),
     & svx(maxns),svy(maxns),svz(maxns),sx(maxns),sy(maxns),sz(maxns),
     & sptype(maxns),sxvox(maxns),syvox(maxns),szvox(maxns),
     & sabsvx(maxns),nsec
      integer intype
      real*8 esrc,eabs,eabsph,param
      common /dpmsrc/ esrc,eabs,eabsph,param,intype

      integer j
      real*8 rng,twopi,pi,mc2,sinthe,phi
      parameter (pi=3.1415926535897932d0,twopi=2.d0*pi)
      parameter (mc2=510.9991d3)

      vz = 2.d0*rng()-1.d0
      sinthe = dsqrt(1.d0-vz*vz)
      phi = twopi*rng()
      vy = sinthe*dsin(phi)
      vx = sinthe*dcos(phi)
      do 10 j=1,2
        nsec = nsec+1
        if (nsec.gt.maxns) then
          write(*,*) 'putann:error: stack is full, enlarge maxns.'
          stop
        endif
        sptype(nsec) = 0
        senerg(nsec) = mc2
        svx(nsec) = vx
        svy(nsec) = vy
        svz(nsec) = vz
        sx(nsec) = x
        sy(nsec) = y
        sz(nsec) = z
        sxvox(nsec) = xvox
        syvox(nsec) = yvox
        szvox(nsec) = zvox
        sabsvx(nsec) = absvox
   10 continue
      svx(nsec) = -svx(nsec)
      svy(nsec) = -svy(nsec)
      svz(nsec) = -svz(nsec)
      end


      real*8 function SBREMS(E,M)
c*******************************************************************
c*    New routine to replace sambre in dpm                         *
c*    for more accurate sampling of Bremsstrahlung energies.       *
c*    Canibalized from PENELOPE.                                   *
c*******************************************************************
      integer m
      real*8 e
C
C     RANDOM SAMPLING OF HARD BREMSSTRAHLUNG EMISSION.
C
      integer maxmat
      parameter (maxmat=5)
      real*8 pi,twopi,rev
      PARAMETER (PI=3.1415926535897932D0, TWOPI=2.0D0*PI)
      PARAMETER (REV=5.1099906D5)
C
      real*8 wcion,wcbre
      common /cutoff/ wcion,wcbre
C
      real*8 f0,bcb
      COMMON/XBRE/F0(MAXMAT,3),BCB(MAXMAT)
C
      real*8 rng
C
      integer ic
      real*8 toten, gam, em, ed, ec, ec1, ec2, emecs, xx, f2le
      real*8 f00, f1, f2, f1c, f2c, f1d, f2d, ed1, bec, bed
      real*8 eps, eps1, pa, pa0, pa1, pa2, pb, pb2
C
      TOTEN=E+REV
      GAM=TOTEN/REV
      PA=BCB(M)*GAM
      EM=E/TOTEN
      IF(EM.GT.9.9999D-1) EM=9.9999D-1
C  ****  REJECTION FUNCTIONS.
      ED=(E-5.0D0*REV)/TOTEN
      EC=wcbre/TOTEN
      EC2=EC*EC
      EMECS=EM*EM-EC2
C  ****  LOW ENERGY CORRECTION.
      XX=E/REV
      F2LE=(4.650D0-F0(M,3)*(6.005D0-F0(M,3)*2.946D0))/XX
     1    -(3.242D1-F0(M,3)*(6.708D1+F0(M,3)*3.906D0))/(1+XX)**2
     2    +(2.033D1+F0(M,3)*(2.338D1-F0(M,3)*7.742D1))/(1+XX)**3
C
      IF(EC.LT.ED) THEN
        IC=2
      ELSE
        IC=1
      ENDIF
      F00=F0(M,IC)+F2LE
      EC1=1.0D0-EC
      BEC=EC/(PA*EC1)
      CALL SCHIFF(BEC,F1,F2)
      F1C=F1+F00
      F2C=(F2+F00)*EC1
      IF(EC.LT.ED) THEN
        F00=F0(M,1)+F2LE
        ED1=1.0D0-ED
        BED=ED/(PA*ED1)
        CALL SCHIFF(BED,F1,F2)
        F1D=F1+F00
        F2D=(F2+F00)*EC1
        IF(F1D.GT.F1C) F1C=F1D
        IF(F2D.GT.F2C) F2C=F2D
      ENDIF
      PA1=EMECS*F1C
      PA2=2.66666666666666D0*DLOG(EM/EC)*F2C
C
    1 CONTINUE
      IF(rng()*(PA1+PA2).GT.PA1) GO TO 2
      EPS=DSQRT(EC2+rng()*EMECS)
      PB=EPS/(PA*(1.0D0-EPS))
      F1=2.0D0-2.0D0*DLOG(1.0D0+PB*PB)
      IF(PB.LT.1.0D-8) THEN
        F1=F1-PB*TWOPI
      ELSE
        F1=F1-4.0D0*PB*DATAN2(1.0D0,PB)
      ENDIF
      IF(EPS.LT.ED) THEN
        F00=F0(M,2)+F2LE
      ELSE
        F00=F0(M,1)+F2LE
      ENDIF
      IF(rng()*F1C.GT.F1+F00) GO TO 1
      GO TO 3
    2 CONTINUE
      EPS=EC*(EM/EC)**rng()
      EPS1=1.0D0-EPS
      PB=EPS/(PA*EPS1)
      PB2=PB*PB
      F1=2.0D0-2.0D0*DLOG(1.0D0+PB2)
      F2=F1-6.666666666666666D-1
      IF(PB.LT.1.0D-8) THEN
        F1=F1-PB*TWOPI
      ELSE
        PA0=4.0D0*PB*DATAN2(1.0D0,PB)
        F1=F1-PA0
        F2=F2+2.0D0*PB2*(4.0D0-PA0-3.0D0*DLOG((1.0D0+PB2)/PB2))
      ENDIF
      F2=0.5D0*(3.0D0*F1-F2)
      IF(EPS.LT.ED) THEN
        F00=F0(M,2)+F2LE
      ELSE
        F00=F0(M,1)+F2LE
      ENDIF
      IF(rng()*F2C.GT.EPS1*(F2+F00)) GO TO 1
    3 sbrems=EPS*TOTEN
      RETURN
      END


      SUBROUTINE SCHIFF(B,F1,F2)
c*******************************************************************
c*    Canibalized from PENELOPE too.                               *
c*******************************************************************
C
C     SCREENING FUNCTIONS F1(B) AND F2(B) IN THE BETHE-HEITLER
C  DIFFERENTIAL CROSS SECTION FOR BREMSSTRAHLUNG EMISSION.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit integer (i-n)
      PARAMETER (PI=3.1415926535897932D0, TWOPI=2.0D0*PI)
      B2=B*B
      F1=2.0D0-2.0D0*DLOG(1.0D0+B2)
      F2=F1-6.666666666666666D-1
      IF(B.LT.1.0D-10) THEN
        F1=F1-TWOPI*B
      ELSE
        A0=4.0D0*B*DATAN2(1.0D0,B)
        F1=F1-A0
        F2=F2+2.0D0*B2*(4.0D0-A0-3.0D0*DLOG((1.0D0+B2)/B2))
      ENDIF
      F2=0.5D0*(3.0D0*F1-F2)
      RETURN
      END


c*** end of file *****************************************************
